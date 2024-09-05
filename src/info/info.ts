/**
 * @packageDocumentation
 * @module info
 */

import assert from 'assert';

import { Parameter, Property } from '../dataset';
import { DisplayTarget, EnvironmentIndexer, Indexes } from '../indexer';
import { binarySearch, getElement, logger, sendWarning } from '../utils';

import * as plotlyStyles from '../map/plotly/plotly-styles';

import Collapse from '../collapse';
import { Slider } from './slider';
import { Table } from './table';

import INFO_SVG from '../static/info.svg';
import * as styles from '../styles';

export function filter<T extends Record<string, Property>>(
    obj: T,
    predicate: (o: Property) => boolean
): Record<string, Property> {
    const result: Record<string, Property> = {};
    for (const key in obj) {
        if (Object.hasOwnProperty.call(obj, key) && predicate(obj[key])) {
            result[key] = obj[key];
        }
    }
    return result;
}

/**
 * Information associated with the current structure or atom
 */
interface Info {
    /** Manual id setting with input[type=number] */
    number: HTMLInputElement;
    /** Slider to select an environment / play the trajectory */
    slider: Slider;
    /** Property table */
    table: Table;
}

/**
 * The {@link EnvironmentInfo} class displays information about structure or
 * environment, using a slider for selection and an hidden by default table
 * displaying all properties.
 */
export class EnvironmentInfo {
    /** Callback used when the user changes one of the sliders value */
    public onchange: (indexes: Indexes) => void = () => {};
    /** delay in ms between "frames" during playback over structures/atoms, default 700 ms */
    public playbackDelay: number = 700;

    private _shadow: ShadowRoot;
    private _root: HTMLElement;
    private _atom?: Info;
    private _structure: Info;
    private _indexer: EnvironmentIndexer;
    private _properties: { [name: string]: Property };
    private _parameters?: { [name: string]: Parameter };
    private _target: 
  ;

    /**
     * Create a new {@link EnvironmentInfo} inside the DOM element with given `id`
     * @param element    HTML element or string 'id' of the element where
     *                   the sliders and tables should live
     * @param properties properties to be displayed
     * @param indexer    {@link EnvironmentIndexer} used to translate indexes from
     *                   environments index to structure/atom indexes
     * @param target     Display target, either atom or structure
     * @param parameters used to describe multidimensional properties
     */
    constructor(
        element: string | HTMLElement,
        properties: { [name: string]: Property },
        indexer: EnvironmentIndexer,
        target: DisplayTarget,
        parameters?: { [name: string]: Parameter }
    ) {
        // Create a host element to attach the shadow DOM
        const containerElement = getElement(element);
        const hostElement = document.createElement('div');
        containerElement.appendChild(hostElement);
        this._shadow = hostElement.attachShadow({ mode: 'open' });

        // Add styles to shadow
        this._shadow.adoptedStyleSheets = [
            styles.bootstrap,
            styles.chemiscope,
            plotlyStyles.globalStyleSheet,
        ];

        this._root = document.createElement('div');
        this._shadow.appendChild(this._root);

        this._indexer = indexer;
        this._target = target;

        // Save properties and parameters to re/create structure and atom elements once target changed
        this._properties = properties;
        this._parameters = parameters;

        // Initialize _structure with a placeholder
        this._structure = {
            number: document.createElement('input'),
            slider: new Slider(document.createElement('div'), 'structure'),
            table: new Table(
                document.createElement('div'),
                'structure',
                'structure',
                properties,
                parameters
            ),
        };

        // Create structure and atom html elements
        this._renderTargetPart();
    }

    /**
     * Change widget display target and adapt the element to the new target
     * @param target display target
     */
    public switchTarget(target: DisplayTarget) {
        // Update widget target
        this._target = target;

        // Update element related to the display target
        this._renderTargetPart();
    }

    /** Show properties for the given `indexes`, and update the sliders values */
    public show(indexes: Indexes): void {
        const previousStructure = this._indexes().structure;
        this._structure.number.value = `${indexes.structure + 1}`;
        this._structure.slider.update(indexes.structure);
        this._structure.table.show(indexes);

        if (indexes.structure !== previousStructure && this._atom !== undefined) {
            const activeAtoms = this._indexer.activeAtoms(indexes.structure);
            this._atom.number.value = `${activeAtoms[0] + 1}`;
            this._atom.slider.reset(activeAtoms);
        }

        if (indexes.atom !== undefined) {
            if (this._atom === undefined) {
                if (indexes.atom !== 0) {
                    throw Error(
                        'Invalid state: got an atomic number to update, but I am displaying only structures'
                    );
                }
            } else {
                this._atom.number.value = `${indexes.atom + 1}`;
                this._atom.slider.update(indexes.atom);
                this._atom.table.show(indexes);
            }
        }
    }

    /**
     * Remove all HTML added by this {@link EnvironmentInfo} in the current document
     */
    public remove(): void {
        this._shadow.host.remove();
    }

    /** Changes the elements based on the target */
    private _renderTargetPart(): void {
        // Construct structure / atom buttons
        this._root.innerHTML = this._getMainHTMLStructure();

        // Filter properties for structure
        const structureProperties = filter(this._properties, (p) => p.target === 'structure');
        this._structure = this._createStructure(structureProperties, this._parameters);

        // Initialize central atom if in atom target
        if (this._target === 'atom') {
            const atomProperties = filter(this._properties, (p) => p.target === 'atom');
            this._atom = this._createAtom(atomProperties, this._parameters);
        }

        // Reset atom if in structure target
        else {
            this._atom = undefined;
        }

        // Initialize the collapse components from their 'data-bs-*' attributes.
        Collapse.initialize(this._root);
    }

    /** Create base html for the structure / atom pannels */
    private _getMainHTMLStructure(): string {
        const atomButton =
            this._target === 'atom'
                ? `
            <div class='btn btn-sm chsp-info-atom-btn' data-bs-toggle='collapse' data-bs-target='#atom' aria-expanded='false' aria-controls='atom'>
                <div class="chsp-info-btns-svg">${INFO_SVG}</div>
                atom <input class='chsp-info-number' type=number value=1 min=1></input>
            </div>`
                : '<div></div>';
        return `
            <div class='accordion chsp-info-tables' id='info-tables'></div>
            <div class='chsp-info-btns'>
                <div class='btn btn-sm chsp-info-structure-btn' data-bs-toggle='collapse' data-bs-target='#structure' aria-expanded='false' aria-controls='structure'>
                    <div class="chsp-info-btns-svg">${INFO_SVG}</div>
                    structure <input class='chsp-info-number' type=number value=1 min=1></input>
                </div>
                ${atomButton}
            </div>`;
    }

    /** Create the structure slider and table */
    private _createStructure(
        properties: { [name: string]: Property },
        parameters?: { [name: string]: Parameter }
    ): Info {
        // Slider
        const slider = new Slider(this._root, 'structure');
        slider.reset(this._indexer.activeStructures());
        slider.startPlayback = (advance) => {
            // Create a delay before playback step run
            setTimeout(() => {
                if (advance()) {
                    // Find the next valid structure index
                    const structure = this._findNextValidIndex(
                        // Get current structure index
                        this._indexes().structure,

                        // Structures number
                        this._indexer.structuresCount(),

                        // Function to get the indexer for a given structure
                        (index) => this._indexer.fromStructure(index, this._target)
                    );

                    // Valid structure was found
                    if (structure !== undefined) {
                        // Update the display with the details of the found structure
                        this.show(this._indexer.fromStructure(structure, this._target)!);
                        this.onchange(this._indexes());

                        // Recursively call startPlayback to continue playback
                        slider.startPlayback(advance);
                    }
                }
            }, this.playbackDelay);
        };
        slider.onchange = () => {
            // Get the current structure from slider
            const structure = this._structure.slider.value();

            // Atom is preset -> update atom-related ui elements
            if (this._atom !== undefined) {
                const activeAtoms = this._indexer.activeAtoms(structure);
                if (activeAtoms.length === 0) {
                    logger.warn(
                        `Cannot change to structure ${
                            structure + 1
                        }, which does not contain any active atoms`
                    );
                    return;
                }

                // Update atom number and slider
                this._atom.number.value = `${activeAtoms[0] + 1}`;
                this._atom.number.max = `${activeAtoms.length}`;
                this._atom.slider.reset(activeAtoms);
            }

            // Get and validate indexes for the current structure
            const indexes = this._indexes();
            assert(indexes !== undefined);

            // Update structure and atom tables
            this._structure.table.show(indexes);
            this._structure.number.value = `${indexes.structure + 1}`;
            if (this._atom !== undefined) {
                this._atom.table.show(indexes);
            }

            // Trigger callback with updated indexes
            this.onchange(indexes);
        };

        // Table
        const tableRoot = this._root.children[0] as HTMLElement;
        assert(tableRoot.tagName.toLowerCase() === 'div');
        const table = new Table(tableRoot, 'structure', 'structure', properties, parameters);

        // Number input
        const number = this._createNumberInput(
            '.chsp-info-structure-btn .chsp-info-number',

            // Max value
            this._indexer.structuresCount().toString(),

            // On change callback
            (value) => {
                // Invalid value -> reset
                if (isNaN(value) || value < 0 || value >= parseInt(number.max, 10)) {
                    number.value = `${this._structure.slider.value() + 1}`;
                }

                // Not an active structure -> reset
                else if (binarySearch(this._indexer.activeStructures(), value) === -1) {
                    number.value = `${this._structure.slider.value() + 1}`;
                }

                // Proceed with the new value
                else {
                    this._structure.slider.update(value);

                    // If atom is present, update its related properties
                    if (this._atom !== undefined) {
                        const activeAtoms = this._indexer.activeAtoms(value);
                        this._atom.number.value = `${activeAtoms[0] + 1}`;
                        this._atom.number.max = `${activeAtoms.length}`;
                        this._atom.slider.reset(activeAtoms);
                    }

                    // Show the updated structure table
                    const indexes = this._indexes();
                    this._structure.table.show(indexes);
                    if (this._atom !== undefined) {
                        this._atom.table.show(indexes);
                    }

                    // Trigger callback with the updated indexes
                    this.onchange(indexes);
                }
            }
        );

        return { number, slider, table };
    }

    /** Create the atom slider and table */
    private _createAtom(
        properties: { [name: string]: Property },
        parameters?: { [name: string]: Parameter }
    ) {
        // Slider
        const slider = new Slider(this._root, 'atom');
        slider.reset(this._indexer.activeAtoms(this._structure.slider.value()));
        slider.startPlayback = (advance) => {
            // Create a delay before playback step run
            setTimeout(() => {
                // Check if the advance callback function allows playback to continue
                if (advance()) {
                    // Get the current indexes
                    const current = this._indexes();
                    assert(current.atom !== undefined);

                    // Get next atom index
                    const atom = this._findNextValidIndex(
                        // Start from the current atom index
                        current.atom,

                        // Atoms number
                        this._indexer.atomsCount(current.structure),

                        // Function to get the indexer for a given atom
                        (index) =>
                            this._indexer.fromStructureAtom(this._target, current.structure, index)
                    );

                    // Valid atom was found
                    if (atom !== undefined) {
                        // Update the display with the details of the found atom
                        this.show(
                            this._indexer.fromStructureAtom(this._target, current.structure, atom)!
                        );
                        this.onchange(this._indexes());

                        // Recursively call startPlayback to continue playback
                        slider.startPlayback(advance);
                    }
                }
            }, this.playbackDelay);
        };
        slider.onchange = () => {
            assert(this._atom !== undefined);

            // Get the current indexes
            const indexes = this._indexes();
            if (indexes === undefined) {
                const structure = this._structure.slider.value();
                const atom = this._atom.slider.value();
                logger.warn(
                    `Environment for atom ${atom} in structure ${structure} is not part of this dataset`
                );
                return;
            }

            // Update the atom table and number input with the new indexes
            assert(indexes.atom !== undefined);
            this._atom.table.show(indexes);
            this._atom.number.value = `${indexes.atom + 1}`;

            // Trigger callback with the updated indexes
            this.onchange(indexes);
        };

        // Table
        const tableRoot = this._root.children[0] as HTMLElement;
        assert(tableRoot.tagName.toLowerCase() === 'div');
        const table = new Table(tableRoot, 'atom', 'atom', properties, parameters);

        // Number input
        const number = this._createNumberInput(
            '.chsp-info-atom-btn .chsp-info-number',

            // Max number
            this._indexer.atomsCount(this._structure.slider.value()).toString(),

            // On change callback
            (value) => {
                assert(this._atom !== undefined);
                const activeAtoms = this._indexer.activeAtoms(this._structure.slider.value());

                // Invalid value -> reset
                if (isNaN(value) || value < 0 || value >= parseInt(number.max, 10)) {
                    // reset to the current slider value if we got an invalid value
                    number.value = `${this._atom.slider.value() + 1}`;
                }

                // Not an active structure -> reset
                else if (binarySearch(activeAtoms, value) === -1) {
                    number.value = `${this._atom.slider.value() + 1}`;
                }

                // Proceed with the new value
                else {
                    this._atom.slider.update(value);

                    // Show the updated structure table
                    const indexes = this._indexes();
                    this._atom.table.show(indexes);

                    // Trigger callback with the updated indexes
                    this.onchange(indexes);
                }
            }
        );

        return { number, slider, table };
    }

    /** Helper function to find the next valid structure/atom */
    private _findNextValidIndex(
        index: number,
        count: number,
        fromIndexer: (index: number) => Indexes | undefined
    ): number | undefined {
        let attempts = 0;
        while (attempts <= count) {
            index = (index + 1) % count;
            const result = fromIndexer(index);
            if (result !== undefined) {
                return index;
            }
            attempts++;
        }
        return undefined;
    }

    /** Helper function to create html number input element */
    private _createNumberInput(
        selector: string,
        max: string,
        handler: (value: number) => void
    ): HTMLInputElement {
        const number = this._root.querySelector(selector) as HTMLInputElement;
        number.max = max;
        number.onchange = () => handler(parseInt(number.value, 10) - 1);
        // Don't collapse the info table when clicking on the input field
        number.onclick = (event) => event.stopPropagation();
        return number;
    }

    /** Get the currently selected structure/atom/environment */
    private _indexes(): Indexes {
        const structure = this._structure.slider.value();
        let indexes;
        if (this._atom !== undefined) {
            const atom = this._atom.slider.value();
            indexes = this._indexer.fromStructureAtom(this._target, structure, atom);
        } else {
            assert(this._target !== 'atom');
            indexes = this._indexer.fromStructureAtom(this._target, structure);
        }
        assert(indexes !== undefined);
        return indexes;
    }
}
