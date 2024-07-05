/**
 * @packageDocumentation
 * @module info
 */

import assert from 'assert';

import { Parameter, Property } from '../dataset';
import { EnvironmentIndexer, Indexes } from '../indexer';
import { binarySearch, getElement, sendWarning } from '../utils';

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
    public onchange: (indexes: Indexes) => void;
    /** delay in ms between "frames" during playback over structures/atoms */
    public playbackDelay: number;

    private _shadow: ShadowRoot;
    private _root: HTMLElement;
    private _atom?: Info;
    private _structure?: Info;
    private _indexer: EnvironmentIndexer;
    private _properties: { [name: string]: Property };
    private _parameters: { [name: string]: Parameter } | undefined;

    /**
     * Create a new {@link EnvironmentInfo} inside the DOM element with given `id`
     * @param element    HTML element or string 'id' of the element where
     *                   the sliders and tables should live
     * @param properties properties to be displayed
     * @param indexer    {@link EnvironmentIndexer} used to translate indexes from
     *                   environments index to structure/atom indexes
     * @param viewer     {@link ViewersGrid} from which we get the playback delay
     */
    constructor(
        element: string | HTMLElement,
        properties: { [name: string]: Property },
        indexer: EnvironmentIndexer,
        parameters?: { [name: string]: Parameter }
    ) {
        const containerElement = getElement(element);

        const hostElement = document.createElement('div');
        containerElement.appendChild(hostElement);

        this._shadow = hostElement.attachShadow({ mode: 'open' });
        this._shadow.adoptedStyleSheets = [
            styles.bootstrap,
            styles.chemiscope,
            plotlyStyles.globalStyleSheet,
        ];

        this._root = document.createElement('div');
        this._shadow.appendChild(this._root);

        this._indexer = indexer;
        this.onchange = () => {};
        // 700 ms between steps is a good default
        this.playbackDelay = 700;

        this._properties = properties;
        this._parameters = parameters;
        this.togglePerAtom();
    }

    /** Proceed with the panels related to the mode */
    public togglePerAtom(): void {
        // Construct atom related html if mode is atom
        const atomButton =
            this._indexer.mode === 'atom'
                ? `<div class='btn btn-sm chsp-info-atom-btn'
                data-bs-toggle='collapse'
                data-bs-target='#atom'
                aria-expanded='false'
                aria-controls='atom'>
                <div class="chsp-info-btns-svg">${INFO_SVG}</div>
                atom <input class='chsp-info-number' type=number value=1 min=1></input>
            </div>`
                : '<div></div>';

        // Construct main HTML structure
        this._root.innerHTML = `
        <div class='accordion chsp-info-tables' id='info-tables'></div>
        <div class='chsp-info-btns'>
            <div class='btn btn-sm chsp-info-structure-btn'
                data-bs-toggle='collapse'
                data-bs-target='#structure'
                aria-expanded='false'
                aria-controls='structure'>
                <div class="chsp-info-btns-svg">${INFO_SVG}</div>
                structure <input class='chsp-info-number' type=number value=1 min=1></input>
            </div>
            ${atomButton}
        </div>`;

        // Filter properties for structure
        const structureProperties = filter(this._properties, (p) => p.target === 'structure');
        this._structure = this._createStructure(structureProperties, this._parameters);

        // Initialize central atom if in atom mode
        if (this._indexer.mode === 'atom') {
            const atomProperties = filter(this._properties, (p) => p.target === 'atom');
            this._atom = this._createAtom(atomProperties, this._parameters);
        }

        // Initialize the collapse components from their 'data-bs-*' attributes.
        Collapse.initialize(this._root);
    }

    /** Show properties for the given `indexes`, and update the sliders values */
    public show(indexes: Indexes): void {
        assert(this._structure !== undefined);
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

    /** Create the structure slider and table */
    private _createStructure(
        properties: { [name: string]: Property },
        parameters?: { [name: string]: Parameter }
    ): Info {
        const slider = new Slider(this._root, 'structure');
        slider.reset(this._indexer.activeStructures());

        const tableRoot = this._root.children[0] as HTMLElement;
        assert(tableRoot.tagName.toLowerCase() === 'div');
        const table = new Table(tableRoot, 'structure', 'structure', properties, parameters);
        slider.startPlayback = (advance) => {
            setTimeout(() => {
                if (advance()) {
                    // find the next structure/atom to display
                    const current = this._indexes();
                    const structuresCount = this._indexer.structuresCount();

                    let iterations = 0;
                    let structure = current.structure;
                    let indexes = undefined;
                    while (indexes === undefined) {
                        structure = (structure + 1) % structuresCount;
                        // try to find the first atom in the structure with
                        // associated data
                        for (let atom = 0; atom < this._indexer.atomsCount(structure); atom++) {
                            indexes = this._indexer.from_structure_atom(structure, atom);
                            if (indexes !== undefined) {
                                break;
                            }
                        }
                        // prevent infinite loop
                        if (iterations === structuresCount) {
                            return;
                        }
                        iterations += 1;
                    }

                    assert(indexes !== undefined);
                    this.show(indexes);
                    this.onchange(indexes);
                    // continue playing until the advance callback returns false
                    slider.startPlayback(advance);
                }
            }, this.playbackDelay);
        };

        slider.onchange = () => {
            assert(this._structure !== undefined);
            const structure = this._structure.slider.value();

            if (this._atom !== undefined) {
                const activeAtoms = this._indexer.activeAtoms(structure);
                if (activeAtoms.length === 0) {
                    sendWarning(
                        `can not change to structure ${
                            structure + 1
                        } which does not contain any active atom`
                    );
                    return;
                }
                this._atom.number.value = `${activeAtoms[0] + 1}`;
                this._atom.number.max = `${activeAtoms.length}`;
                this._atom.slider.reset(activeAtoms);
            }

            const indexes = this._indexes();
            assert(indexes !== undefined);

            this._structure.table.show(indexes);
            this._structure.number.value = `${indexes.structure + 1}`;

            if (this._atom !== undefined) {
                this._atom.table.show(indexes);
            }

            this.onchange(indexes);
        };

        const number = this._root.querySelector(
            '.chsp-info-structure-btn .chsp-info-number'
        ) as HTMLInputElement;
        number.max = this._indexer.structuresCount().toString();

        number.onchange = () => {
            assert(this._structure !== undefined);
            const value = parseInt(number.value, 10) - 1;
            if (isNaN(value) || value < 0 || value >= parseInt(number.max, 10)) {
                // reset to the current slider value if we got an invalid value
                number.value = `${this._structure.slider.value() + 1}`;
            } else if (binarySearch(this._indexer.activeStructures(), value) === -1) {
                // also reset if we got a value which is not an active structure
                number.value = `${this._structure.slider.value() + 1}`;
            } else {
                this._structure.slider.update(value);
                if (this._atom !== undefined) {
                    const activeAtoms = this._indexer.activeAtoms(value);
                    this._atom.number.value = `${activeAtoms[0] + 1}`;
                    this._atom.number.max = `${activeAtoms.length}`;
                    this._atom.slider.reset(activeAtoms);
                }

                const indexes = this._indexes();
                this._structure.table.show(indexes);
                if (this._atom !== undefined) {
                    this._atom.table.show(indexes);
                }

                this.onchange(indexes);
            }
        };

        // Don't collapse the info table when clicking on the input field
        number.onclick = (event) => event.stopPropagation();

        return { number, slider, table };
    }

    /** Create the atom slider and table */
    private _createAtom(
        properties: { [name: string]: Property },
        parameters?: { [name: string]: Parameter }
    ) {
        assert(this._structure !== undefined);
        const slider = new Slider(this._root, 'atom');
        slider.reset(this._indexer.activeAtoms(this._structure.slider.value()));

        slider.startPlayback = (advance) => {
            setTimeout(() => {
                if (advance()) {
                    const current = this._indexes();
                    const structure = current.structure;

                    const atomsCount = this._indexer.atomsCount(structure);
                    let atom = current.atom;
                    assert(atom !== undefined);

                    let indexes = undefined;
                    let iterations = 0;
                    while (indexes === undefined) {
                        atom = (atom + 1) % atomsCount;
                        indexes = this._indexer.from_structure_atom(structure, atom);

                        // prevent infinite loop if the current structure has no
                        // environments
                        if (iterations === atomsCount) {
                            return;
                        }
                        iterations += 1;
                    }

                    this.show(indexes);
                    this.onchange(indexes);
                    // continue playing until the advance callback returns false
                    slider.startPlayback(advance);
                }
            }, this.playbackDelay);
        };

        slider.onchange = () => {
            assert(this._atom !== undefined);
            assert(this._structure !== undefined);
            const indexes = this._indexes();

            if (indexes === undefined) {
                const structure = this._structure.slider.value();
                const atom = this._atom.slider.value();
                sendWarning(
                    `environment for atom ${atom} in structure ${structure} is not part of this dataset`
                );
                return;
            }

            assert(indexes.atom !== undefined);
            this._atom.table.show(indexes);
            this._atom.number.value = `${indexes.atom + 1}`;
            this.onchange(indexes);
        };

        const tableRoot = this._root.children[0] as HTMLElement;
        assert(tableRoot.tagName.toLowerCase() === 'div');
        const table = new Table(tableRoot, 'atom', 'atom', properties, parameters);

        const number = this._root.querySelector(
            '.chsp-info-atom-btn .chsp-info-number'
        ) as HTMLInputElement;
        number.max = this._indexer.atomsCount(this._structure.slider.value()).toString();

        number.onchange = () => {
            assert(this._atom !== undefined);
            assert(this._structure !== undefined);
            const activeAtoms = this._indexer.activeAtoms(this._structure.slider.value());

            const value = parseInt(number.value, 10) - 1;
            if (isNaN(value) || value < 0 || value >= parseInt(number.max, 10)) {
                // reset to the current slider value if we got an invalid value
                number.value = `${this._atom.slider.value() + 1}`;
            } else if (binarySearch(activeAtoms, value) === -1) {
                // also reset if we got a value which is not an active structure
                number.value = `${this._atom.slider.value() + 1}`;
            } else {
                this._atom.slider.update(value);
                const indexes = this._indexes();
                this._atom.table.show(indexes);
                this.onchange(indexes);
            }
        };

        // Don't collapse the info table when clicking on the input field
        number.onclick = (event) => event.stopPropagation();

        return { number, slider, table };
    }

    /** Get the currently selected structure/atom/environment */
    private _indexes(): Indexes {
        assert(this._structure !== undefined);
        const structure = this._structure.slider.value();
        let indexes;
        if (this._atom !== undefined) {
            const atom = this._atom.slider.value();
            indexes = this._indexer.from_structure_atom(structure, atom);
        } else {
            assert(this._indexer.mode === 'structure');
            indexes = this._indexer.from_structure_atom(structure);
        }
        assert(indexes !== undefined);
        return indexes;
    }
}
