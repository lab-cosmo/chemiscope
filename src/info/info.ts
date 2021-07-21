/**
 * @packageDocumentation
 * @module info
 */

import assert from 'assert';

import { Property } from '../dataset';
import { EnvironmentIndexer, Indexes } from '../indexer';
import { binarySearch, generateGUID, getElement, sendWarning } from '../utils';

import { Slider } from './slider';
import { Table } from './table';

import INFO_SVG from '../static/info.svg';

function filter<T extends Record<string, Property>>(
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
 * The [[EnvironmentInfo]] class displays information about structure or
 * environment, using a slider for selection and an hidden by default table
 * displaying all properties.
 */
export class EnvironmentInfo {
    /** Callback used when the user changes one of the sliders value */
    public onchange: (indexes: Indexes) => void;
    /** delay in ms between "frames" during playback over structures/atoms */
    public playbackDelay: number;

    private _root: HTMLElement;
    private _atom?: Info;
    private _structure: Info;
    private _indexer: EnvironmentIndexer;

    /**
     * Create a new [[EnvironmentInfo]] inside the DOM element with given `id`
     * @param element    HTML element or string 'id' of the element where
     *                   the sliders and tables should live
     * @param properties properties to be displayed
     * @param indexer    [[EnvironmentIndexer]] used to translate indexes from
     *                   environments index to structure/atom indexes
     * @param viewer     [[ViewersGrid]] from which we get the playback delay
     */
    constructor(
        element: string | HTMLElement,
        properties: { [name: string]: Property },
        indexer: EnvironmentIndexer
    ) {
        this._root = getElement(element);

        this._indexer = indexer;
        this.onchange = () => {};
        // 700 ms between steps is a good default
        this.playbackDelay = 700;

        const structureId = `chsp-${generateGUID()}`;
        const atomId = `chsp-${generateGUID()}`;

        let atomButton = '<div></div>';
        if (this._indexer.mode === 'atom') {
            atomButton = `
            <div class='btn btn-sm chsp-info-atom-btn'
                data-toggle='collapse'
                data-target='#${atomId}'
                aria-expanded='false'
                aria-controls='${atomId}'>
                <div class="chsp-info-btns-svg">${INFO_SVG}</div>
                atom <input class='chsp-info-number' type=number value=1 min=1></input>
            </div>
            `;
        }

        this._root.innerHTML = `
        <div class='accordion chsp-info-tables' id='info-tables'></div>
        <div class='chsp-info-btns'>
            <div class='btn btn-sm chsp-info-structure-btn'
                data-toggle='collapse'
                data-target='#${structureId}'
                aria-expanded='false'
                aria-controls='${structureId}'>
                <div class="chsp-info-btns-svg">${INFO_SVG}</div>
                structure <input class='chsp-info-number' type=number value=1 min=1></input>

            </div>
            ${atomButton}
        </div>`;

        const structureProperties = filter(properties, (p) => p.target === 'structure');
        this._structure = this._createStructure(structureId, structureProperties);

        if (this._indexer.mode === 'atom') {
            const atomProperties = filter(properties, (p) => p.target === 'atom');
            this._atom = this._createAtom(atomId, atomProperties);
        }
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
     * Remove all HTML added by this [[EnvironmentInfo]] in the current document
     */
    public remove(): void {
        this._root.innerHTML = '';
    }

    /** Create the structure slider and table */
    private _createStructure(id: string, properties: { [name: string]: Property }): Info {
        const slider = new Slider(this._root, 'structure');
        slider.reset(this._indexer.activeStructures());

        const tableRoot = this._root.children[0] as HTMLElement;
        assert(tableRoot.tagName.toLowerCase() === 'div');
        const table = new Table(tableRoot, 'structure', id, properties);

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
    private _createAtom(id: string, properties: { [name: string]: Property }) {
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
        const table = new Table(tableRoot, 'atom', id, properties);

        const number = this._root.querySelector(
            '.chsp-info-atom-btn .chsp-info-number'
        ) as HTMLInputElement;
        number.max = this._indexer.atomsCount(this._structure.slider.value()).toString();

        number.onchange = () => {
            assert(this._atom !== undefined);
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
