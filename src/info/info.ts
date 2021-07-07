/**
 * @packageDocumentation
 * @module info
 */

import assert from 'assert';

import { Property } from '../dataset';
import { EnvironmentIndexer, Indexes } from '../indexer';
import { generateGUID, getElement } from '../utils';

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
    /**
     * Callback fired when the use click the play button of the structure slider.
     * The `advance` callback indicate whether to continue playback or not.
     */
    public startStructurePlayback: (advance: () => boolean) => void;
    /**
     * Callback fired when the use click the play button of the atom slider.
     * The `advance` callback indicate whether to continue playback or not.
     */
    public startAtomPlayback: (advance: () => boolean) => void;

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
        this.startStructurePlayback = () => {};
        this.startAtomPlayback = () => {};

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
        this._structure.number.value = `${indexes.structure + 1}`;
        this._structure.slider.update(indexes.structure);
        this._structure.table.show(indexes);

        if (this._atom !== undefined) {
            const n_atoms = this._indexer.atomsCount(indexes.structure);
            this._atom.number.value = '1';
            this._atom.slider.reset(n_atoms - 1);
            this._structure.table.show(indexes);
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
        const n_structures = this._indexer.structuresCount();
        slider.reset(n_structures - 1);

        const tableRoot = this._root.children[0] as HTMLElement;
        assert(tableRoot.tagName.toLowerCase() === 'div');
        const table = new Table(tableRoot, 'structure', id, properties);

        slider.startPlayback = (advance) => this.startStructurePlayback(advance);
        slider.onchange = () => {
            if (this._atom !== undefined) {
                const n_atoms = this._indexer.atomsCount(this._structure.slider.value());
                this._atom.number.value = '1';
                this._atom.number.max = n_atoms.toString();
                this._atom.slider.reset(n_atoms - 1);
            }

            const indexes = this._indexes();
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
        // Don't collapse the info table when clicking on the input field
        number.onclick = (event) => event.stopPropagation();
        number.onchange = () => {
            const value = parseInt(number.value, 10) - 1;
            if (isNaN(value) || value < 0 || value >= parseInt(number.max, 10)) {
                // reset to the current slider value if we got an invalid value
                number.value = `${this._structure.slider.value() + 1}`;
            } else {
                this._structure.slider.update(value);
                if (this._atom !== undefined) {
                    const n_atoms = this._indexer.atomsCount(value);
                    this._atom.slider.reset(n_atoms - 1);
                    this._atom.number.value = '1';
                    this._atom.number.max = n_atoms.toString();
                }

                const indexes = this._indexes();
                this._structure.table.show(indexes);
                if (this._atom !== undefined) {
                    this._atom.table.show(indexes);
                }

                this.onchange(indexes);
            }
        };

        return { number, slider, table };
    }

    /** Create the atom slider and table */
    private _createAtom(id: string, properties: { [name: string]: Property }) {
        const slider = new Slider(this._root, 'atom');
        const n_atoms = this._indexer.atomsCount(this._structure.slider.value());
        slider.reset(n_atoms - 1);
        slider.startPlayback = (advance) => this.startAtomPlayback(advance);
        slider.onchange = () => {
            assert(this._atom !== undefined);
            const indexes = this._indexes();
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
        number.max = n_atoms.toString();
        // Don't collapse the info table when clicking on the input field
        number.onclick = (event) => event.stopPropagation();
        number.onchange = () => {
            assert(this._atom !== undefined);
            const value = parseInt(number.value, 10) - 1;
            if (isNaN(value) || value < 0 || value >= parseInt(number.max, 10)) {
                // reset to the current slider value if we got an invalid value
                number.value = `${this._atom.slider.value() + 1}`;
            } else {
                this._atom.slider.update(value);
                const indexes = this._indexes();
                this._atom.table.show(indexes);

                this.onchange(indexes);
            }
        };

        return { number, slider, table };
    }

    /** Get the currently selected structure/atom/environment */
    private _indexes(): Indexes {
        const structure = this._structure.slider.value();
        if (this._atom !== undefined) {
            const atom = this._atom.slider.value();
            return this._indexer.from_structure_atom(structure, atom);
        } else {
            assert(this._indexer.mode === 'structure');
            return this._indexer.from_structure_atom(structure);
        }
    }
}
