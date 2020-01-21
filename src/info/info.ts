/**
 * @packageDocumentation
 * @module info
 */

import assert from 'assert';

import {Indexes, EnvironmentIndexer} from '../utils';
import {Property} from '../dataset';

import {Slider} from './slider';
import {Table} from './table';

function filter<T extends Object>(obj: T, predicate: (o: any) => boolean): T {
    const result: any = {};
    for (const key in obj) {
        if (obj.hasOwnProperty(key) && predicate(obj[key])) {
            result[key] = obj[key];
        }
    }

    return result;
};

/** @hidden
 * Generate a random id with 20 characters
 */
function generateId() {
    const arr = new Uint8Array(20);
    window.crypto.getRandomValues(arr);
    return Array.from(arr, (dec) => ('0' + dec.toString(16)).substr(-2)).join('');
}

/**
 * Information associated with the current structure or atom
 */
interface Info {
    /** The button hiding/showing the table, displaying the index of the structure / atom */
    label: HTMLButtonElement;
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
    private _root: HTMLElement;
    private _atom?: Info;
    private _structure: Info;
    private _keepOrientiation: HTMLInputElement;
    private _indexer: EnvironmentIndexer;

    /** Callback used when the user changes one of the sliders value */
    public onchange: (indexes: Indexes, keepOrientation: boolean) => void;

    /**
     * Create a new [[EnvironmentInfo]] inside the DOM element with given `id`
     * @param id         HTML id of the DOM element where the sliders and
     *                   tables should live
     * @param properties properties to be displayed
     * @param indexer    [[EnvironmentIndexer]] used to translate indexes from
     *                   environments index to structure/atom indexes
     */
    constructor(id: string, properties: {[name: string]: Property}, indexer: EnvironmentIndexer) {
        const root = document.getElementById(id);
        if (root === null) {
            throw Error(`could not find HTML element #${id}`);
        }
        this._root = root;
        this._indexer = indexer;
        this.onchange = () => {};

        const structureId = 'skv-' + generateId();
        const atomId = 'skv-' + generateId();
        const settingsID = 'skv-' + generateId();

        let atomButton = '<div></div>';
        if (this._indexer.target === 'atom') {
            atomButton = `
            <button type="button" class="btn btn-sm skv-info-atom-btn"
                data-toggle="collapse"
                data-target="#${atomId}"
                aria-expanded="false"
                aria-controls="${atomId}">
                    atom 1
            </button>
            `;
        }

        this._root.innerHTML = `
        <div class="skv-info">
            <button type="button" class="btn btn-sm skv-info-structure-btn"
                data-toggle="collapse"
                data-target="#${structureId}"
                aria-expanded="false"
                aria-controls="${structureId}">
                    structure 1
            </button>
            ${atomButton}
            <button data-toggle="collapse"
                    data-target="#${settingsID}"
                    class="btn btn-light btn-sm skv-trajectory-settings-btn">
                <div class="skv-hamburger"><div></div><div></div><div></div></div>
            </button>
        </div>
        <div class="collapse skv-trajectory-settings" id="${settingsID}">
            <div class="input-group input-group-sm">
                <div class="input-group-prepend">
                    <label class="input-group-text" for=skv-playback-delay title="playback delay in tenths of seconds" style="cursor: help;">delay</label>
                </div>
                <input id=skv-playback-delay class="form-control" type="number" min=1 value=7>
            </div>

            <div class="custom-control custom-switch skv-keep-orientation">
                <input type="checkbox" class="custom-control-input" id=skv-is-trajectory>
                <label class="custom-control-label" for=skv-is-trajectory title="keep the molecule orientation" style="cursor: help;">trajectory</label>
            </div>
        </div>`;

        this._keepOrientiation = this._root.querySelector('#skv-is-trajectory') as HTMLInputElement;

        this._structure = this._createStructure(structureId, filter(properties, (p) => p.target === 'structure'));

        if (this._indexer.target === 'atom') {
            this._atom = this._createAtom(atomId, filter(properties, (p) => p.target === 'atom'));
        }
    }

    /** Create the structure slider and table */
    private _createStructure(id: string, properties: {[name: string]: Property}): Info {
        const delay = this._root.querySelector('#skv-playback-delay') as HTMLInputElement;

        const slider = new Slider(this._root, 'structure', delay);;
        const n_structures = this._indexer.structuresCount();
        slider.reset(n_structures - 1);

        const table = new Table(this._root, 'structure', id, properties);

        slider.onchange = () => {
            if (this._atom !== undefined) {
                const n_atoms = this._indexer.atomsCount(this._structure.slider.value());
                this._atom.label.innerText = `atom 1`;
                this._atom.slider.reset(n_atoms - 1);
            }

            const indexes = this._indexes();
            this._structure.table.show(indexes);
            this._structure.label.innerText = `structure ${indexes.structure + 1}`;

            if (this._atom !== undefined) {
                this._atom.table.show(indexes);
            }

            this.onchange(indexes, this._keepOrientiation.checked);
        }

        const label = this._root.getElementsByClassName('skv-info-structure-btn')[0] as HTMLButtonElement;

        return { label, slider, table };
    }

    /** Create the atom slider and table */
    private _createAtom(id: string, properties: {[name: string]: Property}) {
        const delay = this._root.querySelector('#skv-playback-delay') as HTMLInputElement;
        const slider = new Slider(this._root, 'atom', delay);
        const n_atoms = this._indexer.atomsCount(this._structure.slider.value());
        slider.reset(n_atoms - 1);
        slider.onchange = () => {
            const indexes = this._indexes();
            this._atom!.table.show(indexes);
            this._atom!.label.innerText = `atom ${indexes.atom! + 1}`;
            this.onchange(indexes, this._keepOrientiation.checked);
        }

        const table = new Table(this._root, 'atom', id, properties);

        const label = this._root.getElementsByClassName('skv-info-atom-btn')[0] as HTMLButtonElement;
        return { label, slider, table };
    }

    /** Show properties for the given `indexes`, and update the sliders values */
    public show(indexes: Indexes) {
        this._structure.label.innerText = `structure ${indexes.structure + 1}`;
        this._structure.slider.update(indexes.structure);
        this._structure.table.show(indexes);

        if (this._atom !== undefined) {
            const n_atoms = this._indexer.atomsCount(indexes.structure);
            this._atom.label.innerText = `atom 1`;
            this._atom.slider.reset(n_atoms - 1);
            this._structure.table.show(indexes);
        }

        if (indexes.atom !== undefined) {
            if (this._atom === undefined) {
                throw Error("Invalid state: got an atomic number to update, but I am displaying only structures")
            } else {
                this._atom.label.innerText = `atom ${indexes.atom + 1}`;
                this._atom.slider.update(indexes.atom);
                this._atom.table.show(indexes);
            }
        }
    }

    /** Get the currently selected structure/atom/environment */
    private _indexes(): Indexes {
        const structure = this._structure.slider.value();
        if (this._atom !== undefined) {
            const atom = this._atom.slider.value();
            return this._indexer.from_structure_atom(structure, atom);
        } else {
            assert(this._indexer.target == 'structure');
            return this._indexer.from_structure_atom(structure);
        }
    }
}
