import assert from 'assert';

import {Indexes, EnvironmentIndexer} from './indexer';
import {Slider} from './slider';
import {Table} from './table';
import {Property} from './dataset';

function filter<T extends Object>(obj: T, predicate: (o: any) => boolean): T {
    const result: any = {};
    for (const key in obj) {
        if (obj.hasOwnProperty(key) && predicate(obj[key])) {
            result[key] = obj[key];
        }
    }

    return result;
};

function generateId() {
    const arr = new Uint8Array(20);
    window.crypto.getRandomValues(arr);
    return Array.from(arr, (dec) => ('0' + dec.toString(16)).substr(-2)).join('');
}

interface Info {
    label: HTMLButtonElement;
    slider: Slider;
    table: Table;
}

/// Display information about structure or environment, using both a slider for
/// selection and an table displaying all properties (hidden by default).
export class EnvironmentInfo {
    private _root: HTMLElement;
    private _atom?: Info;
    private _structure: Info;
    private _keepOrientiation: HTMLInputElement;
    private _indexer: EnvironmentIndexer;
    public onchange: (indexes: Indexes, keepOrientation: boolean) => void;

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

        this._root.innerHTML = `<div class="skv-info">
            <button type="button" class="btn btn-sm skv-info-structure-btn"
                data-toggle="collapse"
                data-target="#${structureId}"
                aria-expanded="false"
                aria-controls="${structureId}">
                    structure 1
            </button>
            ${atomButton}
            <div class="input-group input-group-sm">
                <div class="input-group-prepend">
                    <label class="input-group-text" for=skv-playback-delay title="playback delay in tenths of seconds" style="cursor: help;">delay</label>
                </div>
                <input id=skv-playback-delay class="form-control" type="number" min=1 value=7>
            </div>

            <div class="custom-control custom-switch">
                <input type="checkbox" class="custom-control-input" id=skv-is-trajectory>
                <label class="custom-control-label" for=skv-is-trajectory title="keep the molecule orientation" style="cursor: help;">trajectory</label>
            </div>
        </div>`;

        this._keepOrientiation = this._root.querySelector('#skv-is-trajectory') as HTMLInputElement;

        this._structure = this._createStructure(structureId, filter(properties, (p) => p.target === 'structure'), indexer);

        if (this._indexer.target === 'atom') {
            this._atom = this._createAtom(atomId, filter(properties, (p) => p.target === 'atom'), indexer);
        }
    }

    private _createStructure(id: string, properties: {[name: string]: Property}, indexer: EnvironmentIndexer): Info {
        const delay = this._root.querySelector('#skv-playback-delay') as HTMLInputElement;

        const slider = new Slider(this._root, 'structure', delay);;
        const n_structures = this._indexer.structuresCount();
        slider.reset(n_structures - 1);

        const table = new Table(this._root, 'structure', id, properties, indexer);

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

    private _createAtom(id: string, properties: {[name: string]: Property}, indexer: EnvironmentIndexer) {
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

        const table = new Table(this._root, 'atom', id, properties, indexer);

        const label = this._root.getElementsByClassName('skv-info-atom-btn')[0] as HTMLButtonElement;
        return { label, slider, table };
    }

    /// The environment index changed outside, update the sliders
    public select({structure, atom}: Indexes) {
        this._structure.label.innerText = `structure ${structure + 1}`;
        this._structure.slider.update(structure);
        this._structure.table.show({structure, atom});

        if (this._atom !== undefined) {
            const n_atoms = this._indexer.atomsCount(structure);
            this._atom.label.innerText = `atom 1`;
            this._atom.slider.reset(n_atoms - 1);
            this._structure.table.show({structure, atom});
        }

        if (atom !== undefined) {
            if (this._atom === undefined) {
                throw Error("Invalid state: got an atomic number to update, but I am displaying only structures")
            } else {
                this._atom.label.innerText = `atom ${atom + 1}`;
                this._atom.slider.update(atom);
                this._atom.table.show({structure, atom});
            }
        }
    }

    private _indexes(): Indexes {
        const structure = this._structure.slider.value();
        if (this._atom !== undefined) {
            const atom = this._atom.slider.value();
            return {structure, atom};
        } else {
            assert(this._indexer.target == 'structure');
            return {structure};
        }
    }
}
