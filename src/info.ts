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
    slider: Slider;
    table: Table;
}

/// Display information about structure or environment, using both a slider for
/// selection and an table displaying all properties (hidden by default).
export class EnvironmentInfo {
    private _root: HTMLElement;
    private _atom?: Info;
    private _structure!: Info;
    private _indexer: EnvironmentIndexer;
    public onchange: (indexes: Indexes) => void;

    constructor(id: string, properties: {[name: string]: Property}, indexer: EnvironmentIndexer) {
        const root = document.getElementById(id);
        if (root === null) {
            throw Error(`could not find HTML element #${id}`);
        }
        this._root = root;
        this._root.innerHTML = '';
        this._indexer = indexer;
        this.onchange = () => {};

        this._createStructure(filter(properties, (p) => p.target === 'structure'), indexer);

        if (this._indexer.target === 'atom') {
            this._createAtom(filter(properties, (p) => p.target === 'atom'), indexer);
        }
    }

    private _createStructure(properties: {[name: string]: Property}, indexer: EnvironmentIndexer) {
        const id = 'skv-' + generateId();
        const slider = new Slider('structure', this._root, id);
        const n_structures = this._indexer.structuresCount();
        slider.reset(n_structures - 1);

        const table = new Table(this._root, 'structure', id, properties, indexer);

        slider.onchange = () => {
            if (this._atom !== undefined) {
                const n_atoms = this._indexer.atomsCount(this._structure.slider.value());
                this._atom.slider.reset(n_atoms - 1);
            }

            const indexes = this._indexes();
            this._structure.table.show(indexes);
            if (this._atom !== undefined) {
                this._atom.table.show(indexes);
            }

            this.onchange(indexes);
        }

        this._structure = { slider, table };
    }

    private _createAtom(properties: {[name: string]: Property}, indexer: EnvironmentIndexer) {
        const id = 'skv-' + generateId();
        const slider = new Slider('atom', this._root, id);
        const n_atoms = this._indexer.atomsCount(this._structure.slider.value());
        slider.reset(n_atoms - 1);
        slider.onchange = () => {
            const indexes = this._indexes();
            table.show(indexes)
            this.onchange(indexes);
        }

        const table = new Table(this._root, 'atom', id, properties, indexer);
        this._atom = { slider, table };
    }

    /// The environment index changed outside, update the slider
    public select({structure, atom}: Indexes) {
        this._structure.slider.update(structure);
        this._structure.table.show({structure, atom});

        if (this._atom !== undefined) {
            const n_atoms = this._indexer.atomsCount(structure);
            this._atom.slider.reset(n_atoms - 1);
            this._structure.table.show({structure, atom});
        }

        if (atom !== undefined) {
            if (this._atom === undefined) {
                throw Error("Invalid state: got an atomic number to update, but I am displaying only structures")
            } else {
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
