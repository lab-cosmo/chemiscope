import assert from 'assert';

import {Indexes, EnvironmentIndexer} from './indexer';
import {Slider} from './slider';

export class EnvironmentInfo {
    private _root: HTMLElement;
    private _atom?: Slider;
    private _structure: Slider;
    private _indexer: EnvironmentIndexer;
    public onchange: (indexes: Indexes) => void;

    constructor(id: string, indexer: EnvironmentIndexer) {
        const root = document.getElementById(id);
        if (root === null) {
            throw Error(`could not find HTML element #${id}`);
        }
        this._root = root;
        this._root.innerHTML = '';
        this._indexer = indexer;
        this.onchange = () => {};

        this._structure = new Slider('structure');
        this._root.append(this._structure.root);

        const n_structures = this._indexer.structuresCount();
        this._structure.reset(n_structures - 1);
        this._structure.onchange = () => {
            if (this._atom !== undefined) {
                const n_atoms = this._indexer.atomsCount(this._structure.value());
                this._atom.reset(n_atoms - 1);
            }
            this.onchange(this._indexes());
        }

        if (this._indexer.mode === 'atom') {
            this._atom = new Slider('atom');
            this._root.append(this._atom.root);

            const n_atoms = this._indexer.atomsCount(this._structure.value());
            this._atom.reset(n_atoms - 1);
            this._atom.onchange = () => this.onchange(this._indexes());
        }
    }

    /// The environment index changed outside, update the slider
    public changed({structure, atom}: Indexes) {
        this._structure.update(structure);

        if (this._atom !== undefined) {
            const n_atoms = this._indexer.atomsCount(structure);
            this._atom.reset(n_atoms - 1);
        }

        if (atom !== undefined) {
            if (this._atom === undefined) {
                throw Error("Invalid state: got an atomic number to update, but I am displaying only structures")
            } else {
                this._atom.update(atom);
            }
        }
    }

    private _indexes(): Indexes {
        const structure = this._structure.value();
        if (this._atom !== undefined) {
            const atom = this._atom.value();
            return {structure, atom};
        } else {
            assert(this._indexer.mode == 'structure');
            return {structure};
        }
    }
}
