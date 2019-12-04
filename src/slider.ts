import assert from 'assert';
import {Indexes, EnvironmentIndexer} from './indexer';

function createSlider(name: string): [HTMLElement, HTMLInputElement, HTMLElement] {
    const template = document.createElement('template');
    template.innerHTML = `<div class="input-group input-group-sm">
        <div class="input-group-prepend">
            <label class="input-group-text skv-slider-label">${name}</label>
            <span class="input-group-text"><div class="skv-play-button"></div></span>
        </div>
        <input class="form-control custom-range" type='range' min=0 value=0 step=1></input>
    </div>`;
    const group = template.content.firstChild! as HTMLElement;
    const slider = group.querySelector('input')! as HTMLInputElement;
    const play = group.querySelector('.skv-play-button')! as HTMLElement;
    return [group, slider, play];
}

export class EnvironmentSlider {
    private _root: HTMLElement;
    private _atom?: HTMLInputElement;
    private _structure!: HTMLInputElement;
    private _onchange: (indexes: Indexes) => void;
    private _indexer: EnvironmentIndexer;

    constructor(id: string, indexer: EnvironmentIndexer) {
        const root = document.getElementById(id);
        if (root === null) {
            throw Error(`could not find HTML element #${id}`);
        }
        this._root = root;
        this._root.innerHTML = '';
        this._onchange = () => {};
        this._indexer = indexer;

        this._createStructureSlider();
        if (this._indexer.mode === 'atom') {
            this._createAtomSlider();
        }
    }

    private _createStructureSlider() {
        const [group, slider, play]  = createSlider('structure');
        this._root.append(group);

        this._structure = slider;
        const n_structure = this._indexer.structuresCount();
        this._structure.max = (n_structure - 1).toString();
        this._structure.onchange = () => {
            if (this._atom !== undefined) {
                this._atom.value = "0";
                const n_atoms = this._indexer.atomsCount(parseInt(this._structure.value));
                this._atom.max = (n_atoms - 1).toString();
            }
            this._onchange(this._indexes());
        }

        play.onclick = () => {
            play.classList.toggle('skv-playing');
            this._playStructure();
        }
    }

    private _createAtomSlider() {
        const [group, slider, play]  = createSlider('atom');
        this._root.append(group);

        this._atom = slider;
        const n_atoms = this._indexer.atomsCount(0);
        this._atom.max = (n_atoms - 1).toString();
        this._atom.onchange = () => {
            this._onchange(this._indexes());
        }

        play.onclick = () => {
            play.classList.toggle('skv-playing');
            this._playAtom();
        }
    }

    /// Call the given callback when the user change the slider value
    public onChange(callback: (indexes: Indexes) => void) {
        this._onchange = callback;
    }

    /// The environment index changed outside, update the slider
    public changed({structure, atom}: Indexes) {
        this._structure.value = structure.toString();
        if (atom !== undefined) {
            if (this._atom === undefined) {
                throw Error("Invalid state: got an atomic number to update, but I am displaying only structures")
            } else {
                this._atom.value = atom.toString();
            }
        }
    }

    private _playStructure() {
        const play = this._structure.parentElement!.querySelector('.skv-play-button')!;
        setTimeout(() => {
            if (play.classList.contains('skv-playing')) {
                // change value
                let structure = parseInt(this._structure.value)
                structure = (structure + 1) % this._indexer.structuresCount();
                this._structure.value = structure.toString();
                if (this._atom !== undefined) {
                    this._atom.value = "0";
                    const n_atoms = this._indexer.atomsCount(structure);
                    this._atom.max = (n_atoms - 1).toString();
                }
                this._onchange(this._indexes());

                // contibue playing until 'skv-playing' is no longer there
                this._playStructure();
            }
        }, 750)
    }

    private _playAtom() {
        const play = this._atom!.parentElement!.querySelector('.skv-play-button')!;
        setTimeout(() => {
            if (play.classList.contains('skv-playing')) {
                // change value
                const structure = parseInt(this._structure.value);
                let atom = parseInt(this._atom!.value)
                atom = (atom + 1) % this._indexer.atomsCount(structure);
                this._atom!.value = atom.toString();

                this._onchange(this._indexes());

                // contibue playing until 'skv-playing' is no longer there
                this._playAtom();
            }
        }, 750)
    }

    private _indexes(): Indexes {
        const structure = parseInt(this._structure.value);
        if (this._indexer.mode == 'atom') {
            const atom = parseInt(this._atom!.value);
            return {structure, atom};
        } else {
            assert(this._indexer.mode == 'structure');
            return {structure};
        }
    }
}
