/**
 * @packageDocumentation
 * @module structure
 */

import assert from 'assert';

import {Environment, isPositiveInteger, Structure, UserStructure} from '../dataset';
import {EnvironmentIndexer, Indexes} from '../utils';

import {structure2JSmol} from './jsmol';
import {JSmolWidget} from './widget';

function groupByStructure(n_structures: number, environments?: Environment[]): Environment[][] | undefined {
    if (environments === undefined) {
        return undefined;
    }

    const result: Environment[][] = [];
    for (let i = 0; i < n_structures; i++) {
        result.push([]);
    }

    for (const env of environments) {
        result[env.structure].push(env);
    }

    return result;
}

/**
 * Check that the given object is a structure. Return a string describing the
 * issue with `s` if any, or the empty string if `s` looks like a valid
 * structure.
 */
function checkStructure(s: any): string {
    if (!('size' in s && typeof s.size === 'number' && isPositiveInteger(s.size))) {
        return 'missing "size" in structure';
    }

    for (const key of ['names', 'x', 'y', 'z']) {
        if (!(key in s && s[key].length !== undefined)) {
            return `missing "${name}" in structure`;
        }
        if (s[key].length !== s.size) {
            return `wrong size for "${name}" in structure, expected ${s.size}, got ${s[name].length}`;
        }
    }

    if ('cell' in s) {
        if (s.cell.length !== 9) {
            return '"cell" must be an array of size 9 in structure';
        }
    }

    return '';
}

/**
 * The [[StructureViewer]] class displays a molecule or a crystal in 3D using
 * [JSmol](http://wiki.jmol.org/index.php/JSmol) for rendering.
 */
export class StructureViewer {
    /** Callback used when the user select an environment */
    public onselect: (indexes: Indexes) => void;

    /**
     * Callback used when a new structure should be loaded
     *
     * By default, this assumes that the loaded dataset contains [[Structure]]
     * directly, and returns the data from there. If the loaded dataset contains
     * [[UserStructure]] instead, this callback should be set to transform from
     * data in [[UserStructure.data]] to a [[Structure]].
     *
     * The callback gets two parameter: the structure index (0-based); and the
     * full [[UserStructure]].
     */
    public loadStructure: (index: number, structure: any) => Structure;

    private _widget: JSmolWidget;
    /// Playback delay setting
    private _delay: HTMLInputElement;
    /// List of structures in the dataset
    private _structures: Structure[] | UserStructure[];
    /// Cached string representation of structures
    private _cachedStructures: string[];
    /// Optional list of environments for each structure
    private _environments?: Environment[][];
    private _indexer: EnvironmentIndexer;
    /// index of the currently displayed structure/atom
    private _current: {structure: number; atom?: number};

    /**
     * Create a new [[StructureViewer]] inside the HTML element with the given
     * `id`
     *
     * @param id           HTML id of the DOM element where the viewer should live
     * @param j2sPath      path to the `j2s` files uses by JSmol
     * @param indexer      [[EnvironmentIndexer]] used to translate indexes from
     *                     environments index to structure/atom indexes
     * @param structures   list of structure to display
     * @param environments list of atom-centered environments in the structures,
     *                     used to highlight the selected environment
     */
    constructor(
        id: string,
        j2sPath: string,
        indexer: EnvironmentIndexer,
        structures: Structure[] | UserStructure[],
        environments?: Environment[],
    ) {
        this._widget = new JSmolWidget(id, j2sPath);
        this._delay = document.getElementById(`${this._widget.guid}-playback-delay`) as HTMLInputElement;
        this._structures = structures;
        this._cachedStructures = new Array(structures.length);
        this._environments = groupByStructure(this._structures.length, environments);
        this._indexer = indexer;
        this._current = {structure: -1, atom: -1};

        this.loadStructure = (_, s) => {
            // check that the data does conform to the Structure interface
            if (checkStructure(s) !== '') {
                throw Error(
                    'got custom data for this structure, but no custom loadStructure callback\n' +
                    `the object was ${JSON.stringify(s)}`,
                );
            } else {
                return s;
            }
        };

        this.onselect = () => {};

        this._widget.onselect = (atom: number) => {
            if (this._indexer.mode !== 'atom') {
                return;
            }

            this._widget.highlight(atom);
            this._current.atom = atom;
            // if the viewer is showing a bigger supercell than [1, 1, 1], the
            // atom index can be outside of [0, natoms), so make sure it is
            // inside this range.
            const atom_id = atom % this._widget.natoms()!;
            const indexes = this._indexer.from_structure_atom(this._current.structure, atom_id);
            this.onselect(indexes);
        };
    }

    /**
     * Change the displayed dataset to a new one, without re-creating the
     * viewer itself.
     *
     * @param  indexer      new indexer making the environment index to
     *                      structure/atom pair translation
     * @param  structures   new list of structures to display
     * @param  environments new list of atom centered environments
     */
    public changeDataset(
        indexer: EnvironmentIndexer,
        structures: Structure[] | UserStructure[],
        environments?: Environment[],
    ) {
        this._structures = structures;
        this._cachedStructures = new Array(structures.length);
        this._environments = groupByStructure(this._structures.length, environments);
        this._indexer = indexer;
        this._current = {structure: -1, atom: -1};
    }

    /**
     * Show a new structure, as identified by `indexes`. This will switch to
     * the structure at index `indexes.structure`, and if environments where
     * passed to the constructor and the current display mode is `"atom"`,
     * highlight the atom-centered environment corresponding to `indexes.atom`.
     *
     * @param  indexes         structure / atom pair to display
     */
    public show(indexes: Indexes) {
        if (this._current.structure !== indexes.structure) {
            assert(indexes.structure < this._structures.length);
            const options = {
                packed: false,
                trajectory: true,
            } as any;

            if (this._environments !== undefined) {
                options.environments = this._environments[indexes.structure];
                if (this._indexer.mode === 'atom') {
                    options.highlight = indexes.atom;
                }
            }

            this._widget.load(`inline '${this._structureForJSmol(indexes.structure)}'`, options);
        }

        if (this._indexer.mode === 'atom') {
            if (this._current.atom !== indexes.atom) {
                this._widget.highlight(indexes.atom);
            }
        } else {
            this._widget.highlight(undefined);
        }

        this._current = indexes;
    }

    /**
     * Register a `callback` to compute the placement of the settings modal.
     *
     * The callback gets the current placement of the settings as a
     * [DOMRect](https://developer.mozilla.org/en-US/docs/Web/API/DOMRect),
     * and should return top and left positions in pixels, used with `position:
     * fixed`. The callback is called once, the first time the settings are
     * opened.
     */
    public settingsPlacement(callback: (rect: DOMRect) => {top: number, left: number}) {
        this._widget.settingsPlacement(callback);
    }

    /**
     * Start playing the trajectory of structures in this dataset, until
     * `advance` returns false
     */
    public structurePlayback(advance: () => boolean) {
        setTimeout(() => {
            if (advance()) {
                const structure = (this._current.structure + 1) % this._indexer.structuresCount();
                const indexes = this._indexer.from_structure_atom(structure, 0);
                this.show(indexes);
                this.onselect(indexes);
                // continue playing until the advance callback returns false
                this.structurePlayback(advance);
            }
        }, parseFloat(this._delay.value) * 100);
    }

    /**
     * Start playing the 'trajectory' of atoms in the current structure, until
     * `advance` returns false
     */
    public atomPlayback(advance: () => boolean) {
        setTimeout(() => {
            if (advance()) {
                const structure = this._current.structure;
                const atom = (this._current.atom! + 1) % this._indexer.atomsCount(structure);
                const indexes = this._indexer.from_structure_atom(structure, atom);
                this.show(indexes);
                this.onselect(indexes);
                // continue playing until the advance callback returns false
                this.atomPlayback(advance);
            }
        }, parseFloat(this._delay.value) * 100);
    }

    private _structureForJSmol(index: number): string {
        if (this._cachedStructures[index] === undefined) {
            const s = this.loadStructure(index, this._structures[index]);
            const check = checkStructure(s);
            if (check !== '') {
                throw Error(
                    `got invalid object as structure: ${check}` + '\n' +
                    `the object was ${JSON.stringify(s)}`,
                );
            }
            this._cachedStructures[index] = structure2JSmol(s);
        }
        return this._cachedStructures[index];
    }
}
