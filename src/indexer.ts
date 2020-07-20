/**
 * @packageDocumentation
 * @module utils
 */

import assert from 'assert';
import {Environment, Structure, UserStructure} from './dataset';

/**
 * If a dataset contains both atomic and structure properties, we can only
 * display one kind at the time, indicated by the [[DisplayMode]].
 */
export type DisplayMode = 'structure' | 'atom';

/** Indexes related to a single environment */
export interface Indexes {
    /** The global environment index. */
    environment: number;
    /**
     * Index of the structure which the [[Indexes.environment|environment]]
     * corresponds to.
     */
    structure: number;
    /**
     * Index of the atom in the structure which corresponds to the environment.
     *
     * If we are considering full-structures environments only, this is
     * `undefined`.
     */
    atom?: number;
}

/**
 * [[EnvironmentIndexer]] links environment index and structure/atom indexes
 *
 * Environments can be either full structures or centered on a specific atom.
 * This class makes the link between two representations: a single, global,
 * environment index, used by the map; and the structure/atom pair, used by
 * the structure viewer and the general information panel.
 */
export class EnvironmentIndexer {
    /**
     * Current [[DisplayMode]]. This is useful for datasets that
     * contain both atom-level and structure-level properties.
     */
    public mode: DisplayMode;

    private _structures: Structure[] | UserStructure[];
    private _environments?: Environment[];
    /// Number of atoms before the structure at the same index
    private _atomsBefore: number[];

    /**
     * Create a new [[EnvironmentIndexer]] for the given set of structures and
     * environments.
     *
     * @param mode         should we display atomic or structure properties
     * @param structures   structures used in the current dataset
     * @param environments environments used in the current dataset
     */
    constructor(mode: DisplayMode, structures: Structure[] | UserStructure[], environments?: Environment[]) {
        this.mode = mode;
        this._structures = structures;
        this._environments = environments;

        if (this.mode === 'atom') {
            assert(this._environments !== undefined);
        }

        let count = 0;
        this._atomsBefore = [];
        for (const structure of this._structures) {
            this._atomsBefore.push(count);
            count += structure.size;
        }
    }

    /**
     * Get a full set of indexes from the global environment index
     * @param  environment global index of an environment
     * @return             full [[Indexes]], containing the corresponding
     *                     structure / atom indexes
     */
    public from_environment(environment: number): Indexes {
        if (this.mode === 'structure') {
            return {
                environment: environment,
                structure: environment,
            };
        } else {
            assert(this.mode === 'atom');
            assert(this._environments !== undefined);

            const env = this._environments[environment];
            return {
                atom: env.center,
                environment: environment,
                structure: env.structure,
            };
        }
    }

    /**
     * Get a full set of indexes from the structure/atom indexes
     * @param  structure index of the structure in the full structure list
     * @param  atom      index of the atom in the structure
     * @return           full [[Indexes]], containing the global environment
     *                   index
     */
    public from_structure_atom(structure: number, atom?: number): Indexes {
        if (this.mode === 'structure') {
            assert(atom === undefined || atom === 0);
            return {
                environment: structure,
                structure: structure,
            };
        } else {
            assert(this.mode === 'atom');
            assert(this._environments !== undefined);
            assert(atom !== undefined);

            // assume that environments are ordered structure by structure, then
            // by atom in the structure.
            const environment = this._atomsBefore[structure] + atom;
            assert(this._environments[environment].center === atom);
            assert(this._environments[environment].structure === structure);
            return {
                atom: atom,
                environment: environment,
                structure: structure,
            };
        }
    }

    /** Get the total number of environments currently being displayed */
    public environmentsCount(): number {
        if (this.mode === 'atom') {
            assert(this._environments !== undefined);
            return this._environments.length;
        } else {
            return this._structures.length;
        }
    }

    /** Get the total number of structures we know about */
    public structuresCount(): number {
        return this._structures.length;
    }

    /** Get the total number of atom in the `structure` with given index */
    public atomsCount(structure: number): number {
        return this._structures[structure].size;
    }
}
