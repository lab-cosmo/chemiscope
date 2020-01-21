/**
 * @packageDocumentation
 * @module utils
 */

import assert from 'assert';
import {Structure, Environment, Target} from '../dataset';

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
     * Current [[Target]] being displayed. This is useful for datasets that
     * contain both atom-level and structure-level properties.
     */
    public target: Target;

    private _structures: Structure[];
    private _environments?: Environment[];
    /// Number of atoms before the structure at the same index
    private _atomsBefore: number[];

    /**
     * Create a new [[EnvironmentIndexer]] for the given set of structures and
     * environments.
     *
     * @param target       are we dealing with `'atom'` or `'structure'`
     *                     properties?
     * @param structures   structures used in the current dataset
     * @param environments environments used in the current dataset
     */
    constructor(target: Target, structures: Structure[], environments?: Environment[]) {
        this.target = target;
        this._structures = structures;
        this._environments = environments;

        if (this.target == "atom") {
            assert(this._environments !== undefined);
        }

        let current = 0;
        this._atomsBefore = this._structures.map((structure) => {
            let result = current;
            current += structure.names.length
            return result;
        })
    }

    /**
     * Get a full set of indexes from the global environment index
     * @param  environment global index of an environment
     * @return             full [[Indexes]], containing the corresponding
     *                     structure / atom indexes
     */
    public from_environment(environment: number): Indexes {
        if (this.target == 'structure') {
            return {
                environment: environment,
                structure: environment,
            }
        } else {
            assert(this.target === 'atom');
            assert(this._environments !== undefined);

            const env = this._environments![environment];
            return {
                environment: environment,
                structure: env.structure,
                atom: env.center,
            }
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
        if (this.target == 'structure') {
            assert(atom === undefined || atom === 0);
            return {
                environment: structure,
                structure: structure,
            }
        } else {
            assert(this.target === 'atom');
            assert(this._environments !== undefined);

            // assume that environments are ordered structure by structure, then
            // by atom in the structure.
            const environment = this._atomsBefore[structure] + atom!;
            assert(this._environments![environment].center === atom);
            assert(this._environments![environment].structure === structure);
            return {
                environment: environment,
                structure: structure,
                atom: atom,
            }
        }
    }

    /** Get the total number of environments we know about */
    public environmentsCount(): number {
        if (this._environments !== undefined) {
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
        return this._structures[structure].names.length;
    }
}
