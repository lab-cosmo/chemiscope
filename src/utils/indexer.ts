/**
 * @packageDocumentation
 * @module utils
 */

import assert from 'assert';
import {Structure, Environment, Target} from '../dataset';

/**
 * Indexes of the
 */
export interface Indexes {
    /** The global environment index */
    environment: number;
    /** Index of the structure which the environment corresponds to */
    structure: number;
    /**
     * Index of the atom in the structure which corresponds to the environment.
     *
     * If we are considering full-structures environments only, this is
     * `undefined`.
     */
    atom?: number;
}

/// Map between environment index and structure/atom pairs
export class EnvironmentIndexer {
    public target: Target;

    private _structures: Structure[];
    private _environments?: Environment[];
    /// Number of atoms before the structure at the same index
    private _atomsBefore: number[];

    constructor(target: Target, structures: Structure[], environments?: Environment[]) {
        this.target = target;
        this._structures = structures;
        this._environments = environments;

        let current = 0;
        this._atomsBefore = this._structures.map((structure) => {
            let result = current;
            current += structure.names.length
            return result;
        })
    }

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

    public environmentsCount(): number {
        if (this._environments !== undefined) {
            return this._environments.length;
        } else {
            return this._structures.length;
        }
    }

    public structuresCount(): number {
        return this._structures.length;
    }

    public atomsCount(structure: number): number {
        return this._structures[structure].names.length;
    }
}
