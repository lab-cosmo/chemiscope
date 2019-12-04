import assert from 'assert';
import {Structure, Environment, Mode} from './dataset';

export interface Indexes {
    structure: number;
    atom?: number;
}

/// Map between environment index and structure/atom pairs
export class EnvironmentIndexer {
    public mode: Mode;

    private _structures: Structure[];
    private _environments?: Environment[];
    /// Number of atoms before the structure at the same index
    private _atomsBefore: number[];

    constructor(mode: Mode, structures: Structure[], environments?: Environment[]) {
        this.mode = mode;
        this._structures = structures;
        this._environments = environments;

        let current = 0;
        this._atomsBefore = this._structures.map((structure) => {
            let result = current;
            current += structure.names.length
            return result;
        })
    }

    public indexes(environment: number): Indexes {
        if (this.mode == 'structure') {
            return {
                structure: environment,
            }
        } else {
            assert(this.mode === 'atom');
            assert(this._environments !== undefined);

            const env = this._environments![environment];
            return {
                structure: env.structure,
                atom: env.center,
            }
        }
    }

    public environment(indexes: Indexes): number {
        if (this.mode == 'structure') {
            assert(indexes.atom === undefined || indexes.atom === 0);
            return indexes.structure
        } else {
            assert(this.mode === 'atom');
            assert(this._environments !== undefined);

            // assume that environments are ordered structure by structure, then
            // by atom in the structure.$
            const {structure, atom} = indexes;
            const environment = this._atomsBefore[structure] + atom!;
            assert(this._environments![environment].center === atom);
            assert(this._environments![environment].structure === structure);
            return environment;
        }
    }

    public environmentsCount(): number | undefined {
        if (this._environments !== undefined) {
            return this._environments.length;
        }
    }

    public structuresCount(): number {
        return this._structures.length;
    }

    public atomsCount(structure: number): number {
        return this._structures[structure].names.length;
    }
}
