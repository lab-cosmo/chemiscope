import assert from 'assert';
import JSmolViewer from 'materials-cloud-viewer';

import {structure2JSmol} from './jsmol';
import {Structure, Environment} from './dataset';
import {Indexes, EnvironmentIndexer} from './indexer';

function groupByStructure(n_structures: number, environments?: Environment[]): Environment[][] | undefined {
    if (environments === undefined) {
        return undefined;
    }

    const result: Environment[][] = [];
    for (let i=0; i<n_structures; i++) {
        result.push([]);
    }

    for (const env of environments) {
        result[env.structure].push(env);
    }

    return result;
}

export class StructureViewer {
    private _viewer: JSmolViewer;
    /// List of structures in the dataset
    private _structures: string[];
    /// Optional list of environments for each structure
    private _environments?: Environment[][];
    private _indexer: EnvironmentIndexer;
    /// index of the currently displayed structure/atom
    private _current: Indexes;

    public onselect: (indexes: Indexes) => void;

    constructor(id: string, j2sPath: string, indexer: EnvironmentIndexer, structures: Structure[], environments?: Environment[]) {
        this._viewer = new JSmolViewer(id, j2sPath);
        this._structures = structures.map(structure2JSmol);
        this._environments = groupByStructure(this._structures.length, environments);
        this._indexer = indexer;
        this._current = {structure: -1, atom: -1};
        this.show({structure: 0, atom: 0});

        this.onselect = () => {};

        this._viewer.onselect = (atom: number) => {
            if (this._indexer.target == 'atom') {
                this._viewer.highlight(atom);
            }
            this.onselect({structure: this._current.structure, atom});
        };
    }

    public changeDataset(indexer: EnvironmentIndexer, structures: Structure[], environments?: Environment[]) {
        this._structures = structures.map(structure2JSmol);
        this._environments = groupByStructure(this._structures.length, environments);
        this._indexer = indexer;
        this._current = {structure: -1, atom: -1};
        this.show({structure: 0, atom: 0});
    }

    public show(indexes: Indexes, keepOrientation = false) {
        if (this._current.structure !== indexes.structure) {
            assert(indexes.structure < this._structures.length);
            const options = {
                packed: false,
                keepOrientation: keepOrientation,
            } as any;

            if (this._environments !== undefined) {
                options.environments = this._environments[indexes.structure];
                if (this._indexer.target === 'atom') {
                    options.highlight = indexes.atom;
                }
            }

            this._viewer.load(`inline '${this._structures[indexes.structure]}'`, options);
        }

        if (this._indexer.target === 'atom') {
            if  (this._current.atom != indexes.atom) {
                this._viewer.highlight(indexes.atom)
            }
        } else {
            this._viewer.highlight(undefined);
        }

        this._current = indexes;
    }

    public settingsPlacement(callback: (rect: DOMRect) => {top: number, left: number}) {
        this._viewer.settingsPlacement(callback)
    }
}
