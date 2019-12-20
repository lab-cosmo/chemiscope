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
    /// List of environments for each structure
    private _environments?: Environment[][];
    private _indexer: EnvironmentIndexer;
    /// index of the currently displayed structure/atom
    private _current: Indexes;

    constructor(id: string, j2sPath: string, indexer: EnvironmentIndexer, structures: Structure[], environments?: Environment[]) {
        this._viewer = new JSmolViewer(id, j2sPath);
        this._structures = structures.map(structure2JSmol);
        this._environments = groupByStructure(this._structures.length, environments);
        this._indexer = indexer;
        this._current = {structure: -1, atom: -1};
        this.show({structure: 0, atom: 0});
    }

    public changeDataset(indexer: EnvironmentIndexer, structures: Structure[], environments?: Environment[]) {
        this._structures = structures.map(structure2JSmol);
        this._environments = groupByStructure(this._structures.length, environments);
        this._indexer = indexer;
        this._current = {structure: -1, atom: -1};
        this.show({structure: 0, atom: 0});
    }

    public show(indexes: Indexes) {
        if (this._current.structure !== indexes.structure) {
            assert(indexes.structure < this._structures.length);
            const options = {
                packed: false,
            } as any;

            if (this._environments !== undefined) {
                options.environments = this._environments[indexes.structure];
            }

            this._viewer.load(`inline '${this._structures[indexes.structure]}'`, options);
            // Force atom to be updated
            this._current.atom = -1;
        }

        if (this._indexer.mode === 'atom' && this._current.atom != indexes.atom) {
            this._viewer.highlight(indexes.atom)
        }

        this._current = indexes;
    }

    public settingsPlacement(callback: (rect: DOMRect) => {top: number, left: number}) {
        this._viewer.settingsPlacement(callback)
    }

    public onSelected(callback: (indexes: Indexes) => void) {
        this._viewer.onSelected((atom: number) => {
            callback({structure: this._current.structure, atom});
        })
    }
}
