import assert from 'assert';
import JSmolViewer from 'materials-cloud-viewer';

import {structure2JSmol} from './jsmol';
import {Structure, Environment} from './dataset';
import {Indexes, EnvironmentIndexer} from './indexer';

export class StructureViewer {
    private _viewer: JSmolViewer;
    private _structures: string[];
    private _environments?: Environment[];
    private _indexer: EnvironmentIndexer;
    /// index of the currently displayed structure/atom
    private _current: Indexes;

    constructor(id: string, j2sPath: string, indexer: EnvironmentIndexer, structures: Structure[], environments?: Environment[]) {
        this._viewer = new JSmolViewer(id, j2sPath);
        this._structures = structures.map(structure2JSmol);
        this._environments = environments;
        this._indexer = indexer;
        this._current = {structure: -1, atom: -1};
        this.show({structure: 0, atom: 0});
    }

    public changeDataset(indexer: EnvironmentIndexer, structures: Structure[], environments?: Environment[]) {
        this._structures = structures.map(structure2JSmol);
        this._environments = environments;
        this._indexer = indexer;
        this._current = {structure: -1, atom: -1};
        this.show({structure: 0, atom: 0});
    }

    public show(indexes: Indexes) {
        if (this._current.structure !== indexes.structure) {
            assert(indexes.structure < this._structures.length);
            this._viewer.load(`inline '${this._structures[indexes.structure]}'`);
            // Force atom to be updated
            this._current.atom = -1;
        }

        if (this._indexer.mode === 'atom' && this._current.atom != indexes.atom) {
            this._showEnvironment(this._indexer.environment(indexes));
        }

        this._current = indexes;
    }

    public computeSettingsPlacement(callback: (rect: DOMRect) => {top: number, left: number}) {
        this._viewer.computeSettingsPlacement(callback)
    }

    private _showEnvironment(environment: number) {
        const env = this._environments![environment];
        // Default style
        this._viewer.script("select all;")
        this._viewer.script("wireframe 0.15; spacefill off; dots off;")
        this._viewer.script("color atoms translucent cpk;")

        // Style for the atoms in the environment
        const selection = env.neighbors.map((i) => `@${i + 1}`).join(' or ');
        this._viewer.script(`select ${selection};`)
        this._viewer.script("wireframe 0.15; spacefill 23%; dots off;")
        this._viewer.script("color atoms cpk;")

        // Style for the central atom
        this._viewer.script(`select @${env.center + 1};`)
        this._viewer.script("wireframe off; spacefill 23%; dots 0.6;")
        this._viewer.script("color atoms green;")

        // Reset selection
        this._viewer.script("select all;")
    }
}
