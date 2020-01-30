/**
 * This module contains the [[DefaultVizualizer]] class, which is the default
 * entry point of the code.
 *
 * The default vizualization is organized around three panels: the
 * [[PropertiesMap|map]] (a scatter plot of properties), the
 * [[StructureViewer|structure viewer]], and the general dataset
 * [[EnvironmentInfo|information]]. Each one of these is defined in a
 * separate module.
 *
 * Other organization of the vizualization are possible by using the classes
 * responsible for each sub-panel separately, instead of using the
 * [[DefaultVizualizer]]. In this case, users should make sure to finish the
 * plumbing by setting the right callbacks on each element used.
 *
 * @packageDocumentation
 * @module main
 * @preferred
 */

import {PropertiesMap} from "./map";
import {StructureViewer} from "./structure";
import {EnvironmentInfo} from './info';

import {Dataset, checkDataset} from './dataset';
import {EnvironmentIndexer, addWarningHandler} from './utils';

require('./static/chemiscope.css');

/**
 * Configuration for the [[DefaultVizualizer]]
 */
export interface Config {
    /** Id of the DOM element to use for the [[PropertiesMap|map]] */
    map: string;
    /** Id of the DOM element to use for the [[EnvironmentInfo|environment information]] */
    info: string;
    /** Id of the DOM element to use for the [[StructureViewer|structure viewer]] */
    structure: string;
    /** Path of j2s files, used by JSmol, which is used by the [[StructureViewer|structure viewer]] */
    j2sPath: string;
}

/** @hidden
 * Check if `o` contains all the expected fields to be a [[Config]].
 */
function checkConfig(o: any) {
    if (!('map' in o && typeof o['map'] === 'string')) {
        throw Error("missing 'map' in chemiscope config");
    }

    if (!('info' in o && typeof o['info'] === 'string')) {
        throw Error("missing 'info' in chemiscope config");
    }

    if (!('structure' in o && typeof o['structure'] === 'string')) {
        throw Error("missing 'structure' in chemiscope config");
    }

    if (!('j2sPath' in o && typeof o['j2sPath'] === 'string')) {
        throw Error("missing 'j2sPath' in chemiscope config");
    }
}

// return new Promise((resolve, reject) => {
//     const res = longrunning()
//     resolve('I don\'t actually like fish ' + res)
// })

class DefaultVizualizer {
    public map: PropertiesMap;
    public info: EnvironmentInfo;
    public structure: StructureViewer;

    private _ids: {
        map: string;
        info: string;
        structure: string;
    }
    private _indexer: EnvironmentIndexer;

    private constructor(config: Config, dataset: Dataset) {
        checkConfig(config);
        checkDataset(dataset);

        this._ids = config;
        const mode = (dataset.environments === undefined) ? 'structure' : 'atom';
        this._indexer = new EnvironmentIndexer(mode, dataset.structures, dataset.environments);

        this.map = new PropertiesMap(config.map, dataset.meta.name, this._indexer, dataset.properties);
        this.structure = new StructureViewer(config.structure, config.j2sPath, this._indexer, dataset.structures, dataset.environments);

        this.structure.onselect = (indexes) => {
            this.map.select(indexes);
            this.info.show(indexes);
        };

        this.info = new EnvironmentInfo(config.info, dataset.properties, this._indexer);
        this.info.onchange = (indexes) => {
            this.map.select(indexes);
            this.structure.show(indexes);
        };
        this.info.startStructurePlayback = (advance) => this.structure.structurePlayback(advance);
        this.info.startAtomPlayback = (advance) => this.structure.atomPlayback(advance);

        this.map.onselect = (indexes) => {
            this.info.show(indexes);
            this.structure.show(indexes);
        };
    }

    static load(config: Config, dataset: Dataset): Promise<DefaultVizualizer> {
        return new Promise((resolve, _) => {
            const visualizer = new DefaultVizualizer(config, dataset);
            resolve(visualizer);
        });
    }

    public changeDataset(dataset: Dataset): Promise<void> {
        return new Promise((resolve, _) => {
            checkDataset(dataset);

            const mode = (dataset.environments === undefined) ? 'structure' : 'atom';
            this._indexer = new EnvironmentIndexer(mode, dataset.structures, dataset.environments);

            this.map.changeDataset(dataset.meta.name, this._indexer, dataset.properties);
            this.structure.changeDataset(this._indexer, dataset.structures, dataset.environments);
            this.structure.onselect = (indexes) => {
                this.map.select(indexes);
                this.info.show(indexes);
            };

            this.info = new EnvironmentInfo(this._ids.info, dataset.properties, this._indexer);
            this.info.onchange = (indexes) => {
                this.map.select(indexes);
                this.structure.show(indexes);
            };
            this.info.startStructurePlayback = (advance) => this.structure.structurePlayback(advance);
            this.info.startAtomPlayback = (advance) => this.structure.atomPlayback(advance);


            this.map.onselect = (indexes) => {
                this.info.show(indexes);
                this.structure.show(indexes);
            };
            resolve();
        });
    }
}

export {
    addWarningHandler,
    StructureViewer,
    PropertiesMap,
    EnvironmentInfo,
    EnvironmentIndexer,
    DefaultVizualizer,
};
