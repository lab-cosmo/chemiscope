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
import {EnvironmentIndexer} from './utils';

require('./static/sketchviz.css');

/**
 * Configuration for the [[DefaultVizualizer]]
 */
export interface Config {
    /** Id of the DOM element to use for the [[PropertiesMap|map]] */
    map: string;
    /** Id of the DOM element to use for the [[StructureViewer|structure viewer]] */
    viewer: string;
    /** Path of j2s files, used by JSmol, which is used by the [[StructureViewer|structure viewer]] */
    j2sPath: string;
    /** Id of the DOM element to use for the [[EnvironmentInfo|environment information]] */
    info: string;
}

/** @hidden
 * Check if `o` contains all the expected fields to be a [[Config]].
 */
function checkConfig(o: any) {
    if (!('map' in o && typeof o['map'] === 'string')) {
        throw Error("missing 'map' in Sketchviz.input");
    }

    if (!('viewer' in o && typeof o['viewer'] === 'string')) {
        throw Error("missing 'viewer' in Sketchviz.input");
    }

    if (!('j2sPath' in o && typeof o['j2sPath'] === 'string')) {
        throw Error("missing 'j2sPath' in Sketchviz.input");
    }

    if (!('info' in o && typeof o['info'] === 'string')) {
        throw Error("missing 'info' in Sketchviz.input");
    }
}

class DefaultVizualizer {
    public viewer: StructureViewer;
    public info: EnvironmentInfo;
    public map: PropertiesMap;

    private _ids: {
        viewer: string;
        info: string;
        map: string;
    }
    private _indexer: EnvironmentIndexer;

    constructor(config: Config, dataset: Dataset) {
        checkConfig(config);
        checkDataset(dataset);

        this._ids = config;
        const mode = (dataset.environments === undefined) ? 'structure' : 'atom';
        this._indexer = new EnvironmentIndexer(mode, dataset.structures, dataset.environments);

        this.map = new PropertiesMap(config.map, dataset.meta.name, this._indexer, dataset.properties);
        this.viewer = new StructureViewer(config.viewer, config.j2sPath, this._indexer, dataset.structures, dataset.environments);

        if (mode === 'atom') {
            this.viewer.onselect = (indexes) => {
                this.map.select(indexes);
                this.info.show(indexes);
            };
        }

        this.info = new EnvironmentInfo(config.info, dataset.properties, this._indexer);
        this.info.onchange = (indexes, keepOrientation) => {
            this.map.select(indexes);
            this.viewer.show(indexes, keepOrientation);
        };

        this.map.onselect = (indexes) => {
            this.info.show(indexes);
            this.viewer.show(indexes);
        };
    }

    public changeDataset(dataset: Dataset) {
        checkDataset(dataset);

        const mode = (dataset.environments === undefined) ? 'structure' : 'atom';
        this._indexer = new EnvironmentIndexer(mode, dataset.structures, dataset.environments);

        this.map.changeDataset(dataset.meta.name, this._indexer, dataset.properties);
        this.viewer.changeDataset(this._indexer, dataset.structures, dataset.environments);
        if (mode === 'atom') {
            this.viewer.onselect = (indexes) => {
                this.map.select(indexes);
                this.info.show(indexes);
            };
        } else {
            this.viewer.onselect = () => {};
        }

        this.info = new EnvironmentInfo(this._ids.info, dataset.properties, this._indexer);
        this.info.onchange = (indexes) => {
            this.map.select(indexes);
            this.viewer.show(indexes);
        };

        this.map.onselect = (indexes) => {
            this.info.show(indexes);
            this.viewer.show(indexes);
        };
    }
}

export {
    StructureViewer,
    PropertiesMap,
    EnvironmentInfo,
    EnvironmentIndexer,
    DefaultVizualizer,
};
