/**
 * @packageDocumentation
 * @module main
 */

import {PropertiesMap} from "./map";
import {StructureViewer} from "./structure";
import {EnvironmentInfo} from './info';

import {Dataset, checkDataset} from './dataset';
import {EnvironmentIndexer} from './indexer';

require('./static/sketchviz.css');

interface VizualizerInput extends Dataset {
    plotId: string;
    viewerId: string;
    j2sPath: string;
    infoId: string;
}

function checkInput(o: any) {
    checkDataset(o);

    if (!('plotId' in o && typeof o['plotId'] === 'string')) {
        throw Error("missing 'plotId' in Sketchviz.input");
    }

    if (!('viewerId' in o && typeof o['viewerId'] === 'string')) {
        throw Error("missing 'viewerId' in Sketchviz.input");
    }

    if (!('j2sPath' in o && typeof o['j2sPath'] === 'string')) {
        throw Error("missing 'j2sPath' in Sketchviz.input");
    }

    if (!('infoId' in o && typeof o['infoId'] === 'string')) {
        throw Error("missing 'infoId' in Sketchviz.input");
    }
}


class Vizualizer {
    public plot: PropertiesMap;
    public viewer: StructureViewer;
    public info: EnvironmentInfo;

    private _ids: {
        plot: string;
        viewer: string;
        info: string;
    }
    /// translate between structure/atom <=> environment
    private _indexer: EnvironmentIndexer;

    constructor(input: VizualizerInput) {
        checkInput(input);

        this._ids = {
            plot: input.plotId,
            viewer: input.viewerId,
            info: input.infoId,
        }

        const mode = (input.environments === undefined) ? 'structure' : 'atom';
        this._indexer = new EnvironmentIndexer(mode, input.structures, input.environments);

        this.plot = new PropertiesMap(input.plotId, input.meta.name, mode, input.properties);
        this.viewer = new StructureViewer(input.viewerId, input.j2sPath, this._indexer, input.structures, input.environments);

        if (mode === 'atom') {
            this.viewer.onselect = (indexes) => {
                this.plot.select(this._indexer.environment(indexes));
                this.info.select(indexes);
            };
        }

        this.info = new EnvironmentInfo(input.infoId, input.properties, this._indexer);
        this.info.onchange = (indexes, keepOrientation) => {
            this.plot.select(this._indexer.environment(indexes));
            this.viewer.show(indexes, keepOrientation);
        };

        this.plot.onselect = (environment: number) => {
            const indexes = this._indexer.indexes(environment);
            this.info.select(indexes);
            this.viewer.show(indexes);
        };
    }

    public changeDataset(dataset: Dataset) {
        checkDataset(dataset);

        const mode = (dataset.environments === undefined) ? 'structure' : 'atom';
        this._indexer = new EnvironmentIndexer(mode, dataset.structures, dataset.environments);

        this.plot.changeDataset(dataset.meta.name, mode, dataset.properties);
        this.viewer.changeDataset(this._indexer, dataset.structures, dataset.environments);
        if (mode === 'atom') {
            this.viewer.onselect = (indexes) => {
                this.plot.select(this._indexer.environment(indexes));
                this.info.select(indexes);
            };
        } else {
            this.viewer.onselect = () => {};
        }

        this.info = new EnvironmentInfo(this._ids.info, dataset.properties, this._indexer);
        this.info.onchange = (indexes) => {
            this.plot.select(this._indexer.environment(indexes));
            this.viewer.show(indexes);
        };

        this.plot.onselect = (environment: number) => {
            const indexes = this._indexer.indexes(environment);
            this.info.select(indexes);
            this.viewer.show(indexes);
        };
    }
}

export {
    StructureViewer,
    PropertiesMap,
    EnvironmentInfo,
    EnvironmentIndexer,
    Vizualizer,
};
