import {ScatterPlot} from "./plot";
import {StructureViewer} from "./viewer";
import {EnvironmentSlider} from './slider';
import {PropertiesTable} from './properties';

import {Dataset, checkDataset} from './dataset';
import {EnvironmentIndexer} from './indexer';

require('./static/sketchviz.css');

interface VizualizerInput extends Dataset {
    plotId: string;
    viewerId: string;
    j2sPath: string;
    sliderId: string;
    propertiesId: string;
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

    if (!('sliderId' in o && typeof o['sliderId'] === 'string')) {
        throw Error("missing 'sliderId' in Sketchviz.input");
    }

    if (!('propertiesId' in o && typeof o['propertiesId'] === 'string')) {
        throw Error("missing 'propertiesId' in Sketchviz.input");
    }
}


class Vizualizer {
    public plot: ScatterPlot;
    public viewer: StructureViewer;
    public table: PropertiesTable;
    public slider: EnvironmentSlider;

    private _ids: {
        plot: string;
        viewer: string;
        slider: string;
        table: string;
    }
    /// translate between structure/atom <=> environment
    private _indexer: EnvironmentIndexer;

    constructor(input: VizualizerInput) {
        checkInput(input);

        this._ids = {
            plot: input.plotId,
            viewer: input.viewerId,
            slider: input.sliderId,
            table: input.propertiesId,
        }

        const mode = (input.environments === undefined) ? 'structure' : 'atom';
        this._indexer = new EnvironmentIndexer(mode, input.structures, input.environments);

        this.plot = new ScatterPlot(input.plotId, input.meta.name, mode, input.properties);
        this.viewer = new StructureViewer(input.viewerId, input.j2sPath, this._indexer, input.structures, input.environments);

        this.table = new PropertiesTable(input.propertiesId, input.properties, this._indexer);
        this.slider = new EnvironmentSlider(input.sliderId, this._indexer);
        this.slider.onChange((indexes) => {
            this.plot.select(this._indexer.environment(indexes));
        });

        this.plot.onSelectedUpdate((environment: number) => {
            const indexes = this._indexer.indexes(environment);
            this.slider.changed(indexes);
            this.viewer.show(indexes);
            this.table.show(indexes);
        });
    }

    public changeDataset(dataset: Dataset) {
        checkDataset(dataset);

        const mode = (dataset.environments === undefined) ? 'structure' : 'atom';
        this._indexer = new EnvironmentIndexer(mode, dataset.structures, dataset.environments);

        this.plot.changeDataset(dataset.meta.name, mode, dataset.properties);
        this.viewer.changeDataset(this._indexer, dataset.structures, dataset.environments);

        this.table = new PropertiesTable(this._ids.table, dataset.properties, this._indexer);
        this.slider = new EnvironmentSlider(this._ids.slider, this._indexer);
        this.slider.onChange((indexes) => {
            this.plot.select(this._indexer.environment(indexes));
        });

        this.plot.onSelectedUpdate((environment: number) => {
            const indexes = this._indexer.indexes(environment);
            this.slider.changed(indexes);
            this.viewer.show(indexes);
            this.table.show(indexes);
        });
    }
}

export {
    StructureViewer,
    ScatterPlot,
    EnvironmentSlider,
    EnvironmentIndexer,
    PropertiesTable,
    Vizualizer,
};
