import {ScatterPlot} from "./plot";
import {StructureViewer} from "./viewer";
import {EnvironementSlider} from './slider';
import {PropertiesTable} from './properties';

import {Dataset, checkDataset} from './dataset';

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
    public viewer: StructureViewer;
    public plot: ScatterPlot;
    public slider: EnvironementSlider;
    public table: PropertiesTable;
    private _ids: {
        plot: string;
        viewer: string;
        slider: string;
        table: string;
    }

    constructor(input: VizualizerInput) {
        checkInput(input);

        this._ids = {
            plot: input.plotId,
            viewer: input.viewerId,
            slider: input.sliderId,
            table: input.propertiesId,
        }

        this.plot = new ScatterPlot(input.plotId, input.meta.name, input.properties);
        this.viewer = new StructureViewer(input.viewerId, input.j2sPath, input.structures);
        this.slider = new EnvironementSlider(input.sliderId, input.structures.length - 1);
        this.table = new PropertiesTable(input.propertiesId, input.properties);
        this.plot.onSelectedUpdate((i) => {
            this.slider.changed(i);
            this.viewer.showStructure(i);
            this.table.display(i);
        });
        this.slider.onChange((i) => this.plot.select(i));
    }

    public changeDataset(dataset: Dataset) {
        checkDataset(dataset);

        this.plot.changeDataset(dataset.meta.name, dataset.properties);
        this.viewer.changeDataset(dataset.structures);
        this.slider.changeDataset(dataset.structures.length - 1);
        this.table = new PropertiesTable(this._ids.table, dataset.properties);
    }
}

export {
    StructureViewer,
    ScatterPlot,
    EnvironementSlider,
    PropertiesTable,
    Vizualizer,
};
