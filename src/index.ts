import {ScatterPlot} from "./plot";
import {StructureViewer} from "./viewer";
import {EnvironementSlider} from './slider';
import {PropertiesTable} from './properties';

require('./static/sketchviz.css');

interface SketchvizInput {
    name: string;
    properties: {
        [name: string]: string[] | number[];
    };
    structures: string[];
    mapId: string;
    viewerId: string;
    j2sPath: string;
    sliderID: string;
    propertiesID: string;
}

function checkInput(o: any): void {
    if (!('name' in o && typeof o['name'] === 'string')) {
        throw Error("missing 'name' in Sketchviz.input");
    }

    if (!('properties' in o && typeof o['properties'] === 'object')) {
        throw Error("missing 'properties' in Sketchviz.input");
    }
    for (const key in o['properties']) {
        if (o['properties'][key].length === undefined ||
            !(typeof o['properties'][key][0] === 'string' ||
              typeof o['properties'][key][0] === 'number')) {
                  throw Error("wrong type for 'properties' in Sketchviz.input");
              }
    }

    if (!('structures' in o && o['structures'].length !== undefined && typeof o['structures'][0] === 'string')) {
        throw Error("missing 'structures' in Sketchviz.input");
    }

    // Check that properties and structures have the same size
    const size = o['structures'].length;
    for (const key in o['properties']) {
        if (o['properties'][key].length !== size) {
            throw Error("properties and structures don't have the same size")
        }
    }

    if (!('mapId' in o && typeof o['mapId'] === 'string')) {
        throw Error("missing 'mapId' in Sketchviz.input");
    }

    if (!('viewerId' in o && typeof o['viewerId'] === 'string')) {
        throw Error("missing 'viewerId' in Sketchviz.input");
    }

    if (!('j2sPath' in o && typeof o['j2sPath'] === 'string')) {
        throw Error("missing 'j2sPath' in Sketchviz.input");
    }

    if (!('sliderID' in o && typeof o['sliderID'] === 'string')) {
        throw Error("missing 'sliderID' in Sketchviz.input");
    }

    if (!('propertiesID' in o && typeof o['propertiesID'] === 'string')) {
        throw Error("missing 'propertiesID' in Sketchviz.input");
    }
}


class Vizualizer {
    private _viewer: StructureViewer;
    private _plot: ScatterPlot;
    private _slider: EnvironementSlider;
    private _table: PropertiesTable;

    constructor(input: SketchvizInput) {
        checkInput(input);

        this._plot = new ScatterPlot(input.mapId, input.name, input.properties);
        this._viewer = new StructureViewer(input.viewerId, input.j2sPath, input.structures);
        this._slider = new EnvironementSlider(input.sliderID, input.structures.length - 1);
        this._table = new PropertiesTable(input.propertiesID, input.properties);
        this._plot.onSelectedUpdate((i) => {
            this._slider.changed(i);
            this._viewer.showStructure(i);
            this._table.display(i);
        });
        this._slider.onChange((i) => this._plot.select(i));
    }
}

export {
    StructureViewer,
    ScatterPlot,
    EnvironementSlider,
    PropertiesTable,
    Vizualizer,
};
