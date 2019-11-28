import {Sketchmap} from "./sketchmap";
import {Viewer} from "./viewer";
import {EnvironementSlider} from './slider';

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
}

function setup(input: SketchvizInput) {
    checkInput(input);

    const map = new Sketchmap(input.mapId, input.name, input.properties);
    const viewer = new Viewer(input.viewerId, input.j2sPath, input.structures);
    const slider = new EnvironementSlider(input.sliderID, input.structures.length - 1);
    map.onSelectedUpdate((i) => {
        slider.changed(i);
        viewer.showStructure(i);
    });
    slider.onChange((i) => map.select(i));
}

export {
    setup,
    Sketchmap,
    Viewer,
};
