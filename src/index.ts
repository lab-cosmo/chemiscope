import {Sketchmap} from "./sketchmap";
import {Viewer} from "./viewer";

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
}

function checkInput(o: any): void {
    if (!('name' in o && typeof o['name'] === 'string')) {
        throw Error("missing 'name' in Sketchviz.input");
    }

    if (!('properties' in o && typeof o['properties'] === 'object')) {
        throw Error("missing 'properties' in Sketchviz.input");
    }
    for (const key in o['properties']) {
        if (typeof o['properties'][key] !== 'object' ||
            !(typeof o['properties'][key][0] === 'string' ||
              typeof o['properties'][key][0] === 'number')) {
                  throw Error("wrong type for 'properties' in Sketchviz.input");
              }
    }

    if (!('structures' in o && typeof o['structures'] === 'object' && typeof o['structures'][0] === 'string')) {
        throw Error("missing 'structures' in Sketchviz.input");
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
}

function setup(input: SketchvizInput) {
    checkInput(input);
    const map = new Sketchmap(input.mapId, input.name, input.properties);
    const viewer = new Viewer(input.viewerId, input.j2sPath, input.structures);
    map.onSelectedUpdate((i) => viewer.showStructure(i));
}

export {
    setup,
    Sketchmap,
    Viewer,
};
