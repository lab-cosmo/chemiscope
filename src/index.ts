/**
 * This module contains the [[DefaultVisualizer]] class, which is the default
 * entry point of the code.
 *
 * The default visualization is organized around four panels: the
 * [[MetadataPanel|metadata]] panel; the [[PropertiesMap|map]] (a scatter plot of
 * properties), the [[ViewersGrid|structure viewer]], and the general
 * dataset [[EnvironmentInfo|information]]. Each one of these is defined in a
 * separate module.
 *
 * Other organization of the visualization are possible by using the classes
 * responsible for each sub-panel separately, instead of using the
 * [[DefaultVisualizer]]. In this case, developers should make sure to finish
 * the plumbing by setting the right callbacks on each element used.
 *
 * @packageDocumentation
 * @module main
 * @preferred
 */

import {EnvironmentInfo} from './info';
import {MapPresets, PropertiesMap} from './map';
import {MetadataPanel} from './metadata';
import {StructurePresets, ViewersGrid} from './structure';

import {Dataset, JsObject, Structure, validateDataset} from './dataset';
import {addWarningHandler, EnvironmentIndexer, getNextColor} from './utils';

// tslint:disable-next-line: no-var-requires
require('./static/chemiscope.css');

/**
 * Configuration for the [[DefaultVisualizer]]
 */
export interface Config {
    /** Id of the DOM element to use for the [[MetadataPanel|metadata display]] */
    meta: string;
    map: {
        /** Id of the DOM element to use for the [[PropertiesMap|map]] */
        id: string;
        /** Presets for map settings */
        presets: MapPresets;
    };
    /** Id of the DOM element to use for the [[EnvironmentInfo|environment information]] */
    info: string;
    /** Id of the DOM element to use for the [[ViewersGrid|structure viewer]] */
    structure: {
        id: string;
        presets: StructurePresets;
    };

    /** Path of j2s files, used by JSmol, which is used by the [[StructureViewer|structure viewer]] */
    j2sPath: string;
    /** Custom structure loading callback, used to set [[ViewersGrid.loadStructure]] */
    loadStructure?: (index: number, structure: any) => Structure;
}

/** @hidden
 * Check if `o` contains all the expected fields to be a [[Config]].
 */
function validateConfig(o: JsObject) {
    if (!('meta' in o && typeof o.meta === 'string')) {
        throw Error('missing "meta" key in chemiscope config');
    }

    if (!('map' in o && 'id' in (o.map as JsObject) && typeof (o.map as JsObject).id === 'string')) {
        throw Error('missing "map.id" key in chemiscope config');
    }

    if (!('info' in o && typeof o.info === 'string')) {
        throw Error('missing "info" key in chemiscope config');
    }

    if (!('structure' in o && 'id' in (o.structure as JsObject) && typeof (o.structure  as JsObject).id === 'string')) {
        throw Error('missing "structure.id" key in chemiscope config');
    }

    if (!('j2sPath' in o && typeof o.j2sPath === 'string')) {
        throw Error('missing "j2sPath" key in chemiscope config');
    }

    // from underscore.js / https://stackoverflow.com/a/6000016
    const isFunction = (obj: any) => {
        return !!(obj && obj.constructor && obj.call && obj.apply);
    };

    if ('loadStructure' in o && o.loadStructure !== undefined && !isFunction(o.loadStructure)) {
        throw Error('"loadStructure" should be a function in chemiscope config');
    }
}

/**
 * The default visualization state of chemiscope: three panels (map, structure,
 * info) updating one another when the user interact with any of them.
 */
class DefaultVisualizer {
    /**
     * Load a dataset and create a visualizer.
     *
     * This function returns a `Promise<DefaultVisualizer>` to prevent blocking
     * the browser while everything is loading.
     *
     * @param  config  configuration of the visualizer
     * @param  dataset visualizer input, containing a dataset and optional visualization presets
     * @return         Promise that resolves to a [[DefaultVisualizer]]
     */
    public static load(config: Config, dataset: Dataset): Promise<DefaultVisualizer> {
        return new Promise((resolve, _) => {
            const visualizer = new DefaultVisualizer(config, dataset);
            resolve(visualizer);
        });
    }

    public map: PropertiesMap;
    public info: EnvironmentInfo;
    public meta: MetadataPanel;
    public structure: ViewersGrid;

    private _indexer: EnvironmentIndexer;
    // Stores raw input input so we can output it as JSON later
    private _dataset: Dataset;

    // the constructor is private because the main entry point is the static
    // `load` function
    private constructor(config: Config, dataset: Dataset) {
        validateConfig(config as unknown as JsObject);
        validateDataset(dataset as unknown as JsObject);

        this._dataset = dataset;

        const mode = (dataset.environments === undefined) ? 'structure' : 'atom';
        this._indexer = new EnvironmentIndexer(mode, dataset.structures, dataset.environments);

        this.meta = new MetadataPanel(config.meta, dataset.meta);

        // Structure viewer setup
        this.structure = new ViewersGrid(
            config.structure,
            config.j2sPath,
            this._indexer,
            dataset.structures,
            dataset.environments,
        );

        if (config.loadStructure !== undefined) {
            this.structure.loadStructure = config.loadStructure;
        }

        this.structure.activeChanged = (guid, indexes) => {
            this.map.setActive(guid);
            this.info.show(indexes);
        };

        this.structure.onselect = (indexes) => {
            this.map.select(indexes);
            this.info.show(indexes);
        };

        this.structure.onremove = (guid) => {
            this.map.removeMarker(guid);
        };

        this.structure.oncreate = (guid, color, indexes) => {
            this.map.addMarker(guid, color, indexes);
            this.info.show(indexes);
        };

        // map setup
        this.map = new PropertiesMap(config.map, this._indexer, dataset.properties);

        this.map.onselect = (indexes) => {
            this.info.show(indexes);
            this.structure.show(indexes);
        };

        this.map.activeChanged = (guid, indexes) => {
            this.info.show(indexes);
            this.structure.setActive(guid);
        };

        // information table & slider setup
        this.info = new EnvironmentInfo(config.info, dataset.properties, this._indexer);
        this.info.onchange = (indexes) => {
            this.map.select(indexes);
            this.structure.show(indexes);
        };
        this.info.startStructurePlayback = (advance) => this.structure.structurePlayback(advance);
        this.info.startAtomPlayback = (advance) => this.structure.atomPlayback(advance);

        const initial = {environment: 0, structure: 0, atom: 0};
        this.map.addMarker(this.structure.active, getNextColor([]), initial);
        this.structure.show(initial);
        this.info.show(initial);
    }

    /**
     * Removes all the chemiscope widgets from the DOM
     */
    public remove(): void {
        this.map.remove();
        this.meta.remove();
        this.info.remove();
        this.structure.remove();
    }

    /**
     * Dumps the dataset and settings as a JSON string
     */
    public dump(): string {
        const dumpdata = { ...this._dataset, ...{ presets: {map: this.map.dumpPresets()} } };
        return JSON.stringify(dumpdata);
    }
}

declare const CHEMISCOPE_GIT_VERSION: string;
function version(): string {
    return CHEMISCOPE_GIT_VERSION;
}

export {
    addWarningHandler,
    version,
    MetadataPanel,
    PropertiesMap,
    ViewersGrid,
    EnvironmentInfo,
    EnvironmentIndexer,
    DefaultVisualizer,
};
