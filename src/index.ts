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

import assert from 'assert';

import { EnvironmentIndexer, Indexes } from './indexer';
import { EnvironmentInfo } from './info';
import { PropertiesMap } from './map';
import { MetadataPanel } from './metadata';
import { SavedSettings } from './options';
import { ViewersGrid } from './structure';

import { Dataset, JsObject, Structure, validateDataset } from './dataset';
import { GUID, addWarningHandler, getNextColor } from './utils';

require('./static/chemiscope.css');

/**
 * Configuration for the [[DefaultVisualizer]]
 */
export interface Config {
    /** Id of the DOM element to use for the [[MetadataPanel|metadata display]] */
    meta: string | HTMLElement;
    /** Id of the DOM element to use for the [[PropertiesMap|properties map]] */
    map: string | HTMLElement;
    /** Id of the DOM element to use for the [[EnvironmentInfo|environment information]] */
    info: string | HTMLElement;
    /** Id of the DOM element to use for the [[ViewersGrid|structure viewer]] */
    structure: string | HTMLElement;
    /** Settings for the map & structure viewer */
    settings?: Partial<Settings>;
    /** Custom structure loading callback, used to set [[ViewersGrid.loadStructure]] */
    loadStructure?: (index: number, structure: unknown) => Structure;
}

export interface Settings {
    map: SavedSettings;
    structure: SavedSettings[];
    pinned: number[];
}

/** @hidden
 * Check if `o` contains all the expected fields to be a [[Config]].
 */
function validateConfig(o: JsObject) {
    if (typeof o !== 'object') {
        throw Error('the configuration must be a JavaScript object');
    }

    if (!('meta' in o && (typeof o.meta === 'string' || o.meta instanceof HTMLElement))) {
        throw Error('missing "meta" key in chemiscope configuration');
    }

    if (!('map' in o && (typeof o.map === 'string' || o.map instanceof HTMLElement))) {
        throw Error('missing "map" key in chemiscope configuration');
    }

    if (!('info' in o && (typeof o.info === 'string' || o.info instanceof HTMLElement))) {
        throw Error('missing "info" key in chemiscope configuration');
    }

    if (
        !(
            'structure' in o &&
            (typeof o.structure === 'string' || o.structure instanceof HTMLElement)
        )
    ) {
        throw Error('missing "structure" key in chemiscope configuration');
    }

    if ('settings' in o) {
        if (typeof o.settings !== 'object' || o.settings === null) {
            throw Error('"settings" must be an object in chemiscope configuration');
        }

        validateSettings(o.settings as JsObject);
    }

    // from underscore.js / https://stackoverflow.com/a/6000016
    // eslint-disable-next-line @typescript-eslint/no-explicit-any
    const isFunction = (obj: any) => {
        // eslint-disable-next-line @typescript-eslint/no-unsafe-member-access
        return !!(obj && obj.constructor && obj.call && obj.apply);
    };

    if ('loadStructure' in o && o.loadStructure !== undefined && !isFunction(o.loadStructure)) {
        throw Error('"loadStructure" should be a function in chemiscope config');
    }
}

function validateSettings(settings: JsObject) {
    const checkObject = (key: string, o: unknown) => {
        if (typeof o !== 'object' || o === null) {
            throw Error(`"settings.${key}" must be an object`);
        }
    };

    const checkArray = (key: string, o: unknown) => {
        if (!(Array.isArray(o) && o.length > 0)) {
            throw Error(`"settings.${key}" must be an array containing at least one element`);
        }
    };

    for (const key in settings) {
        if (key === 'map') {
            checkObject(key, settings.map);
        } else if (key === 'structure') {
            checkArray(key, settings.structure);
            for (const value of settings.structure as unknown[]) {
                checkObject('structure entry', value);
            }
        } else if (key === 'pinned') {
            checkArray(key, settings.pinned);
            for (const value of settings.pinned as unknown[]) {
                if (!(Number.isInteger(value) && (value as number) >= 0)) {
                    throw Error('"settings.pinned" must be an array of number');
                }
            }
        } else {
            throw Error(`invalid key "${key}" in settings`);
        }
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
     * @param  dataset visualizer input, containing a dataset and optional visualization settings
     * @return         Promise that resolves to a [[DefaultVisualizer]]
     */
    public static load(config: Config, dataset: Dataset): Promise<DefaultVisualizer> {
        return new Promise((resolve) => {
            const visualizer = new DefaultVisualizer(config, dataset);
            resolve(visualizer);
        });
    }

    public map: PropertiesMap;
    public info: EnvironmentInfo;
    public meta: MetadataPanel;
    public structure: ViewersGrid;

    private _indexer: EnvironmentIndexer;
    // Stores raw input input so we can give it back later
    private _dataset: Dataset;
    // Keep the list of pinned environments around to be able to apply settings
    private _pinned: GUID[];

    // the constructor is private because the main entry point is the static
    // `load` function
    private constructor(config: Config, dataset: Dataset) {
        validateConfig(config as unknown as JsObject);
        validateDataset(dataset as unknown as JsObject);

        this._dataset = dataset;
        this._pinned = [];

        const mode = dataset.environments === undefined ? 'structure' : 'atom';
        this._indexer = new EnvironmentIndexer(mode, dataset.structures, dataset.environments);

        this.meta = new MetadataPanel(config.meta, dataset.meta);

        // Structure viewer setup
        this.structure = new ViewersGrid(
            config.structure,
            this._indexer,
            dataset.structures,
            dataset.environments
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

            // remove the guid from this._pinned
            const index = this._pinned.indexOf(guid);
            assert(index > -1);
            this._pinned.splice(index, 1);
        };

        this.structure.oncreate = (guid, color, indexes) => {
            this.map.addMarker(guid, color, indexes);
            this.info.show(indexes);
            this._pinned.push(guid);
        };

        // map setup
        this.map = new PropertiesMap(
            { element: config.map, settings: getMapSettings(config.settings) },
            this._indexer,
            dataset.properties
        );

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

        let initial: Indexes = { environment: 0, structure: 0, atom: 0 };
        if (config.settings && config.settings.pinned) {
            initial = this._indexer.from_environment(config.settings.pinned[0]);
        }

        const firstGUID = this.structure.active;
        this.map.addMarker(firstGUID, getNextColor([]), initial);
        this._pinned.push(firstGUID);
        this.structure.show(initial);
        this.info.show(initial);

        // setup additional pinned values from the settings if needed
        if (config.settings !== undefined) {
            delete config.settings.map;
            this.applySettings(config.settings);
        }
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
     * Get the current values of settings for all panels in the visualizer
     *
     * @return the viewers settings, suitable to be used with [[applySettings]]
     */
    public saveSettings(): Settings {
        return {
            map: this.map.saveSettings(),
            pinned: this.structure.pinned().map((value) => value.environment),
            structure: this.structure.saveSettings(),
        };
    }

    /**
     * Apply the given settings to all panels in the visualizer
     *
     * @param settings settings for all panels
     */
    public applySettings(settings: Partial<Settings>): void {
        validateSettings(settings);

        if (settings.map !== undefined) {
            this.map.applySettings(settings.map);
        }

        if (settings.pinned !== undefined) {
            // remove all viewers except from the first one and start fresh
            for (const guid of this._pinned.slice(1)) {
                this.map.removeMarker(guid);
                this.structure.removeViewer(guid);
            }
            this._pinned = [this._pinned[0]];

            // Change the first viewer/marker
            assert(settings.pinned.length > 0);
            const indexes = this._indexer.from_environment(settings.pinned[0]);
            this.map.select(indexes);
            this.info.show(indexes);
            this.structure.show(indexes);

            // Create additional viewers as needed
            for (const environment of settings.pinned.slice(1)) {
                const [guid, color] = this.structure.addViewer();
                if (guid === undefined) {
                    throw Error("too many environments in 'pinned' setting");
                }
                const indexes = this._indexer.from_environment(environment);

                this.map.addMarker(guid, color, indexes);
                this.map.setActive(guid);
                this._pinned.push(guid);
                this.info.show(indexes);
                this.structure.setActive(guid);
                this.structure.show(indexes);
            }
        }

        if (settings.structure !== undefined) {
            this.structure.applySettings(settings.structure);
        }
    }

    /**
     * Get the dataset used to create the current visualization
     *
     * If the dataset is using user-specified structures and a loading callback
     * [[Config.loadStructure]]; you can request all structure to be fully
     * resolved and placed inside the dataset.
     *
     * @param  getStructures should all [[UserStructure]] resolved and placed
     *                       inside the dataset?
     * @return the dataset currently visualized
     */
    public dataset(getStructures: boolean = false): Dataset {
        // preserve NaN values in the copy
        const copy = JSON.parse(
            JSON.stringify(this._dataset, (_, value) => {
                // eslint-disable-next-line @typescript-eslint/no-unsafe-return
                return typeof value === 'number' && isNaN(value) ? '***NaN***' : value;
            }),
            // eslint-disable-next-line @typescript-eslint/no-unsafe-return
            (_, value) => (value === '***NaN***' ? NaN : value)
        ) as Dataset;

        if (getStructures) {
            copy.structures = [] as Structure[];
            for (let i = 0; i < this._dataset.structures.length; i++) {
                copy.structures.push(this.structure.loadStructure(i, this._dataset.structures[i]));
            }
        }
        return copy;
    }
}

declare const CHEMISCOPE_GIT_VERSION: string;
function version(): string {
    return CHEMISCOPE_GIT_VERSION;
}

function getMapSettings(settings: Partial<Settings> | undefined): SavedSettings {
    if (settings === undefined) {
        return {};
    }
    if ('map' in settings) {
        const map = settings.map;
        if (typeof map === 'object') {
            if (map === null) {
                throw Error('invalid settings for map, should not be null');
            } else {
                return map;
            }
        } else {
            throw Error(`invalid type '${typeof map}' for map, should be an object`);
        }
    } else {
        return {};
    }
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
