import 'construct-style-sheets-polyfill';
import assert from 'assert';

import { DisplayMode, EnvironmentIndexer, Indexes } from './indexer';
import { EnvironmentInfo } from './info';
import { PropertiesMap } from './map';
import { MetadataPanel } from './metadata';
import { LoadOptions, MoleculeViewer, StructureOptions, ViewersGrid } from './structure';
import { DisplayModeToggle } from './mode-toggle';

import {
    Dataset,
    Environment,
    Metadata,
    Parameter,
    Property,
    Settings,
    Structure,
    Target,
    UserStructure,
} from './dataset';
import { JsObject, getMode, validateDataset } from './dataset';
import {
    GUID,
    PositioningCallback,
    WarningHandler,
    addWarningHandler,
    getNextColor,
    sendWarning,
} from './utils';
import {
    ArrowParameters,
    CustomShapeParameters,
    EllipsoidParameters,
    ShapeParameters,
    SphereParameters,
} from './structure/shapes';

import './static/chemiscope.css';

/**
 * Welcome to Chemiscope's API documentation
 *
 * The default visualization is organized around four panels: the
 * {@link MetadataPanel|metadata} panel; the {@link PropertiesMap|map} (a scatter plot of
 * properties), the {@link ViewersGrid|structure viewer}, and the general
 * dataset {@link EnvironmentInfo|information}. Each one of these is defined in a
 * separate module.
 *
 * Other organization of the visualization are possible by using the classes
 * responsible for each sub-panel separately, instead of using the
 * {@link DefaultVisualizer}. In this case, developers should make sure to finish
 * the plumbing by setting the right callbacks on each element used.
 *
 * @packageDocumentation
 * @module main
 * @preferred
 */

/**
 * Configuration for the {@link DefaultVisualizer}
 */
export interface DefaultConfig {
    /** Id of the DOM element to use to display metadata */
    meta: string | HTMLElement;
    /** Id of the DOM element to use for the properties map */
    map: string | HTMLElement;
    /** Id of the DOM element to use for the environment information */
    info: string | HTMLElement;
    /** Id of the DOM element to use for the structure viewer grid */
    structure: string | HTMLElement;
    /** Custom structure loading callback, used to set {@link ViewersGrid.loadStructure} */
    loadStructure?: (index: number, structure: unknown) => Structure;
    /** Maximum number of structure viewers allowed in {@link ViewersGrid} */
    maxStructureViewers?: number;
}

/** @hidden
 * Check if `o` contains all valid config the expected fields to be a {@link DefaultConfig}.
 */
function validateConfig(o: JsObject, requiredIds: string[]) {
    if (typeof o !== 'object') {
        throw Error('the configuration must be a JavaScript object');
    }

    for (const id of requiredIds) {
        if (!(id in o && (typeof o[id] === 'string' || o[id] instanceof HTMLElement))) {
            throw Error(`missing "${id}" key in chemiscope configuration`);
        }
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

    if ('maxStructureViewers' in o && typeof o.maxStructureViewers !== 'number') {
        throw Error('"maxStructureViewers" should be a number in chemiscope config');
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
        } else if (key === 'mode') {
            if (!['atom', 'structure'].includes(settings.mode as string)) {
                throw Error('"settings.mode" should be either "atom" or "structure"');
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
     * @return         Promise that resolves to a {@link DefaultVisualizer}
     */
    public static load(config: DefaultConfig, dataset: Dataset): Promise<DefaultVisualizer> {
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
    // Display mode of the widgets, structure par default
    private _mode: DisplayMode = 'structure';
    // Stores raw input input so we can give it back later
    private _dataset: Dataset;
    // Keep the list of pinned environments around to be able to apply settings
    private _pinned: GUID[];
    // Controls the display mode (atom or structure)
    private _toggle: DisplayModeToggle | undefined;

    // the constructor is private because the main entry point is the static
    // `load` function
    private constructor(config: DefaultConfig, dataset: Dataset) {
        validateConfig(config as unknown as JsObject, ['meta', 'map', 'info', 'structure']);
        validateDataset(dataset as unknown as JsObject);

        this._dataset = dataset;
        this._pinned = [];

        this._mode = getMode(dataset);
        this._indexer = new EnvironmentIndexer(dataset.structures, dataset.environments);

        this.meta = new MetadataPanel(config.meta, dataset.meta);

        // Structure viewer setup
        this.structure = new ViewersGrid(
            config.structure,
            this._indexer,
            dataset.structures,
            this._mode,
            dataset.properties,
            dataset.environments,
            config.maxStructureViewers
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
            config.map,
            getMapSettings(dataset.settings),
            this._indexer,
            this._mode,
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

        // toggle setup

        // Check if toggle should be visible
        const newMode = this._mode === 'atom' ? 'structure' : 'atom';
        const modeProps = Object.values(dataset.properties).filter((p) => p.target === newMode);
        if (dataset.environments !== undefined && modeProps.length > 1) {
            // Initiate toggle
            this._toggle = new DisplayModeToggle(config.map, this._mode);
            this._toggle.onchange = (mode: DisplayMode) => this._switchMode(mode);
        }

        // information table & slider setup
        this.info = new EnvironmentInfo(
            config.info,
            dataset.properties,
            this._indexer,
            this._mode,
            dataset.parameters
        );
        this.info.onchange = (indexes) => {
            this.map.select(indexes);
            this.structure.show(indexes);
        };

        this.structure.delayChanged = (delay) => {
            this.info.playbackDelay = delay;
        };

        let initial: Indexes = { environment: 0, structure: 0, atom: 0 };

        // if we have sparse environments, make sure to use the first
        // environment actually part of the dataset
        if (dataset.environments !== undefined) {
            initial = this._indexer.fromEnvironment(0, this._mode);
        }

        if (dataset.settings && dataset.settings.pinned) {
            if (
                !Array.isArray(dataset.settings.pinned) ||
                typeof dataset.settings.pinned[0] !== 'number'
            ) {
                throw Error('settings.pinned must be an array of numbers');
            }
            initial = this._indexer.fromEnvironment(dataset.settings.pinned[0], this._mode);
        }

        const firstGUID = this.structure.active;
        this.map.addMarker(firstGUID, getNextColor([]), initial);
        this._pinned.push(firstGUID);
        this.structure.show(initial);
        this.info.show(initial);

        // setup additional pinned values from the settings if needed
        if (dataset.settings !== undefined) {
            delete dataset.settings.map;
            this.applySettings(dataset.settings);
        }
    }

    /**
     * Update all elements to the related display mode
     * @param mode visualisation display mode (atom or structure)
     */
    private _switchMode(mode: DisplayMode) {
        // Check if mode actually new
        if (this._mode !== mode) {
            assert(this._toggle !== undefined);

            // Show loader
            this._toggle.loader(true);

            // Update the current mode
            this._mode = mode;

            // Use setTimeout to ensure the loader is shown before starting async operations
            // eslint-disable-next-line @typescript-eslint/no-misused-promises
            setTimeout(async () => {
                try {
                    // Proceed with EnvironmentInfo
                    this.info.switchMode(this._mode);

                    // Proceed with PropertiesMap
                    await this.map.switchMode(this._mode);

                    // Proceed with ViewersGrid
                    await this.structure.switchMode(this._mode);
                } finally {
                    // Hide loader
                    this._toggle?.loader(false);
                }
            }, 0);
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
        this._toggle?.remove();
    }

    /**
     * Get the current values of settings for all panels in the visualizer
     *
     * @return the viewers settings, suitable to be used with {@link applySettings}
     */
    public saveSettings(): Settings {
        return {
            mode: this._mode,
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
            this.map.applySettings(settings.map as Settings);
        }

        if (settings.pinned !== undefined) {
            if (!Array.isArray(settings.pinned)) {
                throw Error('settings.pinned must be an array');
            }

            // remove all viewers except from the first one and start fresh
            for (const guid of this._pinned.slice(1)) {
                this.map.removeMarker(guid);
                this.structure.removeViewer(guid);
            }
            this._pinned = [this._pinned[0]];

            // Change the first viewer/marker
            assert(settings.pinned.length > 0);
            if (typeof settings.pinned[0] !== 'number') {
                throw Error('settings.pinned must be an array of numbers');
            }

            const indexes = this._indexer.fromEnvironment(settings.pinned[0], this._mode);
            this.map.select(indexes);
            this.info.show(indexes);
            this.structure.show(indexes);

            // Create additional viewers as needed
            for (const environment of settings.pinned.slice(1)) {
                if (typeof environment !== 'number') {
                    throw Error('settings.pinned must be an array of numbers');
                }

                const { guid, color } = this.structure.addViewer();
                if (guid === undefined) {
                    throw Error("too many environments in 'pinned' setting");
                }
                const indexes = this._indexer.fromEnvironment(environment, this._mode);

                this.map.addMarker(guid, color, indexes);
                this.map.setActive(guid);
                this._pinned.push(guid);
                this.info.show(indexes);
                this.structure.setActive(guid);
                this.structure.show(indexes);
            }
        }

        if (settings.structure !== undefined) {
            this.structure.applySettings(settings.structure as Settings[]);
        }
    }

    /**
     * Add the given `callback` to be called whenever a setting changes. The
     * callback will be given the path to the settings as a list of keys; and
     * the new value of the setting.
     *
     * There is currently no way to remove a callback.
     */
    public onSettingChange(callback: (keys: string[], value: unknown) => void): void {
        this.map.onSettingChange((keys, value) => {
            keys = JSON.parse(JSON.stringify(keys)) as string[];
            keys.unshift('map');

            callback(keys, value);
        });

        this.structure.onSettingChange((keys, value) => {
            keys = JSON.parse(JSON.stringify(keys)) as string[];
            keys.unshift('structure');

            callback(keys, value);
        });
    }

    /**
     * Get the dataset used to create the current visualization
     *
     * If the dataset is using user-specified structures and a loading callback
     * {@link DefaultConfig.loadStructure}; you can request all structure to be fully
     * resolved and placed inside the dataset.
     *
     * @param  getStructures should all {@link UserStructure} resolved and placed
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

/**
 * Configuration for the {@link StructureVisualizer}
 */
export interface StructureConfig {
    /** Id of the DOM element to use to display metadata */
    meta: string | HTMLElement;
    /** Id of the DOM element to use for the environment information */
    info: string | HTMLElement;
    /** Id of the DOM element to use for the structure viewer grid */
    structure: string | HTMLElement;
    /** Custom structure loading callback, used to set {@link ViewersGrid.loadStructure} */
    loadStructure?: (index: number, structure: unknown) => Structure;
}

/**
 * A structure-only chemiscope visualizer: two panels (map,
 * info) updating one another when the user interact with any of them.
 */
class StructureVisualizer {
    /**
     * Load a dataset and create a visualizer.
     *
     * This function returns a `Promise<StructureVisualizer>` to prevent blocking
     * the browser while everything is loading.
     *
     * @param  config  configuration of the visualizer
     * @param  dataset visualizer input, containing a dataset and optional visualization settings
     * @return         Promise that resolves to a {@link StructureVisualizer}
     */
    public static load(config: StructureConfig, dataset: Dataset): Promise<StructureVisualizer> {
        return new Promise((resolve) => {
            const visualizer = new StructureVisualizer(config, dataset);
            resolve(visualizer);
        });
    }

    public info: EnvironmentInfo;
    public meta: MetadataPanel;
    public structure: ViewersGrid;

    private _indexer: EnvironmentIndexer;
    private _mode: DisplayMode;

    // the constructor is private because the main entry point is the static
    // `load` function
    private constructor(config: StructureConfig, dataset: Dataset) {
        validateConfig(config as unknown as JsObject, ['meta', 'info', 'structure']);
        validateDataset(dataset as unknown as JsObject);

        this._mode = getMode(dataset);
        this._indexer = new EnvironmentIndexer(dataset.structures, dataset.environments);

        this.meta = new MetadataPanel(config.meta, dataset.meta);

        // Structure viewer setup
        this.structure = new ViewersGrid(
            config.structure,
            this._indexer,
            dataset.structures,
            this._mode,
            dataset.properties,
            dataset.environments,
            1 // only allow one structure
        );

        if (config.loadStructure !== undefined) {
            this.structure.loadStructure = config.loadStructure;
        }

        this.structure.activeChanged = (guid, indexes) => {
            this.info.show(indexes);
        };

        this.structure.onselect = (indexes) => {
            this.info.show(indexes);
        };

        this.structure.oncreate = (guid, color, indexes) => {
            this.info.show(indexes);
        };

        // information table & slider setup
        this.info = new EnvironmentInfo(config.info, dataset.properties, this._indexer, this._mode);
        this.info.onchange = (indexes) => {
            this.structure.show(indexes);
        };

        this.structure.delayChanged = (delay) => {
            this.info.playbackDelay = delay;
        };

        let initial: Indexes = { environment: 0, structure: 0, atom: 0 };
        // if we have sparse environments, make sure to use the first
        // environment actually part of the dataset
        if (dataset.environments !== undefined) {
            initial = this._indexer.fromEnvironment(0, this._mode);
        }

        if (dataset.settings && dataset.settings.pinned) {
            if (
                !Array.isArray(dataset.settings.pinned) ||
                typeof dataset.settings.pinned[0] !== 'number'
            ) {
                throw Error('settings.pinned must be an array of numbers');
            }
            initial = this._indexer.fromEnvironment(dataset.settings.pinned[0], this._mode);
        }

        this.structure.show(initial);
        this.info.show(initial);

        // apply settings if needed
        if (dataset.settings !== undefined) {
            this.applySettings(dataset.settings);
        }
    }

    /**
     * Removes all the chemiscope widgets from the DOM
     */
    public remove(): void {
        this.meta.remove();
        this.info.remove();
        this.structure.remove();
    }

    /**
     * Apply the given `settings` to the structure panels in the visualizer
     */
    public applySettings(settings: Partial<Settings>): void {
        validateSettings(settings);

        if (settings.structure !== undefined) {
            this.structure.applySettings(settings.structure as Settings[]);
        }
    }

    /**
     * Get the current values of settings for all panels in the visualizer
     */
    public saveSettings(): Settings {
        return {
            mode: this._mode,
            pinned: this.structure.pinned().map((value) => value.environment),
            structure: this.structure.saveSettings(),
        };
    }

    /**
     * Add the given `callback` to be called whenever a setting changes. The
     * callback will be given the path to the settings as a list of keys; and
     * the new value of the setting.
     *
     * There is currently no way to remove a callback.
     */
    public onSettingChange(callback: (keys: string[], value: unknown) => void): void {
        this.structure.onSettingChange((keys, value) => {
            keys = JSON.parse(JSON.stringify(keys)) as string[];
            keys.unshift('structure');

            callback(keys, value);
        });
    }
}

/**
 * Configuration for the {@link MapVisualizer}
 */
export interface MapConfig {
    /** Id of the DOM element to use to display metadata */
    meta: string | HTMLElement;
    /** Id of the DOM element to use for the properties map */
    map: string | HTMLElement;
    /** Id of the DOM element to use for the environment information */
    info: string | HTMLElement;
}

/**
 * A map-only visualizer state of chemiscope
 */
class MapVisualizer {
    /**
     * Load a dataset and create a visualizer.
     *
     * This function returns a `Promise<MapVisualizer>` to prevent blocking
     * the browser while everything is loading.
     *
     * @param  config  configuration of the visualizer
     * @param  dataset visualizer input, containing a dataset and optional visualization settings
     * @return         Promise that resolves to a {@link MapVisualizer}
     */
    public static load(config: MapConfig, dataset: Dataset): Promise<MapVisualizer> {
        return new Promise((resolve) => {
            const visualizer = new MapVisualizer(config, dataset);
            resolve(visualizer);
        });
    }

    public info: EnvironmentInfo;
    public map: PropertiesMap;
    public meta: MetadataPanel;

    private _indexer: EnvironmentIndexer;

    // Display mode of the map widget, structure par default
    private _mode: DisplayMode = 'structure';

    // the constructor is private because the main entry point is the static
    // `load` function
    private constructor(config: MapConfig, dataset: Dataset) {
        validateConfig(config as unknown as JsObject, ['meta', 'map', 'info']);
        validateDataset(dataset as unknown as JsObject);

        // We need to create pseudo structures with 0 atoms so that the indexer
        // knows how to associate structure id & point id on the map. The
        // structure id is used by the info panel
        let n_structure;
        for (const key in dataset.properties) {
            const property = dataset.properties[key];
            if (property.target === 'atom') {
                sendWarning('unsupported per-atom property in a map-only viewer');
            } else {
                assert(property.target === 'structure');
                n_structure = property.values.length;
            }
        }

        this._indexer = new EnvironmentIndexer(
            // pseudo-structures to get the EnvironmentIndexer working even
            // though we don't have any structure data
            new Array(n_structure).fill({ size: 1, data: 0 }) as UserStructure[]
        );

        this.meta = new MetadataPanel(config.meta, dataset.meta);

        // map setup
        this.map = new PropertiesMap(
            config.map,
            getMapSettings(dataset.settings),
            this._indexer,
            this._mode,
            dataset.properties
        );

        this.map.onselect = (indexes) => {
            this.info.show(indexes);
        };

        this.map.activeChanged = (guid, indexes) => {
            this.info.show(indexes);
        };

        // information table & slider setup
        this.info = new EnvironmentInfo(config.info, dataset.properties, this._indexer, this._mode);
        this.info.onchange = (indexes) => {
            this.map.select(indexes);
        };

        let initial: Indexes = { environment: 0, structure: 0, atom: 0 };

        // if we have sparse environments, make sure to use the first
        // environment actually part of the dataset
        if (dataset.environments !== undefined) {
            initial = this._indexer.fromEnvironment(0, this._mode);
        }

        if (dataset.settings && dataset.settings.pinned) {
            if (
                !Array.isArray(dataset.settings.pinned) ||
                typeof dataset.settings.pinned[0] !== 'number'
            ) {
                throw Error('settings.pinned must be an array of numbers');
            }
            initial = this._indexer.fromEnvironment(dataset.settings.pinned[0], this._mode);
        }

        this.map.addMarker('map-0' as GUID, 'red', initial);
        this.info.show(initial);
    }

    /**
     * Removes all the chemiscope widgets from the DOM
     */
    public remove(): void {
        this.meta.remove();
        this.info.remove();
        this.map.remove();
    }

    /**
     * Get the current values of settings for all panels in the visualizer
     */
    public saveSettings(): Settings {
        return {
            mode: this._mode,
            map: this.map.saveSettings(),
        };
    }

    /**
     * Apply the given `settings` to the structure panels in the visualizer
     */
    public applySettings(settings: Partial<Settings>): void {
        validateSettings(settings);

        if (settings.map !== undefined) {
            this.map.applySettings(getMapSettings(settings));
        }
    }

    /**
     * Add the given `callback` to be called whenever a setting changes. The
     * callback will be given the path to the settings as a list of keys; and
     * the new value of the setting.
     *
     * There is currently no way to remove a callback.
     */
    public onSettingChange(callback: (keys: string[], value: unknown) => void): void {
        this.map.onSettingChange((keys, value) => {
            keys = JSON.parse(JSON.stringify(keys)) as string[];
            keys.unshift('map');

            callback(keys, value);
        });
    }
}

declare const CHEMISCOPE_GIT_VERSION: string;
/** Get the version of chemiscope as a string */
function version(): string {
    return CHEMISCOPE_GIT_VERSION;
}

function getMapSettings(settings: Partial<Settings> | undefined): Settings {
    if (settings === undefined) {
        return {};
    }

    if ('map' in settings) {
        const map = settings.map as Settings;
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
    // free functions
    addWarningHandler,
    version,
    // dataset definitions
    Dataset,
    Metadata,
    Parameter,
    Property,
    Target,
    Structure,
    UserStructure,
    Environment,
    Settings,
    ShapeParameters,
    SphereParameters,
    EllipsoidParameters,
    ArrowParameters,
    CustomShapeParameters,
    // different panels
    MetadataPanel,
    PropertiesMap,
    ViewersGrid,
    MoleculeViewer,
    StructureOptions,
    LoadOptions,
    EnvironmentInfo,
    EnvironmentIndexer,
    // different visualizers built with the panels
    DefaultVisualizer,
    StructureVisualizer,
    MapVisualizer,
    // miscellaneous helpers
    DisplayMode,
    GUID,
    PositioningCallback,
    Indexes,
    WarningHandler,
};
