/**
 * @packageDocumentation
 * @module structure
 */

import assert from 'assert';

import {
    Environment,
    JsObject,
    Property,
    Settings,
    Structure,
    UserStructure,
    checkStructure,
} from '../dataset';

import { DisplayTarget, EnvironmentIndexer, Indexes } from '../indexer';
import * as styles from '../styles';
import { GUID, PositioningCallback, Warnings, getElement } from '../utils';
import { enumerate, generateGUID, getByID, getFirstKey, getNextColor } from '../utils';

import { MoleculeViewer } from './viewer';

import { filter } from '../info/info';

import CLOSE_SVG from '../static/close.svg';
import DUPLICATE_SVG from '../static/duplicate.svg';
import PNG_SVG from '../static/download-png.svg';
import {
    LoadOptions,
    NumberedEnvironment,
    ViewerGridData,
    downloadURI,
    groupByStructure,
} from './utils';

/**
 * The {@link ViewersGrid} class displays a grid of molecule or a crystal viewers
 * in 3D using {@link MoleculeViewer} for rendering.
 */
export class ViewersGrid {
    /** Callback used when the user select an environment */
    public onselect: (indexes: Indexes) => void;
    /** Callback fired when a viewer is removed from the grid */
    public onremove: (guid: GUID) => void;
    /**
     * Callback fired when a new viewer is created
     *
     * @param guid GUID of the new viewer
     * @param color GUID of the marker indicating the new viewer
     * @param indexes environment showed in the new viewer
     */
    public oncreate: (guid: GUID, color: string, indexes: Indexes) => void;
    /**
     * Callback fired when the active viewer is changed
     *
     * @param guid GUID of the new active viewer
     * @param indexes environment showed in the new active viewer
     */
    public activeChanged: (guid: GUID, indexes: Indexes) => void;
    /**
     * Callback fired when the playback delay for structure/atom playback changed
     *
     * @param delay the new delay value
     */
    public delayChanged: (delay: number) => void;

    /**
     * Callback used when a new structure should be loaded
     *
     * By default, this assumes that the loaded dataset contains {@link Structure}
     * directly, and returns the data from there. If the loaded dataset contains
     * {@link UserStructure} instead, this callback should be set to transform from
     * data in {@link UserStructure.data} to a {@link Structure}.
     *
     * The callback gets two parameter: the structure index (0-based); and the
     * full {@link UserStructure}.
     */
    public loadStructure: (index: number, structure: unknown) => Promise<Structure>;

    public warnings: Warnings;

    /// Shadow root for isolation
    private _shadow: ShadowRoot;
    /// Root element containing the full viewer grid
    private _root: HTMLElement;
    /// List of structures in the dataset
    private _structures: Structure[] | UserStructure[];
    /// List of per-atom properties in the dataset
    private _properties: Record<string, number[]>;
    /// Widget display target (structure or atom)
    private _target: DisplayTarget;
    /// Cached string representation of structures
    private _resolvedStructures: Structure[];
    /// Optional list of environments for each structure
    private _environments?: NumberedEnvironment[][];
    /// Optional list of environments for each structure to save it once we change the target
    private _calculatedEnvironments?: NumberedEnvironment[][];

    /// Maximum number of allowed structure viewers
    private _maxViewers: number;
    /// The indexer translating between environments indexes and structure/atom
    /// indexes
    private _indexer: EnvironmentIndexer;
    /// GUID of the currently active viewer
    private _active: GUID;
    /// Map of cell GUIDS to their color, viewer, and current structure
    private _cellsData: Map<GUID, ViewerGridData>;
    /// Callback used to override all grid viewers' positionSettingsModal
    private _positionSettingsModal?: PositioningCallback;
    /// Store the `onSettingChange` callbacks to be able to add them to new
    /// viewers in the grid
    private _onSettingChangeCallbacks: Array<(keys: string[], value: unknown) => void>;

    /**
     * Create a new {@link ViewersGrid} inside the HTML element with the given
     * `id`
     *
     * @param element      HTML element or string 'id' of the element where
     *                     viewer should live
     * @param indexer      {@link EnvironmentIndexer} used to translate indexes from
     *                     environments index to structure/atom indexes
     * @param structures   list of structure to display
     * @param target       widget display modtargete (structure or atom)
     * @param properties   list of properties
     * @param environments list of atom-centered environments in the structures,
     *                     used to highlight the selected environment
     * @param maxViewers   maximum number of allowed structure viewers
     */
    constructor(
        element: string | HTMLElement,
        indexer: EnvironmentIndexer,
        structures: Structure[] | UserStructure[],
        target: DisplayTarget,
        properties?: { [name: string]: Property },
        environments?: Environment[],
        maxViewers: number = 9,
        warnings: Warnings = new Warnings()
    ) {
        this._target = target;
        this._structures = structures;
        if (properties === undefined) {
            this._properties = {};
        } else {
            // Calculate properties
            const numberProperties = filter(properties, (p) =>
                Object.values(p.values).every((v) => typeof v === 'number')
            );
            const atomProperties = filter(numberProperties, (p) => p.target === 'atom');
            this._properties = Object.fromEntries(
                Object.entries(atomProperties).map(([key, value]) => [
                    key,
                    value.values as number[],
                ])
            );
        }

        this._resolvedStructures = new Array<Structure>(structures.length);
        this._indexer = indexer;

        // Group environments by structure and save the result for future use
        this._calculatedEnvironments = groupByStructure(this._structures, environments);

        // Set the environments to synchronize with the modal and grid
        this._environments = this._calculatedEnvironments;

        // Default structure loading function
        this.loadStructure = (_: unknown, s: unknown): Promise<Structure> => {
            const error = checkStructure(s as JsObject);
            if (error !== '') {
                return Promise.reject(
                    new Error(
                        'got custom data for this structure, but no custom loadStructure callback\n' +
                            `the object was ${JSON.stringify(s)}`
                    )
                );
            }

            return Promise.resolve(s as Structure);
        };

        this.warnings = warnings;

        // Initializes with 1 viewer upon opening.
        this._cellsData = new Map<GUID, ViewerGridData>();
        this.onselect = () => {};
        this.onremove = () => {};
        this.oncreate = () => {};
        this.activeChanged = () => {};
        this.delayChanged = () => {};
        this._onSettingChangeCallbacks = [];

        if (maxViewers > 9) {
            throw Error('chemiscope only supports up to 9 structure viewers in the grid');
        }
        this._maxViewers = maxViewers;

        // Attach a shadow DOM to the host element for isolation
        const containerElement = getElement(element);
        const hostElement = document.createElement('div');
        hostElement.style.setProperty('height', '100%');
        containerElement.appendChild(hostElement);
        this._shadow = hostElement.attachShadow({ mode: 'open' });
        this._shadow.adoptedStyleSheets = [styles.bootstrap, styles.chemiscope];

        // Create a root element inside the shadow DOM to append grid to it
        this._root = document.createElement('div');
        this._root.id = 'grid-root';
        this._root.className = 'chsp-structure-viewer-grid';
        this._shadow.appendChild(this._root);

        this._setupGrid(1);

        this._active = getFirstKey(this._cellsData);
        this.setActive(this._active);
    }

    // eslint-disable-next-line @typescript-eslint/no-unnecessary-type-parameters
    private _getById<T = HTMLElement>(id: string): T {
        return getByID(id, this._shadow);
    }

    /**
     * Get the current active viewer GUID
     */
    public get active(): GUID {
        return this._active;
    }

    /**
     * Get the current active viewer index (position in the grid)
     */
    public get activeIndex(): number {
        return this._getCurrentPositionInGrid(this._active);
    }

    /**
     * Set the current active viewer by index (position in the grid)
     */
    public set activeIndex(index: number) {
        let i = 0;
        for (const guid of this._cellsData.keys()) {
            if (i === index) {
                this.setActive(guid);
                return;
            }
            i++;
        }
        throw Error(`invalid viewer index ${index}`);
    }

    /**
     * Get the PNG data for the active viewer
     */
    public exportActivePNG(): string {
        const data = this._cellsData.get(this._active);
        assert(data !== undefined);
        return data.viewer.exportPNG();
    }

    /**
     * Get the current list of environments showed inside the different viewer
     */
    public pinned(): Indexes[] {
        const result = [];
        for (const data of this._cellsData.values()) {
            result.push(data.current);
        }
        return result;
    }

    /**
     * Show a new structure, as identified by `indexes` in the active viewer.
     *
     * This will switch to the structure at index `indexes.structure`, and if
     * environments where passed to the constructor and the current display
     * target is `'atom'`, highlight the atom-centered environment corresponding
     * to `indexes.atom`.
     *
     * @param  indexes         structure / atom pair to display
     */
    public async show(indexes: Indexes): Promise<void> {
        await this._showInViewer(this.active, indexes).catch((e: unknown) => {
            throw e;
        });
    }

    /**
     * Add a new empty viewer to the grid
     *
     * @return the GUID and color of the new viewer, the GUID is `undefined`
     *         if we already reached the viewer limit.
     */
    public addViewer(): { guid?: GUID; color: string } {
        const newGUIDs = this._setupGrid(this._cellsData.size + 1);
        if (newGUIDs.length === 0) {
            // no new viewer, probably because we already have the maximum
            return { guid: undefined, color: '' };
        }
        assert(newGUIDs.length === 1);
        const newGUID = newGUIDs[0];

        const newData = this._cellsData.get(newGUID);
        assert(newData !== undefined);

        for (const callback of this._onSettingChangeCallbacks) {
            newData.viewer.onSettingChange((keys, value) => {
                // add the index of the viewer in the current grid of viewers
                // at the front of the keys
                let i = 0;
                for (const currentGuid of this._cellsData.keys()) {
                    if (newGUID === currentGuid) {
                        break;
                    }
                    i += 1;
                }
                keys.unshift(i.toString());

                callback(keys, value);
            });
        }

        return { guid: newGUID, color: newData.color };
    }

    /**
     * Removes the viewer with the given `guid` from the viewer grid.
     *
     * @param guid GUID of the viewer to remove
     */
    public removeViewer(guid: GUID): void {
        // don't remove the last viewer
        assert(this._cellsData.size > 1);

        // If we removed the active marker, change the active one
        if (this._active === guid) {
            const newActive = getFirstKey(this._cellsData, guid);

            this.setActive(newActive);

            const data = this._cellsData.get(newActive);
            assert(data !== undefined);
            this.activeChanged(newActive, data.current);
        }

        // remove HTML inserted by the viewer
        const data = this._cellsData.get(guid);
        assert(data !== undefined);
        data.viewer.remove();

        // remove the cell containing the viewer
        const cell = this._getById(`gi-${guid}`);
        cell.remove();

        this._cellsData.delete(guid);
    }

    /**
     * Remove all HTML added by this {@link ViewersGrid} in the current document
     */
    public remove(): void {
        for (const data of this._cellsData.values()) {
            data.viewer.remove();
        }

        this._shadow.host.remove();
    }

    /**
     * Function to set the active viewer for communicating with the map
     * @param guid GUID of the viewer to be active
     */
    public setActive(guid: GUID): void {
        const changeClasses = (guid: string, toggle: boolean) => {
            // change tooltip text in the active marker
            const button = this._getById(`chsp-activate-${guid}`);
            button.classList.toggle('chsp-active-pin', toggle);
            assert(button.parentElement !== null);
            const tooltip = button.parentElement.querySelector('.chsp-tooltip');
            assert(tooltip !== null);
            tooltip.innerHTML = toggle ? 'Active viewer' : 'Choose as active';

            // change style of the cell border
            const cell = this._getById(`gi-${guid}`);
            cell.classList.toggle('chsp-structure-viewer-cell-active', toggle);
        };

        const current = this._cellsData.get(this._active);
        assert(current !== undefined);

        // remove active classes from the previous active viewer
        changeClasses(this._active, false);

        this._active = guid;

        const newData = this._cellsData.get(this._active);
        assert(newData !== undefined);

        // set the right initial value for playback delay
        this.delayChanged(newData.viewer._options.playbackDelay.value);

        changeClasses(this._active, true);
    }

    /**
     * Set a callback to get the initial positioning of the settings modal of
     * the viewers. The same callback is used for all viewers in the grid.
     *
     * The callback gets the current placement of the settings as a
     * [DOMRect](https://developer.mozilla.org/en-US/docs/Web/API/DOMRect), and
     * should return top and left positions in pixels, used with `position:
     * fixed`. The callback is called once, the first time the settings are
     * opened.
     */
    public set positionSettingsModal(value: PositioningCallback) {
        this._positionSettingsModal = value;
        for (const data of this._cellsData.values()) {
            data.viewer.positionSettingsModal = value;
        }
    }

    /**
     * Apply the given saved setting to all viewers in the grid.
     *
     * The settings must be in viewer creation order. If the list
     * of settings contains only one item, then it is applied to
     * the active viewer.
     *
     * @param settings settings for either the active, or all
     *      viewers in the grid
     */
    public applySettings(settings: Settings[]): void {
        if (settings.length === 0) {
            return;
        }

        // applies only to the active cell
        if (settings.length === 1) {
            const data = this._cellsData.get(this._active);
            if (data !== undefined) {
                data.viewer.applySettings(settings[0]);
            }
            return;
        }

        if (settings.length !== this._cellsData.size) {
            throw Error(
                `The number of viewers (${this._cellsData.size}) is different from the structure settings specified (${settings.length})`
            );
        }

        for (const [i, data] of enumerate(this._cellsData.values())) {
            data.viewer.applySettings(settings[i]);
        }
    }

    /**
     * Get the current values of settings for all viewers in the grid.
     *
     * The settings are given in viewer creation order.
     *
     * @return the settings in an array, suitable to be used with {@link applySettings}
     */
    public saveSettings(): Settings[] {
        const settings = [];
        for (const data of this._cellsData.values()) {
            settings.push(data.viewer.saveSettings());
        }
        return settings;
    }

    /**
     * Add the given `callback` to be called whenever a setting changes. The
     * callback will be given the path to the settings as a list of keys; and
     * the new value of the setting.
     *
     * There is currently no way to remove a callback.
     */
    public onSettingChange(callback: (keys: string[], value: unknown) => void): void {
        for (const guid of this._cellsData.keys()) {
            const data = this._cellsData.get(guid);
            assert(data !== undefined);

            data.viewer.onSettingChange((keys: string[], value: unknown) => {
                keys = JSON.parse(JSON.stringify(keys)) as string[];
                keys.unshift(this._getCurrentPositionInGrid(guid).toString());

                callback(keys, value);
            });
        }

        this._onSettingChangeCallbacks.push(callback);
    }

    /**
     * Change widget display target and adapt the element to the new target
     * @param target display target
     */
    public async switchTarget(target: DisplayTarget): Promise<void> {
        // Update widget display target
        this._target = target;

        // Set/reset the environements of grid and viewers based on the target
        this._updateEnvironments(this._target === 'atom');

        // Refresh each cell with the new target
        const refreshPromises = Array.from(this._cellsData.entries()).map(([guid, data]) =>
            this._refreshCell(guid, data, this._target === 'atom')
        );

        // Wait for all refreshes to complete
        await Promise.all(refreshPromises);
    }

    /**
     * Updates the environments based on the target
     * @param envView indicates whether the current dispay target per environments
     */
    private _updateEnvironments(envView: boolean): void {
        // Proceed with target changes by setting/resetting the environments
        this._environments = envView ? this._calculatedEnvironments : undefined;
        if (envView) {
            assert(this._environments !== undefined);
            for (const data of this._cellsData.values()) {
                data.viewer.environments = this._environments[data.current.structure];
            }
        }
    }

    /**
     * Refreshes the cell and updates its data based on the target
     *
     * @param guid id of the cell to refresh
     * @param data indexes associated with the viewer grid
     * @param envView indicates whether the current dispay target is per environments
     */
    private async _refreshCell(guid: GUID, data: ViewerGridData, envView: boolean): Promise<void> {
        // Set/remove atom from indexes based on target
        if (envView) {
            // Check if environment is valid before accessing indexer
            if (data.current.environment >= 0) {
                data.current.atom = this._indexer.fromEnvironment(
                    data.current.environment,
                    this._target
                ).atom;
            } else {
                data.current.atom = undefined;
            }
        } else {
            data.current.atom = undefined;
        }
        this._cellsData.set(guid, data);

        // Recreate stucture settings modal
        const propertyNames = this._properties ? Object.keys(this._properties) : [];
        data.viewer.refreshOptions(envView, propertyNames);

        // Load the viewer with the current indexes
        if (data.current.structure >= 0) {
            await this._loadViewer(data.viewer, data.current.structure, data.current.atom);
        }
        // Highlight the atom if needed. This also trigger a style update for the atoms
        data.viewer.highlight(envView ? data.current.atom : undefined);
    }

    /**
     * Add a new structure viewer to the grid as a copy of the viewer with the
     * `initial` GUID. The new structure viewer is set as the active one.
     *
     * @param  initial GUID of the viewer to duplicate
     * @return the GUID of the new viewer, or `undefined` if we already
     *         reached the viewer limit.
     */
    private async _duplicate(initial: GUID): Promise<ViewerGridData | undefined> {
        const data = this._cellsData.get(initial);
        assert(data !== undefined);

        const newGUID = this.addViewer().guid;
        if (newGUID === undefined) {
            return undefined;
        }
        const newData = this._cellsData.get(newGUID);
        assert(newData !== undefined);

        await this._showInViewer(newGUID, data.current);
        newData.viewer.applySettings(data.viewer.saveSettings());
        this.setActive(newGUID);

        return newData;
    }

    /**
     * Get the structure at the given index, potentially loading them using a
     * user provided loading function.
     *
     * @param  index index of the structure
     * @return       a Structure instance
     */
    private async _structure(index: number): Promise<Structure> {
        if (this._resolvedStructures[index] === undefined) {
            const structure = await this.loadStructure(index, this._structures[index]);

            const check = checkStructure(structure as unknown as JsObject);
            if (check !== '') {
                throw Error(
                    `got invalid object as structure: ${check}
    the object was ${JSON.stringify(structure)}`
                );
            }

            // copy the shapes from the structures, if any
            if ('shapes' in this._structures[index]) {
                structure.shapes = this._structures[index].shapes;
            }

            this._resolvedStructures[index] = structure;
        }
        return this._resolvedStructures[index];
    }

    /**
     * Get the values of the properties for all the atoms in the current structure
     *
     * @param structure index of the current structure
     * @returns
     */
    private _propertiesForStructure(
        structure: number
    ): Record<string, (number | undefined)[]> | undefined {
        const properties: Record<string, (number | undefined)[]> = {};

        // Use _calculatedEnvironments instead of _environments.
        // _environments is undefined when target is 'structure' (to disable environment UI),
        // but we still want to be able to access atomic properties for coloring.
        if (this._calculatedEnvironments !== undefined) {
            const environments = this._calculatedEnvironments[structure];
            for (const name in this._properties) {
                const values = this._properties[name];

                properties[name] = [];
                for (const environment of environments) {
                    if (environment !== undefined) {
                        properties[name].push(values[environment.index]);
                    } else {
                        properties[name].push(undefined);
                    }
                }
            }

            return properties;
        } else {
            return undefined;
        }
    }

    /**
     * Displays the specified structure and atom in the viewer
     *
     * @param guid id of the cell to show
     * @param indexes environment showed in the new viewer
     */
    private async _showInViewer(guid: GUID, indexes: Indexes): Promise<void> {
        // Get current data from the cell
        const data = this._cellsData.get(guid);
        assert(data !== undefined);

        const viewer = data.viewer;
        // If the structure has changed, load the new structure and atom into the viewer
        if (data.current.structure !== indexes.structure) {
            await this._loadViewer(viewer, indexes.structure, indexes.atom);
            data.current = indexes;
        }

        // Toggle the highlight for the atom in the viewer
        viewer.highlight(this._target === 'atom' ? indexes.atom : undefined);

        // Update the current indexes in the data
        data.current = indexes;
    }

    /**
     * Loads a structure and optional atom into the viewer
     *
     * @param viewer viewer instance where the structure will be loaded
     * @param structureIndex index of the structure to be loaded
     * @param atomIndex optional index of the atom to highlight
     */
    private async _loadViewer(
        viewer: MoleculeViewer,
        structureIndex: number,
        atomIndex?: number
    ): Promise<void> {
        const properties = this._propertiesForStructure(structureIndex);
        // Initialize load options with trajectory enabled
        const options: Partial<LoadOptions> = {
            target: this._target,
            trajectory: true,
        };

        // If environments are defined, add them to the options
        if (this._environments !== undefined) {
            options.environments = this._environments[structureIndex];

            // If the target is per atom, add the atom index to the highlight
            if (this._target === 'atom') {
                options.highlight = atomIndex;
            }
        }

        const structure = await this._structure(structureIndex);

        // Load the structure into the viewer
        await new Promise<void>((resolve) => {
            viewer.load(structure, properties, options, () => {
                // Resolve the Promise when the viewer has finished loading
                resolve();
            });
        });
    }

    /**
     * Get an unused color to identify a viewer
     *
     * @return a CSS compatible color name
     */
    private _getNextColor(): string {
        const colors = [];
        for (const data of this._cellsData.values()) {
            colors.push(data.color);
        }
        return getNextColor(colors);
    }

    /**
     * Function to setup the cell in the structure viewer grid.
     * Will generate a GUID string if one does not exist for the cell
     * and instantiate all necessary buttons.
     */
    private _setupCell(cellGUID: GUID, colNum: number, rowNum: number): string {
        const cellId = `gi-${cellGUID}`;
        let cell = this._root.querySelector(`#${cellId}`) as HTMLElement;
        let color = '';

        if (cell === null) {
            cell = document.createElement('div');
            cell.id = cellId;
            cell.classList.add('chsp-structure-viewer-cell', 'grid-item');
            cell.style.gridColumn = `${colNum}`;
            cell.style.gridRow = `${rowNum}`;

            // Do not set `cell.onclick = () => {this.active = GUID}`
            // because it will conflict with the button behavior

            // add a button to set the cell/viewer as the active one
            color = this._getNextColor();
            const template = document.createElement('template');
            template.innerHTML = `<div
                class="chsp-has-tooltip chsp-viewer-button chsp-viewer-action-button"
                style="position: absolute; top: 6px; left: 6px; display: flex">
                    <div id="chsp-activate-${cellGUID}"
                         class="chsp-pin-marker"
                         style="background-color: ${color}; position: relative; margin: auto;"
                    ></div>
                    <span class="chsp-tooltip chsp-tooltip-right">WILL BE FILLED LATER</span>
                </div>`;
            const activate = template.content.firstChild as HTMLElement;
            activate.onclick = () => {
                this.setActive(cellGUID);

                const data = this._cellsData.get(this._active);
                assert(data !== undefined);
                this.activeChanged(cellGUID, data.current);
            };
            cell.appendChild(activate);

            // add a button to remove the cell
            template.innerHTML = `<button
                class="btn btn-light btn-sm chsp-has-tooltip chsp-viewer-button chsp-viewer-action-button"
                style="top: 6px; right: 32px;">
                    <span>${CLOSE_SVG}</span>
                    <span class="chsp-tooltip chsp-tooltip-below">Remove viewer</span>
                </button>`;
            const remove = template.content.firstChild as HTMLElement;
            remove.onclick = () => {
                if (this._cellsData.size === 1) {
                    this.warnings.sendMessage('can not remove the last viewer from the grid');
                    return;
                }

                this.onremove(cellGUID);
                this.removeViewer(cellGUID);
                this._setupGrid(this._cellsData.size);
            };
            cell.appendChild(remove);

            // add a button to duplicate the cell
            template.innerHTML = `<button
                class="btn btn-light btn-sm chsp-has-tooltip chsp-viewer-button chsp-viewer-action-button"
                style="top: 6px; right: 60px;">
                    <span style="position: absolute; top: 2px; right: 5px;">${DUPLICATE_SVG}</span>
                    <span class="chsp-tooltip chsp-tooltip-below">Duplicate viewer</span>
                </button>`;
            const duplicate = template.content.firstChild as HTMLElement;

            duplicate.onclick = async () => {
                const data = await this._duplicate(cellGUID);
                if (data === undefined) {
                    return;
                }
                this.oncreate(this.active, data.color, data.current);
                this.activeChanged(this.active, data.current);
            };
            cell.appendChild(duplicate);

            // add a button to download PNG
            template.innerHTML = `<button
                class="btn btn-light btn-sm chsp-has-tooltip chsp-viewer-button chsp-viewer-action-button"
                style="top: 6px; right: 90px;">
                    <span>${PNG_SVG}</span>
                    <span class="chsp-tooltip chsp-tooltip-below">Download PNG</span>
                </button>`;
            const downloadPNG = template.content.firstChild as HTMLElement;

            downloadPNG.onclick = () => {
                const data = this._cellsData.get(cellGUID);
                assert(data !== undefined);
                const viewer = data.viewer;

                const structID = data.current.structure;
                if (data.current.atom !== undefined) {
                    const atomID = data.current.atom;
                    downloadURI(
                        viewer.exportPNG(),
                        `chemiscope-structure-${structID}-atom-${atomID}.png`
                    );
                } else {
                    downloadURI(viewer.exportPNG(), `chemiscope-structure-${structID}.png`);
                }
            };

            cell.appendChild(downloadPNG);

            this._root.appendChild(cell);
            return color;
        } else {
            cell.style.gridColumn = `${colNum}`;
            cell.style.gridRow = `${rowNum}`;
            return '';
        }
    }

    /**
     * Function to initialize the grid instance for `nViewers` cells and place
     * onto the DOM element mapped in `this._root`. If more cells are needed,
     * this function return the list of new cell GUID
     */
    private _setupGrid(nViewers: number): GUID[] {
        const newGUID = [] as GUID[];
        if (nViewers < 1) {
            this.warnings.sendMessage('Cannot delete last molecular viewer.');
            return newGUID;
        } else if (nViewers > this._maxViewers) {
            this.warnings.sendMessage(
                `Viewer grid cannot contain more than ${this._maxViewers} viewers.`
            );
            return newGUID;
        }

        // Determine best arrangement for the number of viewers requested
        const arrangement = this.bestGridArrangement(nViewers);
        if (this._cellsData.size > nViewers) {
            this.warnings.sendMessage(
                `Eliminating last ${this._cellsData.size - nViewers} viewers.`
            );
            let i = 0;
            for (const guid of this._cellsData.keys()) {
                if (i >= nViewers) {
                    this.removeViewer(guid);
                }
                i += 1;
            }
        }

        // Start at the third row to skip the inc/dec buttons
        let rowNum = 1;
        let colNum = 1;

        const mapKeys = this._cellsData.keys();
        for (let c = 0; c < nViewers; c++) {
            let cellGUID: GUID;
            if (c >= this._cellsData.size) {
                cellGUID = generateGUID();
            } else {
                cellGUID = mapKeys.next().value as GUID;
            }

            let color = this._setupCell(cellGUID, colNum, rowNum);
            if (color === '') {
                color = this._getNextColor();
            }

            colNum++;
            if (colNum > arrangement.columns) {
                rowNum++;
                colNum = 1;
            }

            // add a new cells if necessary
            if (!this._cellsData.has(cellGUID)) {
                const propertiesName = this._properties ? Object.keys(this._properties) : [];
                const viewer = new MoleculeViewer(
                    this._getById(`gi-${cellGUID}`),
                    propertiesName,
                    this.warnings
                );

                viewer.onSettingChange((keys, value) => {
                    if (keys[0] === 'playbackDelay' && this._active === cellGUID) {
                        this.delayChanged(value as number);
                    }
                });

                viewer.onselect = (atom: number) => {
                    if (this._target !== 'atom' || this._active !== cellGUID) {
                        return;
                    }

                    // if the viewer is showing a bigger supercell than [1, 1, 1], the
                    // atom index can be outside of [0, natoms), so make sure it is
                    // inside this range.
                    const data = this._cellsData.get(this._active);
                    assert(data !== undefined);
                    const natoms = viewer.natoms();
                    assert(natoms !== undefined);
                    const indexes = this._indexer.fromStructureAtom(
                        this._target,
                        data.current.structure,
                        atom % natoms
                    );

                    if (indexes === undefined) {
                        // user clicked on an atom which is not in this dataset
                        return;
                    }

                    viewer.highlight(atom);
                    this.onselect(indexes);
                };

                const current = { atom: undefined, structure: -1, environment: -1 };

                this._cellsData.set(cellGUID, {
                    color: color,
                    current: current,
                    viewer: viewer,
                });

                if (this._positionSettingsModal !== undefined) {
                    viewer.positionSettingsModal = this._positionSettingsModal;
                }
                newGUID.push(cellGUID);
            }
        }

        // Force a refresh of the viewer in case the aspect ratio changed
        for (const data of this._cellsData.values()) {
            data.viewer.resize();
        }

        return newGUID;
    }

    /**
     * Function to return *optimal* arrangement of n viewers in a grid.
     */
    private bestGridArrangement(n: number) {
        switch (n) {
            case 1:
            case 2:
                return { rows: n, columns: 1 };
            case 3:
            case 4:
                return { rows: 2, columns: 2 };
            case 5:
            case 6:
                return { rows: 3, columns: 2 };
            case 8:
            case 7:
            case 9:
                return { rows: 3, columns: 3 };
            default:
                throw Error('reached unreachable code: too many viewers in the grid');
        }
    }

    /** get the current position in the viewer grid of a viewer by GUID */
    private _getCurrentPositionInGrid(guid: GUID): number {
        let position = 0;
        for (const currentGuid of this._cellsData.keys()) {
            if (guid === currentGuid) {
                break;
            }
            position += 1;
        }
        return position;
    }
}
