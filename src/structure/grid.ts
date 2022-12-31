/**
 * @packageDocumentation
 * @module structure
 */

import assert from 'assert';

import {
    Environment,
    JsObject,
    Settings,
    Structure,
    UserStructure,
    checkStructure,
} from '../dataset';

import { EnvironmentIndexer, Indexes } from '../indexer';
import * as styles from '../styles';
import { GUID, PositioningCallback, getElement } from '../utils';
import { enumerate, generateGUID, getByID, getFirstKey, getNextColor, sendWarning } from '../utils';

import { LoadOptions, MoleculeViewer } from './viewer';

import CLOSE_SVG from '../static/close.svg';
import DUPLICATE_SVG from '../static/duplicate.svg';
import PNG_SVG from '../static/download-png.svg';

/**
 * Create a list of environments grouped together by structure.
 *
 *
 * This function returns `undefined` if `environment` is undefined, else it
 * returns a list of list of environments, such as `list[0]` contains all
 * environments in structure 0; `list[33]` all environments in structure 33, etc.
 *
 * @param  n_structures Expected number of structures
 * @param  environments Full list of environments
 *
 * @return              The list of environments grouped by structure
 */
function groupByStructure(
    structures: (Structure | UserStructure)[],
    environments?: Environment[]
): Environment[][] | undefined {
    if (environments === undefined) {
        return undefined;
    }

    const result = Array.from({ length: structures.length }).map((_, i) =>
        Array.from({ length: structures[i].size })
    );

    for (const env of environments) {
        result[env.structure][env.center] = env;
    }

    return result as Environment[][];
}

interface WidgetGridData {
    /// the viewer itself
    widget: MoleculeViewer;
    /// color associated with this viewer
    color: string;
    /// set of indexes currently displayed in this viewer
    current: Indexes;
}

/**
 * The {@link ViewersGrid} class displays a grid of molecule or a crystal viewers
 * in 3D using {@link MoleculeViewer} widgets for rendering.
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
    public loadStructure: (index: number, structure: unknown) => Structure;

    /// Shadow root for isolation
    private _shadow: ShadowRoot;
    /// Root element containing the full viewer grid
    private _root: HTMLElement;
    /// List of structures in the dataset
    private _structures: Structure[] | UserStructure[];
    /// Cached string representation of structures
    private _resolvedStructures: Structure[];
    /// Optional list of environments for each structure
    private _environments?: Environment[][];
    /// Maximum number of allowed structure viewers
    private _maxViewers: number;
    /// The indexer translating between environments indexes and structure/atom
    /// indexes
    private _indexer: EnvironmentIndexer;
    /// GUID of the currently active widget
    private _active: GUID;
    /// Map of Widgets GUIDS to their color, widget, and current structure
    private _viewers: Map<GUID, WidgetGridData>;
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
     * @param environments list of atom-centered environments in the structures,
     *                     used to highlight the selected environment
     * @param maxViewers   maximum number of allowed structure viewers
     */
    constructor(
        element: string | HTMLElement,
        indexer: EnvironmentIndexer,
        structures: Structure[] | UserStructure[],
        environments?: Environment[],
        maxViewers: number = 9
    ) {
        this._structures = structures;
        this._resolvedStructures = new Array<Structure>(structures.length);
        this._environments = groupByStructure(this._structures, environments);
        this._indexer = indexer;

        this.loadStructure = (_, s) => {
            // check that the data does conform to the Structure interface
            if (checkStructure(s as JsObject) !== '') {
                throw Error(
                    'got custom data for this structure, but no custom loadStructure callback\n' +
                        `the object was ${JSON.stringify(s)}`
                );
            } else {
                return s as Structure;
            }
        };

        // Initializes with 1 widget upon opening.
        this._viewers = new Map<GUID, WidgetGridData>();
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

        const containerElement = getElement(element);
        const hostElement = document.createElement('div');

        hostElement.style.setProperty('height', '100%');
        containerElement.appendChild(hostElement);

        this._shadow = hostElement.attachShadow({ mode: 'open' });
        this._shadow.adoptedStyleSheets = [styles.bootstrap, styles.chemiscope];

        this._root = document.createElement('div');
        this._root.id = 'grid-root';
        this._root.className = 'chsp-structure-viewer-grid';
        this._shadow.appendChild(this._root);

        this._setupGrid(1);

        this._active = getFirstKey(this._viewers);
        this.setActive(this._active);
    }

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
     * Get the current list of environments showed inside the different viewer
     */
    public pinned(): Indexes[] {
        const result = [];
        for (const data of this._viewers.values()) {
            result.push(data.current);
        }
        return result;
    }

    /**
     * Show a new structure, as identified by `indexes` in the active viewer.
     *
     * This will switch to the structure at index `indexes.structure`, and if
     * environments where passed to the constructor and the current display
     * mode is `'atom'`, highlight the atom-centered environment corresponding
     * to `indexes.atom`.
     *
     * @param  indexes         structure / atom pair to display
     */
    public show(indexes: Indexes): void {
        this._showInViewer(this.active, indexes);
    }

    /**
     * Add a new empty viewer to the grid
     *
     * @return the GUID and color of the new viewer, the GUID is `undefined`
     *         if we already reached the viewer limit.
     */
    public addViewer(): { guid?: GUID; color: string } {
        const newGUIDs = this._setupGrid(this._viewers.size + 1);
        if (newGUIDs.length === 0) {
            // no new widget, probably because we already have MAX_WIDGETS
            return { guid: undefined, color: '' };
        }
        assert(newGUIDs.length === 1);
        const newGUID = newGUIDs[0];

        const newData = this._viewers.get(newGUID);
        assert(newData !== undefined);

        for (const callback of this._onSettingChangeCallbacks) {
            newData.widget.onSettingChange((keys, value) => {
                // add the index of the viewer in the current grid of viewers
                // at the front of the keys
                let i = 0;
                for (const currentGuid of this._viewers.keys()) {
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
        assert(this._viewers.size > 1);

        // If we removed the active marker, change the active one
        if (this._active === guid) {
            this.setActive(getFirstKey(this._viewers, guid));
        }

        // remove HTML inserted by the widget
        const data = this._viewers.get(guid);
        assert(data !== undefined);
        data.widget.remove();

        // remove the cell containing the widget
        const cell = this._getById(`gi-${guid}`);
        cell.remove();

        this._viewers.delete(guid);
    }

    /**
     * Remove all HTML added by this {@link ViewersGrid} in the current document
     */
    public remove(): void {
        for (const data of this._viewers.values()) {
            data.widget.remove();
        }

        this._shadow.host.remove();
    }

    /**
     * Function to set the active widget for communicating with the map
     */
    public setActive(guid: GUID): void {
        const changeClasses = (guid: string, toggle: boolean) => {
            // change tooltip text in the active marker
            const button = this._getById(`chsp-activate-${guid}`);
            button.classList.toggle('chsp-active-structure', toggle);
            assert(button.parentElement !== null);
            const tooltip = button.parentElement.querySelector('.chsp-tooltip');
            assert(tooltip !== null);
            tooltip.innerHTML = toggle ? 'Active viewer' : 'Choose as active';

            // change style of the cell border
            const cell = this._getById(`gi-${guid}`);
            cell.classList.toggle('chsp-structure-viewer-cell-active', toggle);
        };

        const current = this._viewers.get(this._active);
        assert(current !== undefined);

        // remove active classes from the previous active viewer
        changeClasses(this._active, false);

        this._active = guid;

        const newViewer = this._viewers.get(this._active);
        assert(newViewer !== undefined);

        // links playback delay options
        newViewer.widget._options.playbackDelay.onchange.push((value) => {
            this.delayChanged(value);
        });

        // set the right initial value for playback delay
        this.delayChanged(newViewer.widget._options.playbackDelay.value);

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
        for (const viewer of this._viewers.values()) {
            viewer.widget.positionSettingsModal = value;
        }
    }

    /**
     * Apply the given saved setting to all viewers in the grid.
     *
     * The settings must be in viewer creation order.
     *
     * @param settings settings for all viewers in the grid
     */
    public applySettings(settings: Settings[]): void {
        if (settings.length === 0) {
            return;
        }

        assert(settings.length === this._viewers.size);
        for (const [i, data] of enumerate(this._viewers.values())) {
            data.widget.applySettings(settings[i]);
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
        for (const data of this._viewers.values()) {
            settings.push(data.widget.saveSettings());
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
        for (const guid of this._viewers.keys()) {
            const data = this._viewers.get(guid);
            assert(data !== undefined);

            data.widget.onSettingChange((keys: string[], value: unknown) => {
                keys = JSON.parse(JSON.stringify(keys)) as string[];
                keys.unshift(this._getCurrentPositionInGrid(guid).toString());

                callback(keys, value);
            });
        }

        this._onSettingChangeCallbacks.push(callback);
    }

    /**
     * Add a new structure viewer to the grid as a copy of the viewer with the
     * `initial` GUID. The new structure viewer is set as the active one.
     *
     * @param  initial GUID of the viewer to duplicate
     * @return the GUID of the new viewer, or `undefined` if we already
     *         reached the viewer limit.
     */
    private _duplicate(initial: GUID): WidgetGridData | undefined {
        const data = this._viewers.get(initial);
        assert(data !== undefined);

        const newGUID = this.addViewer().guid;
        if (newGUID === undefined) {
            return undefined;
        }
        const newData = this._viewers.get(newGUID);
        assert(newData !== undefined);

        newData.widget.applySettings(data.widget.saveSettings());
        this._showInViewer(newGUID, data.current);
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
    private _structure(index: number): Structure {
        if (this._resolvedStructures[index] === undefined) {
            const s = this.loadStructure(index, this._structures[index]);
            const check = checkStructure(s as unknown as JsObject);
            if (check !== '') {
                throw Error(
                    `got invalid object as structure: ${check}
    the object was ${JSON.stringify(s)}`
                );
            }
            this._resolvedStructures[index] = s;
        }
        return this._resolvedStructures[index];
    }

    private _showInViewer(guid: GUID, indexes: Indexes): void {
        const data = this._viewers.get(guid);
        assert(data !== undefined);

        const widget = data.widget;
        if (data.current.structure !== indexes.structure) {
            const options: Partial<LoadOptions> = {
                trajectory: true,
            };
            assert(indexes.structure < this._structures.length);

            if (this._environments !== undefined) {
                options.environments = this._environments[indexes.structure];
                if (this._indexer.mode === 'atom') {
                    options.highlight = indexes.atom;
                }
            }

            widget.load(this._structure(indexes.structure), options);
            data.current = indexes;
        }

        if (this._indexer.mode === 'atom') {
            if (data.current.atom !== indexes.atom) {
                widget.highlight(indexes.atom);
            }
        } else {
            widget.highlight(undefined);
        }

        data.current = indexes;
    }

    /**
     * Get an unused color to identify a viewer
     *
     * @return a CSS compatible color name
     */
    private _getNextColor(): string {
        const colors = [];
        for (const data of this._viewers.values()) {
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

            // add a button to activate the widget (i.e. set as `active`)
            const template = document.createElement('template');
            color = this._getNextColor();
            template.innerHTML = `<div
                class="chsp-has-tooltip"
                style="position: absolute; top: 7px; right: 160px;">
                    <div id="chsp-activate-${cellGUID}"
                         class="chsp-structure-marker"
                         style="background-color: ${color}; top: 13px; right: 0px;"
                    ></div>
                    <span class="chsp-tooltip">WILL BE FILLED LATER</span>
                </div>`;
            const activate = template.content.firstChild as HTMLElement;
            activate.onclick = () => {
                this.setActive(cellGUID);

                const data = this._viewers.get(this._active);
                assert(data !== undefined);
                this.activeChanged(cellGUID, data.current);
            };
            cell.appendChild(activate);

            // add a button to remove the widget
            template.innerHTML = `<button
                class="btn btn-light btn-sm chsp-has-tooltip chsp-viewer-button chsp-viewer-action-button"
                style="top: 6px; right: 41px;">
                    <span>${CLOSE_SVG}</span>
                    <span class="chsp-tooltip">Remove viewer</span>
                </button>`;
            const remove = template.content.firstChild as HTMLElement;
            remove.onclick = () => {
                if (this._viewers.size === 1) {
                    sendWarning('can not remove the last viewer from the grid');
                    return;
                }

                this.onremove(cellGUID);
                this.removeViewer(cellGUID);
                this._setupGrid(this._viewers.size);
            };
            cell.appendChild(remove);

            // add a button to duplicate the widget
            template.innerHTML = `<button
                class="btn btn-light btn-sm chsp-has-tooltip chsp-viewer-button chsp-viewer-action-button"
                style="top: 6px; right: 76px;">
                    <span>${DUPLICATE_SVG}</span>
                    <span class="chsp-tooltip">Duplicate viewer</span>
                </button>`;
            const duplicate = template.content.firstChild as HTMLElement;

            duplicate.onclick = () => {
                const data = this._duplicate(cellGUID);
                if (data === undefined) {
                    return;
                }
                this.oncreate(this.active, data.color, data.current);
            };
            cell.appendChild(duplicate);

            // add a button to download PNG
            template.innerHTML = `<button
                class="btn btn-light btn-sm chsp-has-tooltip chsp-viewer-button chsp-viewer-action-button"
                style="top: 6px; right: 111px;">
                    <span>${PNG_SVG}</span>
                    <span class="chsp-tooltip">Download PNG</span>
                </button>`;
            const downloadPNG = template.content.firstChild as HTMLElement;

            downloadPNG.onclick = () => {
                const viewer = this._viewers.get(cellGUID);
                assert(viewer !== undefined);
                const widget = viewer.widget;

                const structID = viewer.current.structure;
                if (viewer.current.atom !== undefined) {
                    const atomID = viewer.current.atom;
                    downloadURI(
                        widget.exportPNG(),
                        `chemiscope-structure-${structID + 1}-atom-${atomID + 1}.png`
                    );
                } else {
                    downloadURI(widget.exportPNG(), `chemiscope-structure-${structID + 1}.png`);
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
     * Function to initialize the grid instance for `this._nwidgets` cells and
     * place onto the DOM element mapped in `this._root`. If more cells are
     * needed, this function return the list of new cell GUID
     */
    private _setupGrid(nViewers: number): GUID[] {
        const newGUID = [] as GUID[];
        if (nViewers < 1) {
            sendWarning('Cannot delete last widget.');
            return newGUID;
        } else if (nViewers > this._maxViewers) {
            sendWarning(`Viewer grid cannot contain more than ${this._maxViewers} widgets.`);
            return newGUID;
        }

        // Determine best arrangement for the number of viewers requested
        const arrangement = this.bestGridArrangement(nViewers);
        if (this._viewers.size > nViewers) {
            sendWarning(`Warning: Eliminating last ${this._viewers.size - nViewers} viewers.`);
            let i = 0;
            for (const guid of this._viewers.keys()) {
                if (i >= nViewers) {
                    this.removeViewer(guid);
                }
                i += 1;
            }
        }

        // Start at the third row to skip the inc/dec buttons
        let rowNum = 1;
        let colNum = 1;

        const mapKeys = this._viewers.keys();
        for (let c = 0; c < nViewers; c++) {
            let cellGUID: GUID;
            if (c >= this._viewers.size) {
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

            // add a new widget if necessary
            if (!this._viewers.has(cellGUID)) {
                const widget = new MoleculeViewer(this._getById<HTMLElement>(`gi-${cellGUID}`));

                widget.onselect = (atom: number) => {
                    if (this._indexer.mode !== 'atom' || this._active !== cellGUID) {
                        return;
                    }

                    // if the viewer is showing a bigger supercell than [1, 1, 1], the
                    // atom index can be outside of [0, natoms), so make sure it is
                    // inside this range.
                    const data = this._viewers.get(this._active);
                    assert(data !== undefined);
                    const natoms = widget.natoms();
                    assert(natoms !== undefined);
                    const indexes = this._indexer.from_structure_atom(
                        data.current.structure,
                        atom % natoms
                    );

                    if (indexes === undefined) {
                        // user clicked on an atom which is not in this dataset
                        return;
                    }

                    widget.highlight(atom);
                    this.onselect(indexes);
                };

                const current = { atom: undefined, structure: -1, environment: -1 };

                this._viewers.set(cellGUID, {
                    color: color,
                    current: current,
                    widget: widget,
                });

                if (this._positionSettingsModal !== undefined) {
                    widget.positionSettingsModal = this._positionSettingsModal;
                }
                newGUID.push(cellGUID);
            }
        }

        // Force a refresh of the viewer in case the aspect ratio changed
        for (const widgetData of this._viewers.values()) {
            widgetData.widget.resize();
        }

        return newGUID;
    }

    /**
     * Function to return *optimal* arrangement of n widgets.
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
                throw Error('reached unreachable code: too many viewer in the grid');
        }
    }

    /** get the current position in the viewer grid of a viewer by GUID */
    private _getCurrentPositionInGrid(guid: GUID): number {
        let position = 0;
        for (const currentGuid of this._viewers.keys()) {
            if (guid === currentGuid) {
                break;
            }
            position += 1;
        }
        return position;
    }
}

/**
 *  Creates a download request from a URI
 * @param uri       URI of the image
 * @param name      Name of the downloaded image
 */
function downloadURI(uri: string, name: string) {
    const link = document.createElement('a');
    link.download = name;
    link.href = uri;
    document.body.appendChild(link);
    link.click();
    document.body.removeChild(link);
    link.remove();
}
