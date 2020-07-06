/**
 * @packageDocumentation
 * @module structure
 */

import assert from 'assert';

import {checkStructure, Environment, JsObject, Structure, UserStructure} from '../dataset';

import {EnvironmentIndexer, Indexes} from '../indexer';
import {SettingsPreset} from '../settings';
import {GUID, PositioningCallback} from '../utils';
import {generateGUID, getByID, getFirstKey, getNextColor, sendWarning} from '../utils';

import {structure2JSmol} from './utils';
import {JSmolWidget, LoadOptions} from './widget';

import CLOSE_SVG from '../static/close.svg';
import DUPLICATE_SVG from '../static/duplicate.svg';

const MAX_WIDGETS = 9;

/**
 * Function to return *optimal* arrangement of n widgets.
 * Defaults to a 4 x ? grid. Hardcoded for simplicity.
 */
function bestGridArrangement(n: number) {
    switch (n) {
    case 1:
    case 2:
        return {rows: n, columns: 1};
    case 3:
    case 4:
        return {rows: 2, columns: 2};
    case 5:
    case 6:
        return {rows: 3, columns: 2};
    case 8:
    case 7:
    case 9:
        return {rows: 3, columns: 3};
    default:
        throw Error(`Can not create a grid with more than ${MAX_WIDGETS} elements`);
    }
}

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
function groupByStructure(n_structures: number, environments?: Environment[]): Environment[][] | undefined {
    if (environments === undefined) {
        return undefined;
    }

    const result: Environment[][] = [];
    for (let i = 0; i < n_structures; i++) {
        result.push([]);
    }

    for (const env of environments) {
        result[env.structure].push(env);
    }

    return result;
}

interface WidgetGridData {
    widget: JSmolWidget;
    color: string;
    current: Indexes;
}

/**
 * The [[ViewersGrid]] class displays a grid of molecule or a crystal viewers
 * in 3D using [JSmol](http://wiki.jmol.org/index.php/JSmol) for rendering.
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
     * Callback used when a new structure should be loaded
     *
     * By default, this assumes that the loaded dataset contains [[Structure]]
     * directly, and returns the data from there. If the loaded dataset contains
     * [[UserStructure]] instead, this callback should be set to transform from
     * data in [[UserStructure.data]] to a [[Structure]].
     *
     * The callback gets two parameter: the structure index (0-based); and the
     * full [[UserStructure]].
     */
    public loadStructure: (index: number, structure: unknown) => Structure;

    /// Root element containing the full viewer grid
    private _root: HTMLElement;
    /// Playback delay setting
    private _delay: HTMLInputElement;
    /// List of structures in the dataset
    private _structures: Structure[] | UserStructure[];
    /// Cached string representation of structures
    private _cachedStructures: string[];
    /// Optional list of environments for each structure
    private _environments?: Environment[][];
    /// The indexer translating between environments indexes and structure/atom
    /// indexes
    private _indexer: EnvironmentIndexer;
    // path to the j2s files used by JSmol.
    // saved for instantiating new Widget instances
    private _j2spath: string;
    // storage for the overall presets
    private _presets: SettingsPreset;
    /// GUID of the currently active widget
    private _active: GUID;
    /// Map of Widgets GUIDS to their color, widget, and current structure
    private _viewers: Map<GUID, WidgetGridData>;
    /// Callback used to override all grid viewers' positionSettingsModal
    private _positionSettingsModal?: PositioningCallback;

    /**
     * Create a new [[ViewersGrid]] inside the HTML element with the given
     * `id`
     *
     * @param id           HTML id of the DOM element where the viewer should live
     * @param j2sPath      path to the `j2s` files uses by JSmol
     * @param indexer      [[EnvironmentIndexer]] used to translate indexes from
     *                     environments index to structure/atom indexes
     * @param structures   list of structure to display
     * @param environments list of atom-centered environments in the structures,
     *                     used to highlight the selected environment
     */
    constructor(
        config: { id: string, presets: SettingsPreset },
        j2sPath: string,
        indexer: EnvironmentIndexer,
        structures: Structure[] | UserStructure[],
        environments?: Environment[],
    ) {
        this._structures = structures;
        this._cachedStructures = new Array(structures.length);
        this._environments = groupByStructure(this._structures.length, environments);
        this._indexer = indexer;
        this._j2spath = j2sPath;
        this._presets = config.presets;

        this.loadStructure = (_, s) => {
            // check that the data does conform to the Structure interface
            if (checkStructure(s as JsObject) !== '') {
                throw Error(
                    'got custom data for this structure, but no custom loadStructure callback\n' +
                    `the object was ${JSON.stringify(s)}`,
                );
            } else {
                return s as Structure;
            }
        };

        // Initializes with 1 widget upon opening.
        this._viewers = new Map();
        this.onselect = () => {};
        this.onremove = () => {};
        this.oncreate = () => {};
        this.activeChanged = () => {};

        const root = getByID(config.id);
        this._root = document.createElement('div');
        this._root.id = 'grid-root';
        this._root.className = 'chsp-structure-viewer-grid';
        root.appendChild(this._root);

        this._setupGrid(1);

        this._active = getFirstKey(this._viewers);
        this.setActive(this._active);

        // get the 'delay' setting inside the current widget setting
        this._delay = getByID<HTMLInputElement>(`chsp-${this._active}-playback-delay`);
    }

    /**
     * Get the current active viewer GUID
     */
    public get active(): GUID {
        return this._active;
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
    public show(indexes: Indexes) {
        const data = this._viewers.get(this._active);
        assert(data !== undefined);

        const widget = data.widget;
        if (data.current.structure !== indexes.structure) {
            const options: Partial<LoadOptions> = {
                packed: false,
                trajectory: true,
            };
            assert(indexes.structure < this._structures.length);

            if (this._environments !== undefined) {
                options.environments = this._environments[indexes.structure];
                if (this._indexer.mode === 'atom') {
                    options.highlight = indexes.atom;
                }
            }

            widget.load(`inline '${this._structureForJSmol(indexes.structure)}'`, options);
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
     * Start playing the trajectory of structures in this dataset, until
     * `advance` returns false
     */
    public structurePlayback(advance: () => boolean) {
        setTimeout(() => {
            if (advance()) {
                const widgetData = this._viewers.get(this._active);
                assert(widgetData !== undefined);

                const current = widgetData.current;
                const structure = (current.structure + 1) % this._indexer.structuresCount();
                const indexes = this._indexer.from_structure_atom(structure, 0);
                this.show(indexes);
                this.onselect(indexes);
                // continue playing until the advance callback returns false
                this.structurePlayback(advance);

            }
        }, parseFloat(this._delay.value) * 100);
    }

    /**
     * Start playing the 'trajectory' of atoms in the current structure, until
     * `advance` returns false
     */
    public atomPlayback(advance: () => boolean) {
        setTimeout(() => {
            if (advance()) {
                const widgetData = this._viewers.get(this._active);
                assert(widgetData !== undefined);

                const current = widgetData.current;
                const structure = current.structure;
                const atom = (current.atom! + 1) % this._indexer.atomsCount(structure);
                const indexes = this._indexer.from_structure_atom(structure, atom);
                this.show(indexes);
                this.onselect(indexes);
                // continue playing until the advance callback returns false
                this.atomPlayback(advance);
            }
        }, parseFloat(this._delay.value) * 100);
    }

    /**
     * Remove all HTML added by this [[ViewersGrid]] in the current document
     */
    public remove(): void {
        for (const data of this._viewers.values()) {
            data.widget.remove();
        }

        if (this._root.parentElement !== null) {
            this._root.parentElement.innerHTML = '';
        }
    }

    /**
     * Function to set the active widget for communicating with the map
     */
    public setActive(guid: GUID): void {
        const changeClasses = (toggle: boolean) => {
            assert(this._viewers.has(this._active));
            // change tooltip text in the active marker
            const button = getByID(`chsp-activate-${this._active}`);
            button.classList.toggle('chsp-active-structure', toggle);
            const tooltip = button.parentElement!.querySelector('.chsp-tooltip');
            assert(tooltip !== null);
            tooltip.innerHTML = toggle ? 'this is the active viewer' : 'choose as active';

            // change style of the cell border
            const cell = getByID(`gi-${this._active}`);
            cell.classList.toggle('chsp-structure-viewer-cell-active', toggle);
        };

        // remove active classes from the previous active
        changeClasses(false);

        this._active = guid;
        changeClasses(true);

        this._delay = getByID<HTMLInputElement>(`chsp-${this._active}-playback-delay`);
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

    /*
    * Removes a widget from the structure viewer grid.
    * Will not remove a widget if it is the last one in the structure viewer
    */
    private _removeWidget(guid: GUID) {
        if (this._viewers.size === 1) {
            return;
        }

        // If we removed the active marker, change the active one
        if (this._active === guid) {
            this.setActive(getFirstKey(this._viewers, guid));
        }

        // remove HTML inserted by the widget
        const data = this._viewers.get(guid);
        assert(data !== undefined);
        data.widget.remove();

        // remove the cell containing the widget
        const cell = getByID(`gi-${guid}`);
        cell.remove();

        this._viewers.delete(guid);
    }

    /**
     * Get the structure at the given index in a format JSmol can undertand
     * and load. [[Structure]] already rendered as strings are cached for faster
     * subsequent access.
     * @param  index index of the structure
     * @return       a string that can be passed to JSmol `load INLINE` command
     */
    private _structureForJSmol(index: number): string {
        if (this._cachedStructures[index] === undefined) {
            const s = this.loadStructure(index, this._structures[index]);
            const check = checkStructure(s as unknown as JsObject);
            if (check !== '') {
                throw Error(
                    `got invalid object as structure: ${check}` + '\n' +
                    `the object was ${JSON.stringify(s)}`,
                );
            }
            this._cachedStructures[index] = structure2JSmol(s);
        }
        return this._cachedStructures[index];
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
        let cell = document.getElementById(cellId);
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
                style="position: absolute; top: 9px; right: 115px;">
                    <div id="chsp-activate-${cellGUID}"
                         class="chsp-structure-marker"
                         style="background-color: ${color}; top: 14px; right: 0px;"
                    ></div>
                    <span class="chsp-tooltip">WILL BE FILLED LATER</span>
                </div>`;
            const activate = template.content.firstChild! as HTMLElement;
            activate.onclick = () => {
                this.setActive(cellGUID);

                const data = this._viewers.get(this._active);
                assert(data !== undefined);
                this.activeChanged(cellGUID, data.current);
            };
            cell.appendChild(activate);

            // add a button to remove the widget
            template.innerHTML = `<button
                class="btn btn-light btn-sm chsp-has-tooltip chsp-viewer-button"
                style="top: 8px; right: 40px;">
                    <span>${CLOSE_SVG}</span>
                    <span class="chsp-tooltip">remove viewer</span>
                </button>`;
            const remove = template.content.firstChild! as HTMLElement;
            remove.onclick = () => {
                this.onremove(cellGUID);
                this._removeWidget(cellGUID);
                this._setupGrid(this._viewers.size);
            };
            cell.appendChild(remove);

            // add a button to duplicate the widget
            template.innerHTML = `<button
                class="btn btn-light btn-sm chsp-has-tooltip chsp-viewer-button"
                style="top: 8px; right: 70px;">
                    <span>${DUPLICATE_SVG}</span>
                    <span class="chsp-tooltip">duplicate viewer</span>
                </button>`;
            const duplicate = template.content.firstChild! as HTMLElement;

            duplicate.onclick = () => {
                const newGUID = this._setupGrid(this._viewers.size + 1);
                if (newGUID.length === 0) {
                    // no new widget, probably because we already have MAX_WIDGETS
                    return;
                }
                assert(newGUID.length === 1);
                this.setActive(newGUID[0]);

                const data = this._viewers.get(cellGUID);
                assert(data !== undefined);
                const indexes = this._indexer.from_structure_atom(data.current.structure, data.current.atom);
                this.show(indexes);

                const newData = this._viewers.get(newGUID[0]);
                assert(newData !== undefined);
                this.oncreate(newGUID[0], newData.color, indexes);
            };

            cell.appendChild(duplicate);

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
    private _setupGrid(nwidgets: number): GUID[] {
        const newGUID = [] as GUID[];
        if (nwidgets < 1) {
            sendWarning('Cannot delete last widget.');
            return newGUID;
        } else if (nwidgets > MAX_WIDGETS ) {
            sendWarning(`Viewer grid cannot contain more than ${MAX_WIDGETS} widgets.`);
            return newGUID;
        }

        // Determine best arrangement for nwidgets
        const arrangement = bestGridArrangement(nwidgets);
        if (this._viewers.size > nwidgets) {
            sendWarning(`Warning: Eliminating last ${this._viewers.size - nwidgets} viewers.`);
            let i = 0;
            for (const guid of this._viewers.keys()) {
                if (i >= nwidgets) {
                    this._removeWidget(guid);
                }
                i += 1;
            }
        }

        // Start at the third row to skip the inc/dec buttons
        let rowNum = 1;
        let colNum = 1;

        const mapKeys = this._viewers.keys();
        for (let c = 0; c < nwidgets; c++) {
            let cellGUID: GUID;
            if (c >= this._viewers.size) {
                cellGUID = generateGUID();
            } else {
                cellGUID = mapKeys.next().value;
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
                const widget = new JSmolWidget(
                    `gi-${cellGUID}`,
                    this._j2spath,
                    cellGUID,
                );
                widget.applyPresets(this._presets);

                widget.onselect = (atom: number) => {
                    if (this._indexer.mode !== 'atom' || this._active !== cellGUID) {
                        return;
                    }

                    widget.highlight(atom);

                    // if the viewer is showing a bigger supercell than [1, 1, 1], the
                    // atom index can be outside of [0, natoms), so make sure it is
                    // inside this range.
                    const data = this._viewers.get(this._active);
                    assert(data !== undefined);
                    const atom_id = atom % widget.natoms()!;
                    const indexes = this._indexer.from_structure_atom(data.current.structure, atom_id);
                    this.onselect(indexes);
                };

                const current = {atom: undefined, structure: -1, environment: -1};
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
            widgetData.widget.script('refresh');
        }

        return newGUID;
    }
}