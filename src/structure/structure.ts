/**
 * @packageDocumentation
 * @module structure
 */

import assert from 'assert';

import {Environment, isPositiveInteger, Structure, UserStructure} from '../dataset';
import {EnvironmentIndexer, Indexes, PositioningCallback} from '../utils';
import {generateGUID, getByID, getNextColor, sendWarning} from '../utils';

import {structure2JSmol} from './jsmol';
import {JSmolWidget} from './widget';

import CLOSE_SVG from '../static/close.svg';
import DUPLICATE_SVG from '../static/duplicate.svg';

const MAX_WIDGETS = 9;

/*
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

/**
 * Check that the given object is a structure. Return a string describing the
 * issue with `s` if any, or the empty string if `s` looks like a valid
 * structure.
 */
function checkStructure(s: any): string {
    if (!('size' in s && typeof s.size === 'number' && isPositiveInteger(s.size))) {
        return 'missing "size" in structure';
    }

    for (const key of ['names', 'x', 'y', 'z']) {
        if (!(key in s && s[key].length !== undefined)) {
            return `missing '${name}' in structure`;
        }
        if (s[key].length !== s.size) {
            return `wrong size for '${name}' in structure, expected ${s.size}, got ${s[name].length}`;
        }
    }

    if ('cell' in s) {
        if (s.cell.length !== 9) {
            return '"cell" must be an array of size 9 in structure';
        }
    }

    return '';
}

interface WidgetGridData {
    widget: JSmolWidget;
    color: string;
    current: Indexes;
}

/**
 * The [[StructureViewer]] class displays a molecule or a crystal in 3D using
 * [JSmol](http://wiki.jmol.org/index.php/JSmol) for rendering.
 */
export class StructureViewer {
    /** Callback used when the user select an environment */
    public onselect: (indexes: Indexes, selectedGUID?: string) => void;

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
    public loadStructure: (index: number, structure: any) => Structure;

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
    /// GUID of the Active Widget
    private _active: string;
    /// Map of Widgets GUIDS to their color, widget, and current structure
    private _viewers: Map<string, WidgetGridData>;
    /// Callback used to override all grid viewers' positionSettingsModal
    private _positionSettingsModal?: PositioningCallback;

    /**
     * Create a new [[StructureViewer]] inside the HTML element with the given
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
        id: string,
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

        this.loadStructure = (_, s) => {
            // check that the data does conform to the Structure interface
            if (checkStructure(s) !== '') {
                throw Error(
                    'got custom data for this structure, but no custom loadStructure callback\n' +
                    `the object was ${JSON.stringify(s)}`,
                );
            } else {
                return s;
            }
        };

        // Initializes with 1 widget upon opening.
        this._viewers = new Map();
        this.onselect = () => {};

        const root = getByID(id);
        this._root = document.createElement('div');
        this._root.id = 'grid-root';
        this._root.className = 'chsp-structure-viewer-grid';
        root.appendChild(this._root);

        this._setupGrid(1);
        this.active = this._viewers.keys().next().value;
        this._active = this._viewers.keys().next().value;

        // get the 'delay' setting inside the current widget setting
        this._delay = getByID<HTMLInputElement>(`chsp-${this._active}-playback-delay`);
    }

    /**
     * Show a new structure, as identified by `indexes`. This will switch to
     * the structure at index `indexes.structure`, and if environments where
     * passed to the constructor and the current display mode is `'atom'`,
     * highlight the atom-centered environment corresponding to `indexes.atom`.
     *
     * @param  indexes         structure / atom pair to display
     * @param  selectedGUID   GUID of the widget to show
     */
    public show(indexes: Indexes, selectedGUID: string = this._active) {
        const data = this._viewers.get(selectedGUID);
        assert(data !== undefined);

        const widget = data.widget;
        if (data.current.structure !== indexes.structure) {
            const options = {
                packed: false,
                trajectory: true,
            } as any;
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
        this.active = selectedGUID;
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
                this.onselect(indexes, this._active);
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
                this.onselect(indexes, this._active);
                // continue playing until the advance callback returns false
                this.atomPlayback(advance);
            }
        }, parseFloat(this._delay.value) * 100);
    }

    /**
     * Remove all HTML added by this [[StructureViewer]] in the current document
     */
    public remove(): void {
        for (const data of this._viewers.values()) {
            data.widget.remove();
        }
        this._root.parentElement!.innerHTML = '';
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
     * Get the structure at the given index in a format JSmol can undertand
     * and load. [[Structure]] already rendered as strings are cached for faster
     * subsequent access.
     * @param  index index of the structure
     * @return       a string that can be passed to JSmol `load INLINE` command
     */
    private _structureForJSmol(index: number): string {
        if (this._cachedStructures[index] === undefined) {
            const s = this.loadStructure(index, this._structures[index]);
            const check = checkStructure(s);
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
     * Returns the GUID string corresponding to the selected widget.
     */
    public get active() {
        return this._active;
    }

    /*
     * Function to set the active widget for communicating with the map
     */
    public set active(activeGUID: string) {
        /// this is here to prevent infinite loops
        if (activeGUID !== this._active) {
            let indexes;

            if (this._viewers.has(activeGUID) ) {
                if (this._viewers.has(this._active)) {
                    const oldButton = getByID(`chsp-activate-${this._active}`);
                    oldButton.classList.toggle('chsp-active-structure', false);
                    const oldTooltip = oldButton.parentElement!.querySelector('.chsp-tooltip');
                    assert(oldTooltip !== null);
                    oldTooltip.innerHTML = 'choose as active';

                    const oldCell = getByID(`gi-${this._active}`);
                    oldCell.classList.toggle('chsp-structure-viewer-cell-active', false);
                }
                this._active = activeGUID;
                const newButton = getByID(`chsp-activate-${this._active}`);
                newButton.classList.toggle('chsp-active-structure', true);
                const tooltip = newButton.parentElement!.querySelector('.chsp-tooltip');
                assert(tooltip !== null);
                tooltip.innerHTML = 'this is the active viewer';

                const newCell = getByID(`gi-${this._active}`);
                newCell.classList.toggle('chsp-structure-viewer-cell-active', true);

                const activeWidgetData = this._viewers.get(this._active);
                assert(activeWidgetData !== undefined);

                activeWidgetData.widget.onselect = (atom: number) => {
                    if (this._indexer.mode !== 'atom') {
                        return;
                    }

                    activeWidgetData.widget.highlight(atom);

                    // if the viewer is showing a bigger supercell than [1, 1, 1], the
                    // atom index can be outside of [0, natoms), so make sure it is
                    // inside this range.
                    const current = activeWidgetData.current;
                    const atom_id = atom % activeWidgetData.widget.natoms()!;
                    indexes = this._indexer.from_structure_atom(current.structure, atom_id);
                    this.onselect(indexes, this._active);
                };

                this._delay = getByID<HTMLInputElement>(`chsp-${this._active}-playback-delay`);

                if (this._indexer.mode === 'structure') {
                    indexes = this._indexer.from_structure_atom(activeWidgetData.current.structure);
                } else {
                    const structure = activeWidgetData.current.structure;
                    const atom = activeWidgetData.current.atom;
                    indexes = this._indexer.from_structure_atom(structure, atom);
                }
                this.onselect(indexes, this._active);
            }
        }
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

    /*
     * Function to setup the cell in the structure viewer grid.
     * Will generate a GUID string if one does not exist for the cell
     * and instantiate all necessary buttons.
     */
    private _setupCell(cellGUID: string, colNum: number, rowNum: number): string {
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
            const isActive = template.content.firstChild! as HTMLElement;
            isActive.onclick = () => {this.active = cellGUID; };
            cell.appendChild(isActive);

            // add a button to close the widget
            template.innerHTML = `<button
                class="btn btn-light btn-sm chsp-has-tooltip chsp-viewer-button"
                style="top: 8px; right: 40px;">
                    <span>${CLOSE_SVG}</span>
                    <span class="chsp-tooltip">remove viewer</span>
                </button>`;
            const close = template.content.firstChild! as HTMLElement;
            close.onclick = () => {
                this._removeWidget(cellGUID);
                this._setupGrid(this._viewers.size);
            };
            cell.appendChild(close);

            // add a button to duplicate the widget
            template.innerHTML = `<button
                class="btn btn-light btn-sm chsp-has-tooltip chsp-viewer-button"
                style="top: 8px; right: 70px;">
                    <span>${DUPLICATE_SVG}</span>
                    <span class="chsp-tooltip">duplicate viewer</span>
                </button>`;
            const duplicate = template.content.firstChild! as HTMLElement;

            duplicate.onclick = () => {
                const data = this._viewers.get(cellGUID);
                assert(data !== undefined);

                let index;
                if (this._indexer.mode === 'structure') {
                    index = this._indexer.from_structure_atom(data.current.structure);
                } else {
                    index = this._indexer.from_structure_atom(data.current.structure, data.current.atom);
                }

                const newGuid = this._setupGrid(this._viewers.size + 1);
                assert(newGuid.length === 1);
                this.show(index, newGuid[0]);
                this.onselect(index, newGuid[0]);
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

    /*
     * Function to initialize the grid instance for `this._nwidgets` cells and
     * place onto the DOM element mapped in `this._root`. If more cells are
     * needed, this function return the list of new cell GUID
     */
    private _setupGrid(nwidgets: number): string[] {
        const newGUID = [] as string[];
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
            let cellGUID;
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

                const current = {atom: -1, structure: -1, environment: -1};
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

    /*
     * Removes a widget from the structure viewer grid.
     * The parameter force pertains to removing the *only* widget in the grid,
     * which should only be done when changing datasets.
     */
    private _removeWidget(trashedGUID: string, force: boolean = false) {
        if (this._viewers.size > 1 || force === true) {
            const widgetRoot = getByID(`chsp-${trashedGUID}`);

            this.onselect({structure: -1, environment: -1}, trashedGUID);
            if ( widgetRoot.parentNode !== null) {
                widgetRoot.parentNode.removeChild(widgetRoot);
            }

            const deadCell = getByID(`gi-${trashedGUID}`);
            if (deadCell !== null) {deadCell.remove(); }

            this._viewers.delete(trashedGUID);

            if (this._active === trashedGUID ) {
                if (this._viewers.size > 0) {
                  this.active = this._viewers.keys().next().value;
                } else { this._active = ''; }
            }
        }
    }
}
