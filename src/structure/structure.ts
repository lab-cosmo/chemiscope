/**
 * @packageDocumentation
 * @module structure
 */

import assert from 'assert';

import {Environment, isPositiveInteger, Structure, UserStructure} from '../dataset';
import {EnvironmentIndexer, generateGUID, getByID, getNextColor, Indexes, sendWarning} from '../utils';

import {structure2JSmol} from './jsmol';
import {JSmolWidget} from './widget';

import close_svg from './close.svg';
import duplicate_svg from './duplicate.svg';
import HTML_SETTINGS from './settings.html';

const MAX_WIDGETS = 9;

/*
 * Function to return *optimal* arrangement of n widgets.
 * Defaults to a 4 x ? grid. Hardcoded for simplicity.
 */
function bestArrangement(n: number) {
  switch (n) {
    case 1:
    case 2:
      return {rows: n, columns: 1};
      break;
    case 3:
    case 4:
      return {rows: 2, columns: 2};
      break;
    case 5:
    case 6:
      return {rows: 3, columns: 2};
      break;
    case 8:
    case 7:
    case 9:
      return {rows: 3, columns: 3};
      break;
    case 10:
    case 11:
    case 12:
      return {rows: 4, columns: 3};
      break;
    default:
      return {rows: 4, columns: 4};
      break;

  }
}

/*
 * Needs documentation.
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

    // GUID of the Active Widget
    private _active: string;
    // Map of Widgets GUIDS to their color, widget, and current structure
    private _widgetMap: Map<string, WidgetGridData>;
    // List of Widgets GUIDS in the StructureViewer
    private _GUIDs: string[];
    // Documentation needed.
    private _root: HTMLElement;
    /// Playback delay setting
    private _delay: HTMLInputElement;
    /// List of structures in the dataset
    private _structures: Structure[] | UserStructure[];
    /// Cached string representation of structures
    private _cachedStructures: string[];
    /// Optional list of environments for each structure
    private _environments?: Environment[][];
    // Documentation needed.
    private _indexer: EnvironmentIndexer;
    // path to the j2s files used by JSmol.
    // saved for instantiating new Widget instances
    private _j2spath: string;

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
        this._widgetMap = new Map();
        // this._widgetMap = {};
        this._GUIDs = [];
        this.onselect = () => {};

        // Initialize the _root div as a grid (see function below)
        const root = getByID(id);

        // Initialize the grid root
        this._root = document.createElement('div');
        this._root.id = 'grid-root';
        this._root.className = 'chsp-structure-viewer-grid';

        root.appendChild(this._root);

        this._setupGrid(1);
        this.active = this._GUIDs[0];
        this._active = this._GUIDs[0];
        this._delay = this._setDelay();

    }
    /**
     * Show a new structure, as identified by `indexes`. This will switch to
     * the structure at index `indexes.structure`, and if environments where
     * passed to the constructor and the current display mode is `'atom'`,
     * highlight the atom-centered environment corresponding to `indexes.atom`.
     *
     * @param  indexes         structure / atom pair to display
     * @param  wi              GUID of the widget to show
     */
    public show(indexes: Indexes, selectedGUID: string = this._active) {
        const widgetData = this._widgetMap.get(selectedGUID);
        if (widgetData !== undefined) {
          const widget = widgetData.widget;
          if (widgetData.current.structure !== indexes.structure) {
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
              if (widgetData.current.atom !== indexes.atom) {
                  widget.highlight(indexes.atom);
              }
          } else {
              widget.highlight(undefined);
          }

          widgetData.current = indexes;
          this.active = selectedGUID;
      }
    }
    /**
     * Register a `callback` to compute the placement of the settings modal.
     *
     * The callback gets the current placement of the settings as a
     * [DOMRect](https://developer.mozilla.org/en-US/docs/Web/API/DOMRect),
     * and should return top and left positions in pixels, used with `position:
     * fixed`. The callback is called once, the first time the settings are
     * opened.
     */
    public settingsPlacement(callback: (rect: DOMRect) => {top: number, left: number}) {
        for (const guid of this._GUIDs) {
          const widgetData = this._widgetMap.get(guid);
          if (widgetData !== undefined ) {
            widgetData.widget.settingsPlacement(callback);
          }
        }
    }
    /**
     * Start playing the trajectory of structures in this dataset, until
     * `advance` returns false
     */
    public structurePlayback(advance: () => boolean) {
        setTimeout(() => {
            if (advance()) {
                const widgetData = this._widgetMap.get(this._active);
                if (widgetData !== undefined) {
                  const current = widgetData.current;
                  const structure = (current.structure + 1) % this._indexer.structuresCount();
                  const indexes = this._indexer.from_structure_atom(structure, 0);
                  this.show(indexes);
                  this.onselect(indexes, this._active);
                  // continue playing until the advance callback returns false
                  this.structurePlayback(advance);
                }
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
                const widgetData = this._widgetMap.get(this._active);
                if (widgetData !== undefined) {
                  const current = widgetData.current;
                  const structure = current.structure;
                  const atom = (current.atom! + 1) % this._indexer.atomsCount(structure);
                  const indexes = this._indexer.from_structure_atom(structure, atom);
                  this.show(indexes);
                  this.onselect(indexes, this._active);
                  // continue playing until the advance callback returns false
                  this.atomPlayback(advance);
                }
            }
        }, parseFloat(this._delay.value) * 100);
    }

    /**
     * Remove all HTML added by this [[StructureViewer]] in the current document
     */
    public remove(): void {
        for (const data of this._widgetMap.values()) {
            data.widget.remove();
        }
        this._root.parentElement!.innerHTML = '';
    }

    /**
     * Needs documentation
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
    public get active() {return this._active; }
    /*
     * Function to set the active widget for communicating with the map
     */
    public set active(activeGUID: string) {

      /// this is here to prevent infinite loops
      if (activeGUID !== this._active) {
        const newActiveIdx = this._GUIDs.indexOf(activeGUID);
        let indexes;

        if (newActiveIdx >= 0) {
          if (this._active !== '' && this._active !== undefined) {
            const oldButton = getByID(`chsp-activate-${this._active}`);
            oldButton.classList.toggle('chsp-inactive-structure-marker', true);
            oldButton.classList.toggle('chsp-active-structure-marker', false);
            oldButton.innerHTML = `<span class="tooltiptext">Choose as active.</span>`;
          }
          this._active = activeGUID;
          const newButton = getByID(`chsp-activate-${this._active}`);
          newButton.classList.toggle('chsp-inactive-structure-marker', false);
          newButton.classList.toggle('chsp-active-structure-marker', true);
          newButton.innerHTML = `<span class="tooltiptext">This is the active button</span>`;

          const activeWidgetData = this._widgetMap.get(this._active);
          if (activeWidgetData !== undefined) {
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
          }
        }
        const widgetData = this._widgetMap.get(this._active);
        if (widgetData !== undefined) {
          if (this._indexer.mode === 'structure') {
            indexes = this._indexer.from_structure_atom(widgetData.current.structure);
            } else {
            indexes = this._indexer.from_structure_atom(widgetData.current.structure,
                                                        widgetData.current.atom);
            }
          this.onselect(indexes, this._active);
        }
      }
    }
    /*
     * Needs documentation
     */
    private _setDelay(): HTMLInputElement {
        return getByID(`chsp-${this._active}-playback-delay`) as HTMLInputElement;
    }
    /*
     * Function to resize the grid to this._GUIDs.length + inc widgets
     */
    private _resizeGridInc(inc: number) {
      if (this._GUIDs.length + inc < 1) {
        sendWarning('Cannot delete last widget.');
      } else {
        if (this._GUIDs.length + inc > MAX_WIDGETS ) {
        sendWarning(`Widget grid cannot contain more than ${MAX_WIDGETS} widgets.`);
        } else {
        this._setupGrid(this._GUIDs.length + inc);
        }
      }
    }
    /*
     * Function to setup the cell in the structure viewer grid.
     * Will generate a GUID string if one does not exist for the cell
     * and instantiate all necessary buttons.
     */
    private _setupCell(cellNo: number, colNum: number, rowNum: number) {
      if ( cellNo >= this._GUIDs.length) {
        this._GUIDs.push(generateGUID());
      }

      const cellGUID = this._GUIDs[cellNo];
      const cellId = `gi-${cellGUID}`;
      let cell = document.getElementById(cellId);

      if (cell === null) {
        cell = document.createElement('div');
        cell.id = cellId;
        cell.classList.add('chsp-structure-viewer-cell', 'grid-item');
        cell.style.gridColumn = `${colNum}`;
        cell.style.gridRow = `${rowNum}`;

        // Do not set `cell.onclick = () => {this.active = GUID}`
        // because it will conflict with the button behavior

        // add a button to activate the widget (i.e. set as `active`)
        const activeFlag = document.createElement('button');
        activeFlag.classList.add('chsp-inactive-structure-marker');
        activeFlag.id = `chsp-activate-${cellGUID}`;
        const colors = [];
        for (const guid of this._GUIDs) {
          const widgetData = this._widgetMap.get(guid);
          if (widgetData !== undefined) {
            colors.push(widgetData.color);
          }
        }
        const color = getNextColor(colors);
        activeFlag.style.backgroundColor = color;
        activeFlag.onclick = () => {this.active = cellGUID; };
        activeFlag.innerHTML = `<span class="tooltiptext">Choose as active</span>`;
        cell.appendChild(activeFlag);

        // add a button to close the widget
        const close = document.createElement('button');
        close.classList.add('chsp-close-widget-button', 'btn', 'btn-light', 'btn-sm');
        close.id = `chsp-close-widget-button-${cellGUID}`;
        close.onclick = () => {this._removeWidget(cellGUID); this._resizeGridInc(0); };
        close.innerHTML = `<span class="tooltiptext">Close widget</span><object>${close_svg}</object>`;
        cell.appendChild(close);

        // add a button to duplicate the widget
        const duplicate = document.createElement('button');
        duplicate.classList.add('chsp-duplicate-widget-button', 'btn', 'btn-light', 'btn-sm');
        duplicate.id = `chsp-duplicate-widget-button-${cellGUID}`;
        duplicate.onclick = () => {
                const widgetData = this._widgetMap.get(cellGUID);
                if (widgetData !== undefined) {
                  this._resizeGridInc(1);
                  const myCurrent = widgetData.current;
                  let index;
                  if (this._indexer.mode === 'structure') {
                    index = this._indexer.from_structure_atom(myCurrent.structure);
                    } else {
                    index = this._indexer.from_structure_atom(myCurrent.structure,
                                                              myCurrent.atom);
                    }
                  this.show(index, this._GUIDs[this._GUIDs.length - 1]);
                  this.onselect(index, this._GUIDs[this._GUIDs.length - 1]);
                }
                };
        duplicate.innerHTML = `<span class="tooltiptext">Add duplicate widget</span><object>${duplicate_svg}</object>`;
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
     * Function to initialize the grid instance for `this._nwidgets` cells and place
     * onto the DOM element mapped in `this._root`
     */
    private _setupGrid(nwidgets: number) {

      if (nwidgets <= 0) {
        throw Error('Cannot initialize a grid with 0 widgets.');
      }

      // Determine best arrangement for nwidgets
      const arrangement = bestArrangement(nwidgets);

      if (this._GUIDs.length > nwidgets) {
        sendWarning(`Warning: Eliminating last ${this._GUIDs.length - nwidgets} widgets.`);
        const wl = this._GUIDs.length;
        for (let i = nwidgets; i < wl; i++) {
          this._removeWidget(this._GUIDs[i]);
        }
      }

      // Start at the third row to skip the inc/dec buttons
      let rowNum = 1;
      let colNum = 1;

      for (let c = 0; c < nwidgets; c++) {
        let color = this._setupCell(c, colNum, rowNum);
        if (color === '') {
          const colors = [];
          for (const guid of this._GUIDs) {
            const widgetData = this._widgetMap.get(guid);
            if (widgetData !== undefined) {
              colors.push(widgetData.color);
            }
          }
          color = getNextColor(colors);
        }

        colNum++;
        if (colNum > arrangement.columns) {
            rowNum++;
            colNum = 1;
          }

        // add a new widget if necessary
        const widgetData = this._widgetMap.get(this._GUIDs[c]);
        if (widgetData === undefined) {
          const widget = new JSmolWidget(`gi-${this._GUIDs[c]}`,
                                             this._j2spath,
                                             'chsp-' + this._GUIDs[c]);
          const current = {atom: -1, structure: -1, environment: -1};
          this._widgetMap.set(this._GUIDs[c], {color: color,
                                             current: current,
                                             widget: widget,
                                            });
          const index = {atom: 0, structure: 0, environment: 0};
          this.show(index, this._GUIDs[c]);
          this.onselect(index, this._GUIDs[c]);
        } else {
          widgetData.widget.script('refresh');
        }
      }
    }
    /*
     * Removes a widget from the structure viewer grid.
     * The parameter force pertains to removing the *only* widget in the grid,
     * which should only be done when changing datasets.
     */
    private _removeWidget(trashedGUID: string, force: boolean = false) {

      if ( this._GUIDs.length > 1 || force === true) {
        const widgetRoot = getByID(`chsp-${trashedGUID}`);

        this.onselect({structure: -1, environment: -1}, trashedGUID);
        if ( widgetRoot.parentNode !== null) {
          widgetRoot.parentNode.removeChild(widgetRoot);
        }

        const deadCell = getByID(`gi-${trashedGUID}`);
        if (deadCell !== null) {deadCell.remove(); }

        const cellNo = this._GUIDs.indexOf(trashedGUID);
        this._widgetMap.delete(trashedGUID);
        this._GUIDs.splice(cellNo, 1);

        if (this._active === trashedGUID ) {
          if ( this._widgetMap.get(this._GUIDs[0]) !== undefined) {
            this.active = this._GUIDs[0];
          } else {
            this._active = '';
          }
        }
      }
    }

}
