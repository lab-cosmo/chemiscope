/**
 * @packageDocumentation
 * @module structure
 */

import assert from 'assert';

import { JSmolApplet, JmolObject } from 'jsmol';

import { SavedSettings } from '../options';
import { generateGUID, getByID } from '../utils';
import { PositioningCallback } from '../utils';

import { StructureOptions } from './options';

require('../static/chemiscope.css');

function invertedZoom(): boolean {
    // both safari and chrome include `Safari/` in their user agent, Firefox
    // does not
    const safariOrChrome = window.navigator.userAgent.includes('Safari/');
    // TODO: update this condition when ARM mac comes out
    if (safariOrChrome && window.navigator.platform === 'MacIntel') {
        return true;
    } else {
        return false;
    }
}

/** @hidden
 * Create a stylesheet in the main `document` with the given `rules`
 */
function createStyleSheet(rules: string[]): CSSStyleSheet {
    const style = document.createElement('style');
    style.type = 'text/css';
    document.head.appendChild(style);
    const sheet = style.sheet as CSSStyleSheet;
    for (const rule of rules) {
        sheet.insertRule(rule);
    }
    return sheet;
}

/**
 * @hidden
 * Get the "555" style supercell description. This allows to put the extra
 * cells both on negative and positive sides.
 */
function supercell_555(supercell: [number, number, number]): string {
    const central = [
        Math.floor(supercell[0] / 2),
        Math.floor(supercell[1] / 2),
        Math.floor(supercell[2] / 2),
    ];

    const first = [5 - central[0], 5 - central[1], 5 - central[2]];

    const second = [
        5 + supercell[0] - (central[0] + 1),
        5 + supercell[1] - (central[1] + 1),
        5 + supercell[2] - (central[2] + 1),
    ];

    return `{${first[0]}${first[1]}${first[2]} ${second[0]}${second[1]}${second[2]} 1}`;
}

/** A spherical atom-centered environment */
export interface Environment {
    /**
     * Cutoff radius of the environment.
     *
     * Atoms inside the sphere centered on the `center` atom with this radius
     * are part of the environment.
     */
    cutoff: number;
}

/** Position of the setting in the page, using 'position: fixed'. */
export interface Position {
    /** number of pixels from the top of the page */
    top: number;
    /** number of pixels from the left of the page */
    left: number;
}

/** Possible options passed to `JSmolWidget.load` */
export interface LoadOptions {
    /** Supercell to display (default: [1, 1, 1]) */
    supercell: [number, number, number];
    /** Should we display a packed cell (default: false) */
    packed: boolean;
    /** Should preserve we the current camera orientation (default: false) */
    keepOrientation: boolean;
    /** Are we loading a file part of a trajectory (default: false) */
    trajectory: boolean;
    /** List of atom-centered environments */
    environments: Environment[];
    /**
     * Index of the environment to highlight, this is only considered if
     * `environments` is set.
     */
    highlight: number;
}

const DEFAULT_SERVER_URL = 'https://chemapps.stolaf.edu/jmol/jsmol/php/jsmol.php';

/**
 * A light wrapper around JSmol, displaying a single structure. This is kept as
 * a separate class from the [[ViewersGrid]] to be easier to use outside of
 * chemiscope.
 */
export class JSmolWidget {
    /** callback called when a new atom is clicked on */
    public onselect: (atom: number) => void;
    /**
     * Callback to get the initial positioning of the settings modal.
     *
     * The callback is called once, the first time the settings are opened.
     */
    public positionSettingsModal: PositioningCallback;

    /**
     * Unique identifier of this viewer.
     *
     * All HTML elements created by this class use this ID to ensure unicity.
     */
    public guid: string;

    /// The HTML element serving as root element for the viewer
    private _root: HTMLElement;
    /// Reference to the global Jmol variable
    private _Jmol: JmolObject;
    /// Reference to the JSmol applet to be updated with script calls
    private _applet: JSmolApplet;
    /// Representation options from the HTML side
    private _options: StructureOptions;
    /// The supercell used to intialize the viewer
    private _initialSupercell?: [number, number, number];
    // button to reset the environment cutoff to its original value
    private _resetEnvCutof!: HTMLButtonElement;
    // button to reset reset the supercell
    private _resetSupercell!: HTMLButtonElement;
    /// Show some information on the currently displayed cell to the user
    private _cellInfo: HTMLElement;
    /// Dynamic CSS used to hide options related to the unit cell if there is no
    /// unit cell in the current structure
    private _noCellStyle: CSSStyleSheet;
    /// Dynamic CSS used to hide options related to environments if there are
    /// environments for the current structure
    private _noEnvsStyle: CSSStyleSheet;
    /// callback called on successfull (re)loading of a file
    private _loadedCallback?: () => void;
    /// Number of atoms in the currenly loaded structure 1x1x1 supercell.
    ///
    /// If another supercell is displayed, there are na x nb x nc x natoms atoms
    /// in the displayed structure
    ///
    /// This is set to `undefined` until a structure is actually loaded
    private _natoms?: number;
    /// List of atom-centered environments for the current structure
    private _environments?: Environment[];
    /// Index of the central atom to highlight, if any. This requires
    /// `_environments` to be defined.
    private _highlighted?: number;
    /// caching the last selected atom to prevent calling onselect multiple
    /// time with the same atom number
    private _lastSelected: number;

    /**
     * Create a new JSmolWidget inside the HTML DOM element with the given `id`.
     *
     * @param id HTML element id inside which the viewer will be created
     * @param j2sPath path where j2s files can be loaded by Jmol
     * @param guid (optional) unique identifier for the widget
     * @param serverURL URL where to find `jsmol.php`
     */
    constructor(
        id: string,
        j2sPath: string,
        guid?: string,
        serverURL: string = DEFAULT_SERVER_URL
    ) {
        if (window.Jmol === undefined) {
            throw Error('Jmol is required, load it from your favorite source');
        }

        if (guid === undefined) {
            guid = generateGUID();
        }

        // add a 'chsp-' prefic to ensure the id start with letter. It looks like
        // if the id start with a number (2134950-ffff-4879-82d8-5c9f81dd00ab)
        // then bootstrap code linking modal button to the modal fails ¯\_(ツ)_/¯
        this.guid = 'chsp-' + guid;
        this._lastSelected = -1;

        // store a reference to the global Jmol and the HTML root
        this._Jmol = window.Jmol;
        // Do not insert new applets automatically
        this._Jmol.setDocument(false);

        this._root = document.createElement('div');
        const root = getByID(id);
        root.appendChild(this._root);

        this._root.style.position = 'relative';
        this._root.id = this.guid;
        this._root.style.width = '100%';
        this._root.style.height = '100%';

        this._cellInfo = document.createElement('span');
        this._cellInfo.classList.add(
            'chsp-cell-info',
            'chsp-hide-if-no-cell',
            'badge',
            'badge-light'
        );
        this._root.appendChild(this._cellInfo);

        this._options = new StructureOptions(this._root, this.guid, (rect) =>
            this.positionSettingsModal(rect)
        );
        this._connectSettings();

        this._applet = this._createApplet(j2sPath, serverURL);
        if (!invertedZoom()) {
            // Invert (cf -1 below) wheel zoom direction to match the one in the
            // map. _DELTAY is replaced by 1 or -1 depending on the wheel/scroll
            // direction.
            //
            // For some reason, safari & chrome on OsX are already reversed,
            // and binding a callback to WHEEL though jsmol is very slow, so
            // don't reverse in this case
            this.script('bind "WHEEL" "zoom *@{1.15 ** (-1 * _DELTAY)}";');
        }

        this._loadedCallback = undefined;
        // create a global function with unique name and install it as callback
        // for JSmol. This function will then call the callback inside this
        // instance of the viewer
        const name: string = this.guid.replace(/-/g, '_') + '_loaded_callback';
        const window_as_map = (window as unknown) as Record<string, unknown>;

        window_as_map[name] = (
            _a: unknown,
            _b: unknown,
            _c: unknown,
            _d: unknown,
            _e: unknown,
            status: unknown
        ) => {
            // eslint-disable-next-line @typescript-eslint/no-unsafe-call, @typescript-eslint/no-explicit-any, @typescript-eslint/no-unsafe-member-access
            if ((status as any).valueOf() === 3 && this._loadedCallback !== undefined) {
                this._loadedCallback();
            }
        };
        this.script(`set LoadStructCallback "${name}"`);

        this._noCellStyle = createStyleSheet([
            `#${this.guid} .chsp-hide-if-no-cell { display: none; }`,
            `#${this.guid}-settings .chsp-hide-if-no-cell { display: none; }`,
        ]);

        this._noEnvsStyle = createStyleSheet([
            `#${this.guid}-settings .chsp-hide-if-no-environments { display: none; }`,
        ]);

        // By default, position the modal for settings on top of the viewer,
        // centered horizontally
        this.positionSettingsModal = (rect: DOMRect) => {
            const rootRect = this._root.getBoundingClientRect();
            return {
                left: rootRect.left + rootRect.width / 2 - rect.width / 2,
                top: rootRect.top + 20,
            };
        };

        this.onselect = () => {};
        this._root.onclick = () => {
            // picked is 0-based
            const picked = this.evaluate('_atomPicked');

            if (picked === '') {
                return;
            } else {
                assert(typeof picked === 'number');
            }

            if (picked === this._lastSelected) {
                return;
            } else {
                this._lastSelected = picked;
            }

            this.onselect(picked);
        };
    }

    /**
     * Remove all HTML added by this [[JSmolWidget]] in the current document
     */
    public remove(): void {
        if (this._root.parentElement !== null) {
            this._root.parentElement.innerHTML = '';
        }
        this._options.remove();
    }

    /**
     * Get the number of atoms in the structure, or `undefined` if no structure
     * is currenly loaded
     *
     * @return the number of atoms in the currenly loaded structure
     */
    public natoms(): number | undefined {
        return this._natoms;
    }

    /**
     * Run the given `commands` for this viewer. `load` commands should use the
     * [[JSmolWidget.load|corresponding function]].
     */
    public script(commands: string): void {
        if (commands.includes('load ')) {
            throw Error('Do not use JSmolWidget.script to load a structure, but JSmolWidget.load');
        }
        this._Jmol.script(this._applet, commands);
    }

    /**
     * Evaluate the given commands using JSmol and return the corresponding
     * value. This calls `Jmol.evaluateVar` behind the scenes.
     */
    public evaluate(commands: string): unknown {
        return this._Jmol.evaluateVar(this._applet, commands);
    }

    /**
     * Highlight a given `atom` in the current structure.
     *
     * If a supercell larger than [1, 1, 1] is currently displayed, this
     * function accept indexes larger than the result of `natoms()`, and will
     * then highlight atoms outside of the central image.
     *
     * @param atom index of the central atom in the environment to show,
     *             or `undefined` to disable highlighting.
     */
    public highlight(environment?: number): void {
        this._changeHighlighted(environment);
        this._updateState();
    }

    /**
     * Load `data` in this viewer.
     *
     * This function can take any value used in a `load XXX` script in Jmol:
     * inline data, path to a file, $\<molecules\>, etc.
     *
     * @param data molecule to load, in a format supported by Jmol [load
     *             command](https://chemapps.stolaf.edu/jmol/docs/#load).
     * @param options options for the new structure
     */
    public load(data: string, options: Partial<LoadOptions> = {}): void {
        if (data === "''" || data === '""') {
            throw Error('invalid use of JSmolWidget.load to reload data');
        }

        this._natoms = undefined;

        let supercell: [number, number, number];
        if (options.supercell === undefined) {
            // keep pre-existing supercell settings, default to [1, 1, 1] from
            // settings.html
            supercell = [
                this._options.supercell[0].value,
                this._options.supercell[1].value,
                this._options.supercell[2].value,
            ];
        } else {
            supercell = options.supercell;
        }

        if (options.packed !== undefined) {
            this._options.packedCell.value = options.packed;
        }

        if (options === undefined) {
            this._environments = undefined;
        } else {
            this._environments = options.environments;
            if (this._environments !== undefined) {
                this._options.environments.activated.value = true;
            }
        }

        let keepOrientation: boolean;
        if (options.keepOrientation === undefined) {
            // keep pre-existting settings if any
            keepOrientation = this._options.keepOrientation.value;
        } else {
            keepOrientation = options.keepOrientation;
        }

        const trajectoryOptions = getByID(`${this.guid}-trajectory-settings-group`);
        if (options.trajectory === undefined || !options.trajectory) {
            trajectoryOptions.style.display = 'none';
        } else {
            trajectoryOptions.style.display = 'block';
        }

        if (this._environments === undefined) {
            this._noEnvsStyle.disabled = false;
            this._changeHighlighted(undefined);
        } else {
            if (options.packed !== undefined && options.packed) {
                throw Error('Can not have both packed cell and environments');
            }

            this._noEnvsStyle.disabled = true;
            if (options.highlight !== undefined) {
                this._changeHighlighted(options.highlight);
            } else {
                this._changeHighlighted(0);
            }
        }

        const [a, b, c] = supercell;
        if (this._initialSupercell === undefined) {
            this._initialSupercell = supercell;
            this._resetSupercell.innerHTML = `reset ${a}x${b}x${c} supercell`;
        }
        this._options.supercell[0].value = a;
        this._options.supercell[1].value = b;
        this._options.supercell[2].value = c;
        this._setCellInfo();

        /// JSmol overwite the Error object, using `.caller` to get stacktraces.
        /// This fails in modern browsers if any code other code try to throw an
        /// error, since the JSmol Error object will try to get the stack here
        /// too. So we save and reset the global Error object when loading a new
        /// structure.
        const savedError = window.Error;

        this._loadedCallback = () => {
            window.Error = savedError;

            if (this.evaluate('@{unitcell("conventional")}') === 'ERROR') {
                this._noCellStyle.disabled = false;
            } else {
                this._noCellStyle.disabled = true;
            }

            const repeat =
                this._options.supercell[0].value *
                this._options.supercell[1].value *
                this._options.supercell[2].value;
            this._natoms = parseInt(this.evaluate('{*}.size') as string, 10) / repeat;
            if (this._environments !== undefined && this._environments.length !== this._natoms) {
                let message = 'invalid number of environments for this structure ';
                message += `got ${this._environments.length} environments, there are ${this._natoms} atoms`;
                // We can not throw an error here, since it seems to be caught
                // by JSmol.
                // eslint-disable-next-line no-console
                console.error(message);
            }

            // once we are done, disable the callback to prevent it from firing
            // on restyle/reload
            this._loadedCallback = undefined;

            // make sure any change to the setting before the structure is fully
            // loaded are applied
            this._updateState();
        };

        this._do_load(data, keepOrientation);
    }

    /**
     * Applies saved settings, possibly filling in with default values
     */
    public applySettings(settings: SavedSettings): void {
        this._options.applySettings(settings);
    }

    /**
     * Save the values of the current settings in a way that an be used with
     * [[applySettings]] or saved to JSON.
     */
    public saveSettings(): SavedSettings {
        return this._options.saveSettings();
    }

    private _reload() {
        if (!this._loaded()) {
            return;
        }
        this._do_load('""', true);
    }

    private _do_load(data: string, keepOrientation: boolean) {
        if (data.includes(';')) {
            throw Error("invalid ';' in  JSmolWidget.load");
        }
        if (data.includes('packed')) {
            throw Error("invalid 'packed' in  JSmolWidget.load");
        }

        const packed = this._options.packedCell.value ? ' packed' : '';

        const saveOrientation = keepOrientation ? `save orientation "${this.guid}"` : '';
        const restoreOrientation = keepOrientation ? `restore orientation "${this.guid}"` : '';

        this._setCellInfo();

        const supercell = supercell_555([
            this._options.supercell[0].value,
            this._options.supercell[1].value,
            this._options.supercell[2].value,
        ]);

        const commands = `
            ${saveOrientation};
            load ${data} ${supercell} ${packed};
            ${this._updateStateCommands()};
            ${restoreOrientation};
        `;
        this._Jmol.script(this._applet, commands);
    }

    private _changeHighlighted(environment?: number) {
        if (environment !== undefined) {
            const repeat =
                this._options.supercell[0].value *
                this._options.supercell[1].value *
                this._options.supercell[2].value;
            if (this._natoms !== undefined && environment >= this._natoms * repeat) {
                let message = 'selected environment is out of bounds: ';
                message += `got ${environment}, we have ${
                    this._natoms * repeat
                } atoms in the current structure`;
                throw Error(message);
            }

            if (this._environments === undefined) {
                throw Error('no environments defined for the current structure');
            }
        }

        this._highlighted = environment;

        if (this._highlighted === undefined) {
            this._enableEnvironmentSettings(false);
            this._options.environments.cutoff.value = 0;
        } else {
            this._enableEnvironmentSettings(true);

            // keep user defined cutoff, if any
            if (this._options.environments.cutoff.value <= 0) {
                this._options.environments.cutoff.value = this._currentDefaultCutoff();
            }
        }
    }

    private _createApplet(j2sPath: string, serverURL: string): JSmolApplet {
        const width = this._root.clientWidth;
        const height = this._root.clientHeight;
        if (width === 0 || height === 0) {
            assert(this._root.parentElement !== null);
            const parentId = this._root.parentElement.id;
            // eslint-disable-next-line no-console
            console.error(
                `JSmolWidget: width (=${width}px) or heigh (=${height}px) ` +
                    `of #${parentId} is zero, you will not see the molecules`
            );
        }

        const INITIAL_SCRIPT = `
            // use anti-aliasing
            set antialiasdisplay;
            // remove jmol logo
            set frank off;
            // use the smaller of height/width when setting zoom level
            set zoomlarge false;
            // Allow sending script commands while moveto is executing
            set waitformoveto off;
            hide off;
        `;
        // create the main Jmol applet
        const applet = this._Jmol.getApplet(this.guid, {
            height: '100%',
            width: '100%',

            disableInitialConsole: true,
            disableJ2SLoadMonitor: true,
            j2sPath: j2sPath,
            script: INITIAL_SCRIPT,
            serverURL: serverURL,
            use: 'HTML5',
            zIndexBase: 1,
        });

        const div = document.createElement('div');
        div.style.height = '100%';
        div.style.width = '100%';
        this._root.appendChild(div);
        div.innerHTML = this._Jmol.getAppletHtml(applet);
        // Jmol rely on this script being implicitly executed, but this is not
        // the case when using innerHTML (compared to jquery .html()). So let's
        // manually execute it
        applet._cover(false);

        return applet;
    }

    private _connectSettings() {
        // recursively bind the right update function to HTMLSetting `onchange`
        const bindUpdateState = (object: Record<string, unknown>) => {
            for (const key in object) {
                if (key.startsWith('_')) {
                    continue;
                }
                const setting = object[key] as Record<string, unknown>;
                if ('onchange' in setting) {
                    // if any setting changes, update the full state of JSmol
                    // this is not the most efficient thing to do, but some
                    // settings have dependencies on one another.
                    setting.onchange = () => this._updateState();
                } else {
                    bindUpdateState(setting);
                }
            }
        };
        bindUpdateState((this._options as unknown) as Record<string, unknown>);

        // For changes to the cell, we have to reload the structure.
        // This override the function set above
        this._options.packedCell.onchange = () => this._reload();
        this._options.supercell[0].onchange = () => this._reload();
        this._options.supercell[1].onchange = () => this._reload();
        this._options.supercell[2].onchange = () => this._reload();

        // Deal with activation/de-activation of environments
        this._options.environments.activated.onchange = (value) => {
            this._enableEnvironmentSettings(value);
            this._updateState();
        };

        // Setup various buttons
        this._resetEnvCutof = getByID<HTMLButtonElement>(`${this.guid}-env-reset`);
        this._resetEnvCutof.onclick = () => {
            this._options.environments.cutoff.value = this._currentDefaultCutoff();
            this._updateState();
        };

        const alignX = getByID<HTMLButtonElement>(`${this.guid}-align-x`);
        alignX.onclick = () => this.script('moveto 1 axis x');

        const alignY = getByID<HTMLButtonElement>(`${this.guid}-align-y`);
        alignY.onclick = () => this.script('moveto 1 axis y');

        const alignZ = getByID<HTMLButtonElement>(`${this.guid}-align-z`);
        alignZ.onclick = () => this.script('moveto 1 axis z');

        const alignA = getByID<HTMLButtonElement>(`${this.guid}-align-a`);
        alignA.onclick = () => this.script('moveto 1 axis a');

        const alignB = getByID<HTMLButtonElement>(`${this.guid}-align-b`);
        alignB.onclick = () => this.script('moveto 1 axis b');

        const alignC = getByID<HTMLButtonElement>(`${this.guid}-align-c`);
        alignC.onclick = () => this.script('moveto 1 axis c');

        this._resetSupercell = getByID<HTMLButtonElement>(`${this.guid}-reset-supercell`);
        this._resetSupercell.onclick = () => {
            assert(this._initialSupercell !== undefined);
            this._options.supercell[0].value = this._initialSupercell[0];
            this._options.supercell[1].value = this._initialSupercell[1];
            this._options.supercell[2].value = this._initialSupercell[2];
            this._reload();
        };
    }

    private _updateState() {
        if (!this._loaded()) {
            return;
        }
        this.script(this._updateStateCommands());
    }

    private _showAxesCommands(axes: string): string {
        switch (axes) {
            case 'xyz':
                return `
                axes off;
                draw xaxis '>X' vector {0 0 0} {2 0 0} color red width 0.15;
                draw yaxis '>Y' vector {0 0 0} {0 2 0} color green width 0.15;
                draw zaxis '>Z' vector {0 0 0} {0 0 2} color blue width 0.15;
            `;
            case 'abc':
                return `
                draw xaxis delete; draw yaxis delete; draw zaxis delete;
                set axesUnitcell ON; axes 5;
            `;
            case 'off':
                return `
                draw xaxis delete; draw yaxis delete; draw zaxis delete;
                axes off;
            `;
            default:
                throw Error(`unkown axes selected: '${axes}'`);
        }
    }

    /// Get the commands to run to update the visualization state
    private _updateStateCommands(): string {
        const settings = this._options;
        const wireframe = `wireframe ${settings.bonds.value ? '0.15' : 'off'}`;
        const spacefill = `spacefill ${settings.spaceFilling.value ? '80%' : '23%'}`;

        let commands = '';
        if (this._highlighted === undefined || !settings.environments.activated.value) {
            commands += 'select all;';
            commands += 'centerAt average;';
            commands += 'hide none;';
            commands += 'color atoms cpk; color atoms opaque;';
            commands += `dots off; ${wireframe}; ${spacefill};`;
        } else {
            // center of the environment (or structure)
            if (settings.environments.center.value) {
                commands += `select @${this._highlighted + 1};`;
            } else {
                commands += 'select all;';
            }
            commands += 'centerAt average;';

            // Atoms not in the environment
            commands += 'select all;';
            commands += this._backgroundStyle();

            // atoms in the environment
            const cutoff = settings.environments.cutoff.value;
            commands += `select within(${cutoff}, false, @${this._highlighted + 1});`;
            commands += 'hide hidden and not selected;';
            commands += 'color atoms cpk; color atoms opaque;';
            commands += `dots off; ${wireframe}; ${spacefill};`;

            // central atom
            commands += `select @${this._highlighted + 1};`;
            commands += 'hide hidden and not selected;';
            commands += 'color atoms cpk; color atoms opaque;';
            commands += `wireframe off; ${spacefill};`;
            commands += 'color dots green; dots 0.6;';
        }

        commands += `
            select all;
            label ${settings.atomLabels.value ? '%a' : 'off'};
            unitcell ${settings.unitCell.value ? '2' : 'off'};
            spin ${settings.rotation.value ? 'on' : 'off'};
        `;

        return commands + this._showAxesCommands(settings.axes.value);
    }

    private _setCellInfo() {
        const a = this._options.supercell[0].value;
        const b = this._options.supercell[1].value;
        const c = this._options.supercell[2].value;

        this._cellInfo.innerText = this._options.packedCell.value ? 'packed' : 'standard';
        if (a !== 1 || b !== 1 || c !== 1) {
            this._cellInfo.innerText += ` ${a}x${b}x${c} supercell`;
        } else {
            this._cellInfo.innerText += ' cell';
        }
    }

    /// Background *atoms*
    private _backgroundStyle(): string {
        let commands = '';

        const wireframe = `wireframe ${this._options.bonds.value ? '0.15' : 'off'}`;
        const style = this._options.environments.bgStyle.value;
        if (style === 'licorice') {
            commands += `hide none; dots off; ${wireframe}; spacefill off;`;
        } else if (style === 'ball-stick') {
            commands += `hide none; dots off; ${wireframe}; spacefill 23%;`;
        } else if (style === 'hide') {
            commands += 'hide *;';
        } else {
            throw Error(`invalid background atoms style '${style}'`);
        }

        const color = this._options.environments.bgColor.value;
        if (color === 'CPK') {
            commands += 'color atoms cpk;';
            commands += 'color atoms translucent 0.8;';
        } else if (color === 'grey') {
            commands += 'color atoms [xa0a0a0];';
            commands += 'color atoms translucent 0.8;';
        } else {
            throw Error(`invalid background atoms color '${color}'`);
        }

        return commands;
    }

    private _enableEnvironmentSettings(show: boolean) {
        if (this._resetEnvCutof.disabled === show) {
            this._toggleEnvironmentSettings();
        }
    }

    private _toggleEnvironmentSettings() {
        const toggleGroup = getByID(`${this.guid}-env-activated`);
        assert(toggleGroup.parentElement !== null);
        const toggle = toggleGroup.parentElement.lastChild;
        assert(toggle !== null);

        const reset = this._resetEnvCutof;
        if (reset.disabled) {
            reset.disabled = false;
            toggle.nodeValue = 'Disable';

            this._options.environments.cutoff.enable();
            this._options.environments.bgStyle.enable();
            this._options.environments.bgColor.enable();

            // Can not have both environments and packed cell
            this._options.packedCell.value = false;
            this._options.packedCell.disable();
            this._setCellInfo();
        } else {
            reset.disabled = true;
            toggle.nodeValue = 'Enable';

            this._options.environments.cutoff.disable();
            this._options.environments.bgStyle.disable();
            this._options.environments.bgColor.disable();

            this._options.packedCell.enable();
            this._setCellInfo();
        }
    }

    /// Get the current displayed environment, if any. This function deals with
    /// JSmol adding new atoms when loading a supercell, and will always return
    /// a value smaller than the number of atoms in the structure passed to
    /// `JSmolWidget.load`.
    private _environment(): number | undefined {
        if (this._highlighted === undefined || this._environments === undefined) {
            return undefined;
        }

        // `this._natoms ?? 1` is needed because this code might be called
        // before the loading of the structure finishes.
        //
        // In this case, this._highlighted is assumed to be inside [0,
        // this._environments.length). This should be fine since
        // this._highlighted is only set in two places: in this.highlight, and
        // the user should give valid input, or through this.onSelected, and
        // then the structure should be fully loaded
        const natoms = this._natoms ?? 1;

        // % is used to deal with supercell, which create atom indexe
        // outside of [0, this._environments.length)
        return this._highlighted % natoms;
    }

    /// Get the default cutoff for the currently displayed environment
    private _currentDefaultCutoff(): number {
        const i = this._environment();
        if (i === undefined) {
            throw Error('no environments defined when calling _currentCutoff');
        } else {
            assert(this._environments !== undefined);
            return this._environments[i].cutoff;
        }
    }

    /// do we have at least one loaded structure?
    private _loaded(): boolean {
        return this._natoms !== undefined;
    }
}
