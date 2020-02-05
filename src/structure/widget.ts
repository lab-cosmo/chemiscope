/**
 * @packageDocumentation
 * @module structure
 */

import assert from 'assert';

import {JmolObject, JSmolApplet} from 'jsmol';

import {makeDraggable, getByID, generateGUID} from '../utils';

import HTML_SETTINGS from './settings.html';

require('../static/chemiscope.css');


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
    return sheet
}

/**
 * @hidden
 * Get the "555" style supercell description. This allows to put the extra
 * cells both on negative and positive sides.
 */
function supercell_555(supercell: [number, number, number]): string {
    const supercellCenter = [
        Math.floor(supercell[0] / 2),
        Math.floor(supercell[1] / 2),
        Math.floor(supercell[2] / 2),
    ];

    const first = [
        5 - supercellCenter[0],
        5 - supercellCenter[1],
        5 - supercellCenter[2],
    ];

    const second = [
        5 + supercell[0] - (supercellCenter[0] + 1),
        5 + supercell[1] - (supercellCenter[1] + 1),
        5 + supercell[2] - (supercellCenter[2] + 1),
    ];

    return `{${first[0]}${first[1]}${first[2]} ${second[0]}${second[1]}${second[2]} 1}`
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
export interface LoadOption {
    /** Supercell to display (default: [1, 1, 1]) */
    supercell: [number, number, number];
    /** Should we display a packed cell (default: false) */
    packed: boolean;
    /** Should we the current camera orientation (default: false) */
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

/**
 * A light wrapper around JSmol, displaying a single structure. This is kept as
 * a separate class from the [[StructureViewer]] to be easier to use outside of
 * chemiscope.
 */
export class JSmolWidget {
    /// The HTML element serving as root element for the viewer
    private _root: HTMLElement;
    /// Reference to the global Jmol variable
    private _Jmol: JmolObject;
    /// Reference to the JSmol applet to be updated with script calls
    private _applet!: JSmolApplet;
    /// Representation options from the HTML side
    private _options!: {
        // show bonds if checked
        bonds: HTMLInputElement,
        // use space filling representation if checked
        spaceFilling: HTMLInputElement;
        // show atoms labels if checked
        atomLabels: HTMLInputElement;
        // show unit cell information and lines if checked
        unitCell: HTMLInputElement;
        // spin the represenation if checked
        rotation: HTMLInputElement;
        // which axis system to use (none, xyz, abc)
        axes: HTMLSelectElement;
        // options related to environments
        environments: {
            // the cutoff value for spherical environments
            cutoff: HTMLInputElement;
            // which style for atoms not in the environment
            bgStyle: HTMLSelectElement;
            // which colors for atoms not in the environment
            bgColor: HTMLSelectElement;
            // reset the cutoff to its original value
            reset: HTMLButtonElement;
            // enable/disable environments display
            toggle: HTMLButtonElement;
        }
        // keep the orientation constant when loading a new structure if checked
        keepOrientation: HTMLInputElement;
    };
    /// Supercell related options and data
    private _supercell!: {
        /// number of repetitions in the `a` direction for the supercell
        repeat_a: HTMLInputElement;
        /// number of repetitions in the `b` direction for the supercell
        repeat_b: HTMLInputElement;
        /// number of repetitions in the `c` direction for the supercell
        repeat_c: HTMLInputElement;
        /// button used to reset the supercell
        reset: HTMLButtonElement;
        /// Is the current unit cell displayed as a packed cell?
        packed: HTMLInputElement;
        /// The supercell used to intialize the viewer
        initial?: [number, number, number];
    }
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
    /// callback to get the initial positioning of the settings modal. The
    /// callback gets the current placement of the settings as a DOMRect, and
    /// should return top and left positions in pixels, used with
    /// `position: fixed`
    private _settingsPlacement: (rect: DOMRect) => {top: number, left: number};
    /// Number of atoms in the currenly loaded structure 1x1x1 supercell.
    /// If another supercell is displayed, there are na x nb x nc x natoms atoms
    /// in the displayed structure
    private _natoms?: number;
    /// List of atom-centered environments for the current structure
    private _environments?: Environment[];
    /// Index of the central atom to highlight, if any
    private _highlighted?: number;
    /// last selected atom
    private _lastSelected: number;

    /** callback called when a new atom is clicked on */
    public onselect: (atom: number) => void;

    /**
     * Unique identifier of this viewer.
     *
     * All HTML elements created by this class use this ID to ensure unicity.
     */
    public guid: string;

    /**
     * Create a new JSmolWidget inside the HTML DOM element with the given `id`.
     *
     * @param id HTML element id inside which the viewer will be created
     * @param j2sPath path where j2s files can be loaded by Jmol
     * @param serverURL URL where to find `jsmol.php`
     */
    constructor(id: string, j2sPath: string, serverURL: string = "https://chemapps.stolaf.edu/jmol/jsmol/php/jsmol.php") {
        if (window.Jmol === undefined) {
            throw Error('Jmol is required, load it from your favorite source');
        }

        // add a 'chsp-' prefic to ensure the id start with letter. It looks like
        // if the id start with a number (2134950-ffff-4879-82d8-5c9f81dd00ab)
        // then bootstrap code linking modal button to the modal fails ¯\_(ツ)_/¯
        this.guid = 'chsp-' + generateGUID();
        this._lastSelected = -1;

        // store a reference to the global Jmol and the HTML root
        this._Jmol = window.Jmol;
        // Do not insert new applets automatically
        this._Jmol.setDocument(false);

        const root = document.getElementById(id);
        if (root === null) {
            throw Error(`Could not get the "#${id}" element`);
        } else {
            this._root = document.createElement("div");
            root.appendChild(this._root);

            this._root.style.position = "relative";
            this._root.id = this.guid;
            this._root.style.width = "100%";
            this._root.style.height = "100%";
        }

        this._cellInfo = document.createElement("span");
        this._cellInfo.classList.add('chsp-cell-info', 'chsp-hide-if-no-cell', 'badge', 'badge-light');
        this._root.appendChild(this._cellInfo);

        this._createOptions();
        this._createApplet(j2sPath, serverURL);

        this._loadedCallback = undefined;
        // create a global function with unique name and install it as callback
        // for JSmol. This function will then call the callback inside this
        // instance of the viewer
        const name: string = this.guid.replace(/-/g, '_') + '_loaded_callback'
        const window_as_map = window as { [key: string]: any };
        window_as_map[name] = (_a: any, _b: any, _c: any, _d: any, _e: any, status: any) => {
            if (status.valueOf() === 3 && this._loadedCallback !== undefined) {
                this._loadedCallback();
            }
        }
        this.script(`set LoadStructCallback "${name}"`)

        this._noCellStyle = createStyleSheet([
            `#${this.guid} .chsp-hide-if-no-cell { display: none; }`,
            `#${this.guid}-settings .chsp-hide-if-no-cell { display: none; }`
        ])

        this._noEnvsStyle = createStyleSheet([
            `#${this.guid}-settings .chsp-hide-if-no-environments { display: none; }`
        ])

        // By default, position the modal for settings on top of the viewer,
        // centered horizontally
        this._settingsPlacement = (rect: DOMRect) => {
            const rootRect = this._root.getBoundingClientRect();
            return {
                top: rootRect.top + 20,
                left: rootRect.left + rootRect.width / 2 - rect.width / 2,
            }
        }

        this.onselect = () => {};
        this._root.onclick = () => {
            // picked is 0-based
            let picked = this.evaluate('_atomPicked');

            if (picked === "") {
                return;
            } else {
                assert(typeof picked === "number");
            }

            if (picked === this._lastSelected) {
                return;
            } else {
                this._lastSelected = picked;
            }

            this.onselect(picked);
        }
    }

    /**
     * Get the number of atoms in the structure, or `undefined` if no structure
     * is currenly loaded
     *
     * @return the number of atoms in the currenly loaded structure
     */
    public natoms(): number | undefined {
        return this._natoms
    }

    /**
     * Run the given `commands` for this viewer. `load` commands should use the
     * [[JSmolWidget.load|corresponding function]].
     */
    public script(commands: string) {
        if (commands.includes('load ')) {
            throw Error('Do not use JSmolWidget.script to load a structure, but JSmolWidget.load');
        }
        this._Jmol.script(this._applet, commands);
    }

    /**
     * Evaluate the given commands using JSmol and return the corresponding
     * value. This calls `Jmol.evaluateVar` behind the scenes.
     */
    public evaluate(commands: string): any {
        return this._Jmol.evaluateVar(this._applet, commands);
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
    public settingsPlacement(callback: (rect: DOMRect) => Position) {
        this._settingsPlacement = callback;
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
    public highlight(environment?: number) {
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
    public load(data: string, options: Partial<LoadOption> = {}) {
        if (data === "''" || data === '""') {
            throw Error('invalid use of JSmolWidget.load to reload data');
        }

        this._natoms = undefined;

        let supercell: [number, number, number];
        if (options.supercell === undefined) {
            // keep pre-existing supercell settings, default to [1, 1, 1] from
            // settings.html
            supercell = [
                parseInt(this._supercell.repeat_a.value),
                parseInt(this._supercell.repeat_b.value),
                parseInt(this._supercell.repeat_c.value),
            ];
        } else {
            supercell = options.supercell;
        }

        if (options.packed !== undefined) {
            this._supercell.packed.checked = options.packed;
        }

        if (options === undefined) {
            this._environments = undefined;
        } else {
            this._environments = options.environments;
        }

        let keepOrientation: boolean;
        if (options.keepOrientation === undefined) {
            // keep pre-existting settings if any
            keepOrientation = this._options.keepOrientation.checked;
        } else {
            keepOrientation = options.keepOrientation
        }

        const trajectoryOptions = getByID<HTMLElement>(`${this.guid}-trajectory-settings-group`);
        if (options.trajectory === undefined || !options.trajectory) {
            trajectoryOptions.style.display = "none";
        } else {
            trajectoryOptions.style.display = "block";
        }

        if (this._environments === undefined) {
            this._noEnvsStyle.disabled = false;
            this._changeHighlighted(undefined);
        } else {
            if (options.packed !== undefined && options.packed) {
                throw Error("Can not have both packed cell and environments")
            }

            this._noEnvsStyle.disabled = true;
            if (options.highlight !== undefined) {
                this._changeHighlighted(options.highlight);
            } else {
                this._changeHighlighted(0);
            }
        }

        const [a, b, c] = supercell;
        if (this._supercell.initial === undefined) {
            this._supercell.initial = supercell;
            this._supercell.reset.innerHTML = `reset ${a}x${b}x${c} supercell`;
        }
        this._supercell.repeat_a.value = a.toString();
        this._supercell.repeat_b.value = b.toString();
        this._supercell.repeat_c.value = c.toString();
        this._setCellInfo();

        /// JSmol overwite the Error object, using `.caller` to get stacktraces.
        /// This fails in modern browsers if any code other code try to throw an
        /// error, since the JSmol Error object will try to get the stack here
        /// too. So we save and reset the global Error object when loading a new
        /// structure.
        const savedError = window.Error;

        this._loadedCallback = () => {
            window.Error = savedError;

            if (this.evaluate("@{unitcell(\"conventional\")}") === "ERROR") {
                this._noCellStyle.disabled = false;
            } else {
                this._noCellStyle.disabled = true;
            }

            const supercell = parseInt(this._supercell.repeat_a.value) *
                              parseInt(this._supercell.repeat_b.value) *
                              parseInt(this._supercell.repeat_c.value);
            this._natoms = parseInt(this.evaluate("{*}.size")) / supercell;
            if (this._environments !== undefined && this._environments.length !== this._natoms) {
                let message = "invalid number of environments for this structure ";
                message += `got ${this._environments.length} environments, there are ${this._natoms} atoms`
                // We can not throw an error here, since it seems to be caught
                // by JSmol.
                console.error(message);
            }

            // once we are done, disable the callback to prevent it from firing
            // on restyle/reload
            this._loadedCallback = undefined;
        }

        this._do_load(data, keepOrientation);
    }

    private _reload() {
        this._do_load("''", true);
    }

    private _do_load(data: string, keepOrientation: boolean) {
        if (data.includes(';')) {
            throw Error('invalid \';\' in  JSmolWidget.load');
        }
        if (data.includes('packed')) {
            throw Error('invalid \'packed\' in  JSmolWidget.load');
        }

        const packed = this._supercell.packed.checked ? ' packed' : '';

        const saveOrientation = keepOrientation ? `save orientation "${this.guid}"`: '';
        const restoreOrientation = keepOrientation ? `restore orientation "${this.guid}"`: '';

        this._setCellInfo();

        const supercell = supercell_555([
            parseInt(this._supercell.repeat_a.value),
            parseInt(this._supercell.repeat_b.value),
            parseInt(this._supercell.repeat_c.value)
        ]);
        this._Jmol.script(this._applet, `
            ${saveOrientation};
            load ${data} ${supercell} ${packed};
            ${this._updateStateCommands()};
            ${restoreOrientation};
        `);
    }

    private _changeHighlighted(environment?: number) {
        if (environment !== undefined) {
            const supercell = parseInt(this._supercell.repeat_a.value) *
                              parseInt(this._supercell.repeat_b.value) *
                              parseInt(this._supercell.repeat_c.value);
            if (this._natoms !== undefined && environment >= this._natoms * supercell) {
                let message = 'selected environment is out of bounds: ';
                message += `got ${environment}, we have ${this._natoms * supercell} atoms in the current structure`;
                throw Error(message);
            }

            if (this._environments === undefined) {
                throw Error("no environments defined for the current structure");
            }
        }

        this._highlighted = environment;

        if (this._highlighted === undefined) {
            this._enableEnvironmentSettings(false);
            this._options.environments.cutoff.value = "";
        } else {
            this._enableEnvironmentSettings(true);

            // keep user defined cutoff, if any
            if (this._options.environments.cutoff.value === "") {
                this._options.environments.cutoff.value = this._currentDefaultCutoff().toString();
            }
        }
    }

    private _createApplet(j2sPath: string, serverURL: string) {
        const width = this._root.clientWidth;
        const height = this._root.clientHeight;
        if (width === 0 || height === 0) {
            const parentId = this._root.parentElement!.id;
            console.error(
                `JSmolWidget: width (=${width}px) or heigh (=${height}px) ` +
                `of #${parentId} is zero, you will not see the molecules`
            );
        }

        // create the main Jmol applet
        this._applet = this._Jmol.getApplet(this.guid, {
            use: "HTML5",
            script: `
                // use anti-aliasing
                set antialiasdisplay;
                // remove jmol logo
                set frank off;
                // use the smaller of height/width when setting zoom level
                set zoomlarge false;
                // Allow sending script commands while moveto is executing
                set waitformoveto off;
                hide off;
            `,
            disableJ2SLoadMonitor: true,
            disableInitialConsole: true,
            zIndexBase: 1,
            j2sPath: j2sPath,
            serverURL: serverURL,
            width: "100%",
            height: "100%",
        });

        const div = document.createElement("div");
        div.style.height = "100%";
        div.style.width = "100%";
        this._root.appendChild(div);
        div.innerHTML = this._Jmol.getAppletHtml(this._applet);
        // Jmol rely on this script being implicitly executed, but this is not
        // the case when using innerHTML (compared to jquery .html()). So let's
        // manually execute it
        this._applet._cover(false);
    }

    private _createOptions() {
        // use HTML5 template to generate a DOM object from an HTML string
        const template = document.createElement('template');
        template.innerHTML = `<button
            class="btn btn-light btn-sm chsp-open-viewer-settings"
            data-target="#${this.guid}-settings"
            data-toggle="modal">
                <div class="chsp-hamburger"><div></div><div></div><div></div></div>
            </button>`;
        const openSettings = template.content.firstChild!;
        this._root.append(openSettings);

        // replace id to ensure they are unique even if we have mulitple viewers
        // on a single page
        template.innerHTML = HTML_SETTINGS
            .replace(/id=(.*?) /g, (_: string, id: string) => `id=${this.guid}-${id} `)
            .replace(/for=(.*?) /g, (_: string, id: string) => `for=${this.guid}-${id} `);

        const modal = template.content.firstChild!;
        document.body.appendChild(modal);

        const modalDialog = modal.childNodes[1]! as HTMLElement
        if (!modalDialog.classList.contains('modal-dialog')) {
            throw Error("internal error: missing modal-dialog class")
        }

        // Position modal near the actual viewer
        openSettings.addEventListener('click', () => {
            // only set style once, on first open, and keep previous position
            // on next open to keep the 'draged-to' position
            if (modalDialog.getAttribute('data-initial-modal-positions-set') === null) {
                modalDialog.setAttribute('data-initial-modal-positions-set', "true");

                // display: block to ensure modalDialog.offsetWidth is non-zero
                (modalDialog.parentNode as HTMLElement).style.display = 'block';

                const {top, left} = this._settingsPlacement(modalDialog.getBoundingClientRect());

                // set width first, since setting position can influence it
                modalDialog.style.width = `${modalDialog.offsetWidth}px`;
                // unset margins when using position: fixed
                modalDialog.style.margin = '0';
                modalDialog.style.position = 'fixed';
                modalDialog.style.top = `${top}px`;
                modalDialog.style.left = `${left}px`;
            }
        })

        // make the settings modal draggable
        makeDraggable(modalDialog, ".modal-header");

        // get the HTML options and keep them around for use in _updateState
        this._options = {
            bonds:        getByID<HTMLInputElement>(`${this.guid}-bonds`),
            spaceFilling: getByID<HTMLInputElement>(`${this.guid}-space-filling`),
            atomLabels:   getByID<HTMLInputElement>(`${this.guid}-atom-labels`),
            unitCell:     getByID<HTMLInputElement>(`${this.guid}-unit-cell`),
            rotation:     getByID<HTMLInputElement>(`${this.guid}-rotation`),
            axes:         getByID<HTMLSelectElement>(`${this.guid}-axes`),
            environments: {
                cutoff:  getByID<HTMLInputElement>(`${this.guid}-env-cutoff`),
                bgColor: getByID<HTMLSelectElement>(`${this.guid}-env-bg-color`),
                bgStyle: getByID<HTMLSelectElement>(`${this.guid}-env-bg-style`),
                reset:   getByID<HTMLButtonElement>(`${this.guid}-env-reset`),
                toggle:  getByID<HTMLButtonElement>(`${this.guid}-env-toggle`),
            },
            keepOrientation: getByID<HTMLInputElement>(`${this.guid}-keep-orientation`),
        };

        for (const key in this._options) {
            if (key !== 'environments') {
                Reflect.get(this._options, key).onchange = () => this._updateState();
            }
        }

        for (const key in this._options.environments) {
            Reflect.get(this._options.environments, key).onchange = () => this._updateState();
        }

        this._options.environments.reset.onclick = () => {
            this._options.environments.cutoff.value = this._currentDefaultCutoff().toString();
            this._updateState();
        }

        this._options.environments.toggle.onclick = () => {
            this._toggleEnvironmentSettings();
            this._updateState();
        }

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

        // cell handling: we need to reload to change the supercell or the
        // packed status
        this._supercell = {
            repeat_a: getByID<HTMLInputElement>(`${this.guid}-supercell-a`),
            repeat_b: getByID<HTMLInputElement>(`${this.guid}-supercell-b`),
            repeat_c: getByID<HTMLInputElement>(`${this.guid}-supercell-c`),
            packed: getByID<HTMLInputElement>(`${this.guid}-packed-cell`),
            reset: getByID<HTMLButtonElement>(`${this.guid}-reset-supercell`),
        }

        for (const key in this._supercell) {
            if (key !== 'reset') {
                Reflect.get(this._supercell, key).onchange = () => this._reload();
            }
        }

        this._supercell.reset.onclick = () => {
            if (this._supercell.initial === undefined) {
                throw Error("internal bug: this._supercell.initial is undefined")
            }

            const [a, b, c] = this._supercell.initial;
            this._supercell.repeat_a.value = a.toString();
            this._supercell.repeat_b.value = b.toString();
            this._supercell.repeat_c.value = c.toString();
            this._reload();
        }
    }

    private _updateState() {
        this.script(this._updateStateCommands());
    }

    private _showAxesCommands(axes: string): string {
        switch (axes) {
        case "xyz":
            return `
                axes off;
                draw xaxis '>X' vector {0 0 0} {2 0 0} color red width 0.15;
                draw yaxis '>Y' vector {0 0 0} {0 2 0} color green width 0.15;
                draw zaxis '>Z' vector {0 0 0} {0 0 2} color blue width 0.15;
            `;
        case "abc":
            return `
                draw xaxis delete; draw yaxis delete; draw zaxis delete;
                set axesUnitcell ON; axes 5;
            `;
        case "off":
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
        const options = this._options;
        const wireframe = `wireframe ${options.bonds.checked ? '0.15' : 'off'}`;
        const spacefill = `spacefill ${options.spaceFilling.checked ? '80%' : '23%'}`;

        let commands = '';
        if (this._highlighted === undefined || this._options.environments.cutoff.disabled) {
            commands += 'select all;';
            commands += 'hide none;';
            commands += 'color atoms cpk; color atoms opaque;';
            commands += `dots off; ${wireframe}; ${spacefill};`;
        } else {
            // Atoms not in the environment
            commands += 'select all;'
            commands += this._backgroundStyle();

            // atoms in the environment
            const cutoff = parseFloat(this._options.environments.cutoff.value);
            commands += `select within(${cutoff}, false, @${this._highlighted + 1});`;
            commands += 'hide hidden and not selected;';
            commands += 'color atoms cpk; color atoms opaque;';
            commands += `dots off; ${wireframe}; ${spacefill};`;

            // central atom
            commands += `select @${this._highlighted + 1};`
            commands += 'hide hidden and not selected;';
            commands += 'color atoms cpk; color atoms opaque;';
            commands += `wireframe off; ${spacefill};`
            commands += 'color dots green; dots 0.6;';
        }

        commands += `
            select all;
            label ${options.atomLabels.checked ? '%a' : 'off'};
            unitcell ${options.unitCell.checked ? '2' : 'off'};
            spin ${options.rotation.checked ? 'on' : 'off'};
        `;

        return commands + this._showAxesCommands(options.axes.value)
    }

    private _setCellInfo() {
        const a = this._supercell.repeat_a.value;
        const b = this._supercell.repeat_b.value;
        const c = this._supercell.repeat_c.value;

        this._cellInfo.innerText = this._supercell.packed.checked ? 'packed' : 'standard';
        if (a !== '1' || b !== '1' || c !== '1') {
            this._cellInfo.innerText += ` ${a}x${b}x${c} supercell`;
        } else {
            this._cellInfo.innerText += ' cell';
        }
    }

    private _backgroundStyle(): string {
        let commands = '';

        const wireframe = `wireframe ${this._options.bonds.checked ? '0.15' : 'off'}`;
        const style = this._options.environments.bgStyle.value;
        if (style === "licorice") {
            commands += `hide none; dots off; ${wireframe}; spacefill off;`;
        } else if (style === "ball-stick") {
            commands += `hide none; dots off; ${wireframe}; spacefill 23%;`;
        } else if (style === "hide") {
            commands += 'hide *;';
        } else {
            throw Error(`invalid background atoms style '${style}'`)
        }

        const color = this._options.environments.bgColor.value;
        if (color === "CPK") {
            commands += 'color atoms cpk;';
            commands += 'color atoms translucent 0.8;';
        } else if (color === "grey") {
            commands += 'color atoms [xa0a0a0];';
            commands += 'color atoms translucent 0.8;';
        } else {
            throw Error(`invalid background atoms color '${color}'`)
        }

        return commands;
    }

    private _enableEnvironmentSettings(show: boolean) {
        if (this._options.environments.reset.disabled === show) {
            this._toggleEnvironmentSettings();
        }
    }

    private _toggleEnvironmentSettings() {
        const reset = this._options.environments.reset;
        const toggle = this._options.environments.toggle;
        if (reset.disabled) {
            reset.disabled = false;
            reset.classList.toggle('btn-primary');
            reset.classList.toggle('btn-secondary');

            toggle.classList.toggle('btn-primary');
            toggle.classList.toggle('btn-secondary');
            toggle.innerText = 'Disable';

            this._options.environments.cutoff.disabled = false;
            this._options.environments.bgStyle.disabled = false;
            this._options.environments.bgColor.disabled = false;

            // Can not have both environments and packed cell
            this._supercell.packed.checked = false;
            this._supercell.packed.disabled = true;

            this._setCellInfo();
        } else {
            reset.disabled = true;
            reset.classList.toggle('btn-primary');
            reset.classList.toggle('btn-secondary');

            toggle.classList.toggle('btn-primary');
            toggle.classList.toggle('btn-secondary');
            toggle.innerText = 'Enable';

            this._options.environments.cutoff.disabled = true;
            this._options.environments.bgStyle.disabled = true;
            this._options.environments.bgColor.disabled = true;

            this._supercell.packed.disabled = false;
            this._setCellInfo();
        }
    }

    /// Get the current displayed environment, if any. This function deals with
    /// JSmol adding new atoms when loading a supercell, and will always return
    /// a value smaller than the number of atoms in the structure passed to
    /// `JSmolWidget.load`.
    private _environment(): number | undefined {
        if (this._highlighted === undefined ||  this._environments === undefined) {
            return undefined
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
            throw Error("no environments defined when calling _currentCutoff")
        } else {
            return this._environments![i].cutoff;
        }
    }
}
