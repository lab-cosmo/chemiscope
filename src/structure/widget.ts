/**
 * @packageDocumentation
 * @module structure
 */

import assert from 'assert';

import { default as $3Dmol } from './3dmol';

import { SavedSettings } from '../options';
import { generateGUID, getByID, sendWarning, unreachable } from '../utils';
import { PositioningCallback } from '../utils';
import { Structure } from '../dataset';

import { StructureOptions } from './options';
import { $3DmolStructure } from './utils';

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
    return sheet;
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

/** Possible options passed to `JSmolWidget.load` */
export interface LoadOptions {
    /** Supercell to display (default: [1, 1, 1]) */
    supercell: [number, number, number];
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

export class MoleculeViewer {
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

    private _viewer: $3Dmol.GLViewer;
    private _current?: {
        structure: Structure;
        model: $3Dmol.GLModel;
    };
    private _highlighted?: {
        center: number;
        model: $3Dmol.GLModel;
    };

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

    /// Dynamic CSS used to hide options as needed
    private _styles: {
        /// hide options related to unit cell if there is no unit cell in the
        /// current structure
        noCell: CSSStyleSheet;
        /// hide options related to environments if there are no environments
        /// for the current structure
        noEnvs: CSSStyleSheet;
    };
    /// List of atom-centered environments for the current structure
    private _environments?: Environment[];

    /**
     * Create a new JSmolWidget inside the HTML DOM element with the given `id`.
     *
     * @param id HTML element id inside which the viewer will be created
     * @param j2sPath path where j2s files can be loaded by Jmol
     * @param guid (optional) unique identifier for the widget
     * @param serverURL URL where to find `jsmol.php`
     */
    constructor(id: string, guid?: string) {
        if (guid === undefined) {
            guid = generateGUID();
        }

        // add a 'chsp-' prefic to ensure the id start with letter. It looks like
        // if the id start with a number (2134950-ffff-4879-82d8-5c9f81dd00ab)
        // then bootstrap code linking modal button to the modal fails ¯\_(ツ)_/¯
        this.guid = 'chsp-' + guid;

        this._root = document.createElement('div');
        const root = getByID(id);
        root.appendChild(this._root);

        this._root.style.position = 'relative';
        this._root.id = this.guid;
        this._root.style.width = '100%';
        this._root.style.height = '100%';

        const viewer = $3Dmol.createViewer(this._root, {
            antialias: true,
            defaultcolors: $3Dmol.elementColors.Jmol,
            disableFog: true,
            orthographic: true,
        });
        if (viewer === undefined) {
            throw Error('unable to create WebGL canvas');
        }
        this._viewer = viewer;
        this.onselect = () => {};

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

        this._connectOptions();

        this._styles = {
            noCell: createStyleSheet([
                `#${this.guid} .chsp-hide-if-no-cell { display: none; }`,
                `#${this.guid}-settings .chsp-hide-if-no-cell { display: none; }`,
            ]),
            noEnvs: createStyleSheet([
                `#${this.guid}-settings .chsp-hide-if-no-environments { display: none; }`,
            ]),
        };

        // By default, position the modal for settings on top of the viewer,
        // centered horizontally
        this.positionSettingsModal = (rect: DOMRect) => {
            const rootRect = this._root.getBoundingClientRect();
            return {
                left: rootRect.left + rootRect.width / 2 - rect.width / 2,
                top: rootRect.top + 20,
            };
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

    public resize(): void {
        this._viewer.resize();
    }

    /**
     * Get the number of atoms in the structure, or `undefined` if no structure
     * is currenly loaded
     *
     * @return the number of atoms in the currenly loaded structure
     */
    public natoms(): number | undefined {
        if (this._current === undefined) {
            return undefined;
        } else {
            return this._current.structure.size;
        }
    }

    /**
     * Load the given `structure` in this viewer.
     *
     * @param structure structure to load
     * @param options options for the new structure
     */
    public load(structure: Structure, options: Partial<LoadOptions> = {}): void {
        // if the canvas size changed since last structure, make sure we update
        // everything
        this.resize();

        // Deal with loading options
        this._environments = options.environments;

        // TODO
        let keepOrientation: boolean;
        if (options.keepOrientation === undefined) {
            // keep pre-existting settings if any
            keepOrientation = this._options.keepOrientation.value;
        } else {
            keepOrientation = options.keepOrientation;
        }

        let a, b, c;
        if (options.supercell === undefined) {
            // keep pre-existing supercell settings, default to [1, 1, 1] from
            // settings.html
            [a, b, c] = [
                this._options.supercell[0].value,
                this._options.supercell[1].value,
                this._options.supercell[2].value,
            ];
        } else {
            [a, b, c] = options.supercell;
            this._options.supercell[0].value = a;
            this._options.supercell[1].value = b;
            this._options.supercell[2].value = c;
        }
        if (this._initialSupercell === undefined) {
            this._initialSupercell = [a, b, c];
            this._resetSupercell.innerHTML = `reset ${a}x${b}x${c} supercell`;
        }

        if (structure.cell === undefined) {
            this._styles.noCell.disabled = false;
        } else {
            this._styles.noCell.disabled = true;
        }
        this._showSupercellInfo();

        const trajectoryOptions = getByID(`${this.guid}-trajectory-settings-group`);
        if (options.trajectory === undefined || !options.trajectory) {
            trajectoryOptions.style.display = 'none';
        } else {
            trajectoryOptions.style.display = 'block';
        }

        // Actual loading
        this._viewer.removeAllModels();
        this._current = {
            model: this._viewer.addModel(),
            structure: structure,
        };
        $3DmolStructure(this._current.model, structure);
        this._current.model.addAtomSpecs(['index']);

        this._viewer.setClickable({}, true, (atom) => this._atomSelected(atom));

        if (this._environments === undefined) {
            this._styles.noEnvs.disabled = false;
            this._changeHighlighted(undefined);
        } else {
            this._styles.noEnvs.disabled = true;
            this._changeHighlighted(options.highlight === undefined ? 0 : options.highlight);
        }

        this._updateStyle();

        if (!keepOrientation) {
            this._resetView(this._options.environments.center.value);
        }

        this._viewer.render();
    }

    /**
     * Highlight a given `atom` in the current structure.
     *
     * If a supercell larger than [1, 1, 1] is currently displayed, this
     * function accept indexes larger than the result of `natoms()`, and will
     * then highlight atoms outside of the central image.
     *
     * @param center index of the central atom in the environment to show,
     *               or `undefined` to disable highlighting.
     */
    public highlight(center?: number): void {
        this._changeHighlighted(center);
        this._updateStyle();

        const centerView = this._options.environments.center.value;
        if (this._highlighted !== undefined && centerView) {
            this._resetView(centerView);
        }

        this._viewer.render();
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

    private _connectOptions(): void {
        const restyleAndRender = () => {
            this._updateStyle();
            this._viewer.render();
        };

        this._options.atomLabels.onchange = restyleAndRender;
        this._options.spaceFilling.onchange = restyleAndRender;
        this._options.bonds.onchange = restyleAndRender;

        this._options.axes.onchange = () => sendWarning('TODO: axes settings');

        this._options.rotation.onchange = (rotate) => {
            if (rotate) {
                this._viewer.spin('vy');
            } else {
                this._viewer.spin(false);
            }
        };

        this._options.unitCell.onchange = (add) => {
            if (this._current === undefined) {
                return;
            }

            if (add) {
                this._viewer.addUnitCell(this._current.model, {
                    box: { color: 'black' },
                    astyle: { hidden: true },
                    bstyle: { hidden: true },
                    cstyle: { hidden: true },
                });
            } else {
                this._viewer.removeUnitCell(this._current.model);
            }
            this._viewer.render();
        };

        const changedSuperCell = () => {
            this._showSupercellInfo();
            if (this._current === undefined) {
                return;
            }

            // remove the atoms previously added by `replicateUnitCell`
            const atoms = this._current.model.selectedAtoms({});
            this._current.model.removeAtoms(atoms.splice(this._current.structure.size));

            this._viewer.replicateUnitCell(
                this._options.supercell[0].value,
                this._options.supercell[1].value,
                this._options.supercell[2].value,
                this._current.model
            );

            if (this._highlighted !== undefined) {
                this._changeHighlighted(this._highlighted.center);
            }
            this._updateStyle();
            this._viewer.render();
        };
        this._options.supercell[0].onchange = changedSuperCell;
        this._options.supercell[1].onchange = changedSuperCell;
        this._options.supercell[2].onchange = changedSuperCell;

        this._options.environments.bgColor.onchange = restyleAndRender;
        this._options.environments.bgStyle.onchange = restyleAndRender;
        this._options.environments.cutoff.onchange = () => {
            if (this._highlighted !== undefined) {
                this._changeHighlighted(this._highlighted.center);
            }
            restyleAndRender();
        };
        this._options.environments.center.onchange = (center) => {
            this._resetView(center);
            this._viewer.render();
        };

        // Deal with activation/de-activation of environments
        this._options.environments.activated.onchange = (value) => {
            this._enableEnvironmentSettings(value);
            restyleAndRender();
        };

        // Setup various buttons
        this._resetEnvCutof = getByID<HTMLButtonElement>(`${this.guid}-env-reset`);
        this._resetEnvCutof.onclick = () => {
            this._options.environments.cutoff.value = this._currentDefaultCutoff();
            restyleAndRender();
        };

        const alignX = getByID<HTMLButtonElement>(`${this.guid}-align-x`);
        alignX.onclick = () => sendWarning('TODO: camera change');

        const alignY = getByID<HTMLButtonElement>(`${this.guid}-align-y`);
        alignY.onclick = () => sendWarning('TODO: camera change');

        const alignZ = getByID<HTMLButtonElement>(`${this.guid}-align-z`);
        alignZ.onclick = () => sendWarning('TODO: camera change');

        const alignA = getByID<HTMLButtonElement>(`${this.guid}-align-a`);
        alignA.onclick = () => sendWarning('TODO: camera change');

        const alignB = getByID<HTMLButtonElement>(`${this.guid}-align-b`);
        alignB.onclick = () => sendWarning('TODO: camera change');

        const alignC = getByID<HTMLButtonElement>(`${this.guid}-align-c`);
        alignC.onclick = () => sendWarning('TODO: camera change');

        this._resetSupercell = getByID<HTMLButtonElement>(`${this.guid}-reset-supercell`);
        this._resetSupercell.onclick = () => {
            assert(this._initialSupercell !== undefined);
            this._options.supercell[0].value = this._initialSupercell[0];
            this._options.supercell[1].value = this._initialSupercell[1];
            this._options.supercell[2].value = this._initialSupercell[2];
        };
    }

    private _atomSelected(atom: Partial<$3Dmol.AtomSpec>): void {
        // use atom.serial instead of atom.index to ensure we are getting the
        // id of the atom inside the central cell when using a supercell
        assert(atom.serial !== undefined);

        if (this._environments !== undefined) {
            this.highlight(atom.serial);
        }

        this.onselect(atom.serial);
    }

    private _updateStyle(): void {
        if (this._current === undefined) {
            return;
        }

        if (this._highlighted === undefined) {
            this._current.model.setStyle({}, this._mainStyle());
            this._viewer.render();
            return;
        }

        this._current.model.setStyle({}, this._hiddenStyle());
        this._current.model.setStyle(
            { index: this._highlighted.center },
            this._centralStyle(),
            /* add: */ true
        );

        this._highlighted.model.setStyle({}, this._mainStyle());
    }

    private _mainStyle(): Partial<$3Dmol.AtomStyleSpec> {
        const style: Partial<$3Dmol.AtomStyleSpec> = {
            sphere: {
                scale: this._options.spaceFilling.value ? 1.0 : 0.3,
            },
        };
        if (this._options.bonds.value) {
            style.stick = {
                radius: 0.15,
            };
        }

        return style;
    }

    private _hiddenStyle(): Partial<$3Dmol.AtomStyleSpec> {
        const style: Partial<$3Dmol.AtomStyleSpec> = {};

        const bgStyle = this._options.environments.bgStyle.value;
        if (bgStyle === 'hide') {
            // nothing to do
        } else if (bgStyle === 'licorice' || bgStyle === 'ball-stick') {
            style.stick = {
                radius: 0.15,
                opacity: 0.7,
                hidden: !this._options.bonds.value,
            };

            if (bgStyle === 'ball-stick') {
                style.sphere = {
                    scale: 0.3,
                    opacity: 0.7,
                };
            }
        } else {
            unreachable();
        }

        const bgColor = this._options.environments.bgColor.value;
        if (bgColor === 'CPK') {
            // nothing to do
        } else if (bgColor === 'grey') {
            if (style.stick !== undefined) {
                style.stick.color = 'grey';
            }

            if (style.sphere !== undefined) {
                style.sphere.color = 'grey';
            }
        } else {
            unreachable();
        }

        return style;
    }

    private _centralStyle(): Partial<$3Dmol.AtomStyleSpec> {
        return {
            sphere: {
                scale: 0.4,
                color: 'green',
                opacity: 0.7,
            },
        };
    }

    private _showSupercellInfo(): void {
        const a = this._options.supercell[0].value;
        const b = this._options.supercell[1].value;
        const c = this._options.supercell[2].value;

        if (a !== 1 || b !== 1 || c !== 1) {
            this._cellInfo.innerText = `${a}x${b}x${c} supercell`;
        } else {
            this._cellInfo.innerText = '';
        }
    }

    private _enableEnvironmentSettings(show: boolean): void {
        if (this._resetEnvCutof.disabled === show) {
            this._toggleEnvironmentSettings();
        }
    }

    private _toggleEnvironmentSettings(): void {
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
        } else {
            reset.disabled = true;
            toggle.nodeValue = 'Enable';

            this._options.environments.cutoff.disable();
            this._options.environments.bgStyle.disable();
            this._options.environments.bgColor.disable();
        }
    }

    /// Get the default cutoff for the currently displayed environment
    private _currentDefaultCutoff(): number {
        if (this._highlighted === undefined) {
            throw Error('no central environments defined when calling _currentCutoff');
        } else {
            assert(this._environments !== undefined);
            return this._environments[this._highlighted.center].cutoff;
        }
    }

    private _changeHighlighted(center?: number): void {
        if (this._highlighted !== undefined) {
            this._viewer.removeModel(this._highlighted.model);
        }

        if (center === undefined) {
            this._enableEnvironmentSettings(false);
            this._options.environments.cutoff.value = 0;
            this._highlighted = undefined;
        } else {
            this._enableEnvironmentSettings(true);
            // keep user defined cutoff, if any
            if (this._options.environments.cutoff.value <= 0) {
                this._options.environments.cutoff.value = this._currentDefaultCutoff();
            }

            // We need to create a separate model to have different opacity
            // in the background & highlighted atoms
            // https://github.com/3dmol/3Dmol.js/issues/166
            // prettier-ignore
            const selection = { or: [
                    { index: center },
                    { within: { distance: this._options.environments.cutoff.value, sel: { index: center } } },
                ],
            };
            this._highlighted = {
                model: this._viewer.createModelFrom(selection),
                center: center,
            };
        }
    }

    private _resetView(center: boolean): void {
        if (center && this._highlighted !== undefined) {
            this._viewer.zoomTo({ serial: this._highlighted.center });
        } else {
            this._viewer.zoomTo();
        }

        this._viewer.zoom(2);
        this._viewer.setSlab(-1000, 1000);
    }
}
