/**
 * @packageDocumentation
 * @module structure
 */

import assert from 'assert';

import { default as $3Dmol } from './3dmol';
import { assignBonds } from './3dmol/assignBonds';

import { SavedSettings } from '../options';
import { generateGUID, getByID, getElement, unreachable } from '../utils';
import { PositioningCallback } from '../utils';
import { Structure } from '../dataset';

import { StructureOptions } from './options';

require('../static/chemiscope.css');

/** @hidden
 * Create a stylesheet in the main `document` with the given `rules`
 */
function createStyleSheet(rules: string[]): CSSStyleSheet {
    const style = document.createElement('style');
    document.head.appendChild(style);
    const sheet = style.sheet as CSSStyleSheet;
    for (const rule of rules) {
        sheet.insertRule(rule);
    }
    return sheet;
}

/**
 * Add data from the `structure` to the `model`
 *
 * @param model 3Dmol GLModel that will contain structure data
 * @param structure the structure to convert
 */
function setup3DmolStructure(model: $3Dmol.GLModel, structure: Structure): void {
    if (structure.cell !== undefined) {
        const cell = structure.cell;
        // prettier-ignore
        const matrix = new $3Dmol.Matrix3(
            cell[0], cell[3], cell[6],
            cell[1], cell[4], cell[7],
            cell[2], cell[5], cell[8]
        );
        model.setCrystMatrix(matrix);
    }

    const atoms = [];
    for (let i = 0; i < structure.size; i++) {
        const x = structure.x[i];
        const y = structure.y[i];
        const z = structure.z[i];
        atoms.push({
            serial: i,
            elem: structure.names[i],
            x: x,
            y: y,
            z: z,
        });
    }

    model.addAtoms(atoms);
}

interface LabeledArrow {
    label: $3Dmol.Label;
    arrow: $3Dmol.GLShape;
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

/** Possible options passed to `MoleculeViewer.load` */
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
     * All HTML elements created by this class use this ID to ensure uniqueness.
     */
    public guid: string;

    /// The HTML element serving as root element for the viewer
    private _root: HTMLElement;

    /// Instance of the 3Dmol viewer
    private _viewer: $3Dmol.GLViewer;
    /// Currently displayed structure, if any
    private _current?: {
        /// the structure itself
        structure: Structure;
        /// Corresponding 3Dmol Model
        model: $3Dmol.GLModel;
        /// Atomic labels
        atomLabels: $3Dmol.Label[];
    };
    /// Currently highlighted environment, if any
    private _highlighted?: {
        /// index of the central atom
        center: number;
        /// 3Dmol Model containing only atoms inside the spherical cutoff
        model: $3Dmol.GLModel;
    };
    /// Axes shown, if any
    private _axes?: [LabeledArrow, LabeledArrow, LabeledArrow];

    /// Representation options from the HTML side
    private _options: StructureOptions;
    /// The supercell used to initialize the viewer
    private _initialSupercell?: [number, number, number];
    // button to reset the environment cutoff to its original value
    private _resetEnvCutoff!: HTMLButtonElement;
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
     * Create a new `MoleculeViewer` inside the HTML DOM element with the given `id`.
     *
     * @param element HTML element or HTML id of the DOM element
     *                where the viewer will be created
     * @param guid (optional) unique identifier for the widget
     */
    constructor(element: string | HTMLElement, guid?: string) {
        if (guid === undefined) {
            guid = generateGUID();
        }

        // add a 'chsp-' prefix to ensure the id start with letter. It looks like
        // if the id start with a number (2134950-ffff-4879-82d8-5c9f81dd00ab)
        // then bootstrap code linking modal button to the modal fails ¯\_(ツ)_/¯
        this.guid = 'chsp-' + guid;

        this._root = document.createElement('div');

        const root = getElement(element);
        root.appendChild(this._root);

        this._root.style.position = 'relative';
        this._root.id = this.guid;
        this._root.style.width = '100%';
        this._root.style.height = '100%';

        const viewer = $3Dmol.createViewer(this._root, {
            antialias: true,
            defaultcolors: $3Dmol.elementColors.Jmol,
            backgroundColor: '0xffffff',
            backgroundAlpha: 0,
            disableFog: true,
            orthographic: false,
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
     * Remove all HTML added by this [[MoleculeViewer]] in the current document
     */
    public remove(): void {
        if (this._root.parentElement !== null) {
            this._root.parentElement.innerHTML = '';
        }
        this._options.remove();
    }

    /**
     * Resize the 3Dmol viewer/canvas to the size of the HTML element containing
     * this widget.
     */
    public resize(): void {
        this._viewer.resize();
    }

    /**
     * Get the number of atoms in the structure, or `undefined` if no structure
     * is currently loaded
     *
     * @return the number of atoms in the currently loaded structure
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

        let keepOrientation: boolean;
        if (options.keepOrientation === undefined) {
            // keep pre-existing settings if any
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

        if (this._options.unitCell.value && this._current !== undefined) {
            this._viewer.removeUnitCell(this._current.model);
        }

        // unload previous structure
        this._viewer.removeAllModels();
        if (this._current !== undefined) {
            for (const label of this._current.atomLabels) {
                this._viewer.removeLabel(label);
            }
        }

        // load new structure
        this._current = {
            model: this._viewer.addModel(),
            structure: structure,
            atomLabels: [],
        };
        setup3DmolStructure(this._current.model, structure);
        this._viewer.replicateUnitCell(
            this._options.supercell[0].value,
            this._options.supercell[1].value,
            this._options.supercell[2].value,
            this._current.model
        );

        if (this._options.unitCell.value) {
            this._viewer.addUnitCell(this._current.model, {
                box: { color: 'black' },
                astyle: { hidden: true },
                bstyle: { hidden: true },
                cstyle: { hidden: true },
            });
        }
        assignBonds(this._current.model.selectedAtoms({}) as $3Dmol.AtomSpec[]);
        this._current.model.addAtomSpecs(['index']);

        this._viewer.setClickable({}, true, (atom) => this._selectAtom(atom));

        if (this._environments === undefined) {
            this._styles.noEnvs.disabled = false;
            this._changeHighlighted(undefined);
        } else {
            this._styles.noEnvs.disabled = true;
            this._changeHighlighted(options.highlight === undefined ? 0 : options.highlight);
        }

        this._updateStyle();

        const centerView = this._options.environments.center.value;
        if (!keepOrientation || centerView) {
            this._resetView(centerView);
        }

        // make sure to reset axes/labels when the structure changes
        this._options.axes.onchange(this._options.axes.value, 'JS');
        this._options.atomLabels.onchange(this._options.atomLabels.value, 'JS');

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

    /**
     * Returns a PNG screenshot of the viewer as a URI string
     */
    public exportPNG(): string {
        return this._viewer.pngURI();
    }

    private _connectOptions(): void {
        const restyleAndRender = () => {
            this._updateStyle();
            this._viewer.render();
        };

        this._options.spaceFilling.onchange = restyleAndRender;
        this._options.bonds.onchange = restyleAndRender;

        this._options.atomLabels.onchange = (showLabels) => {
            if (this._current === undefined) {
                return;
            }

            // move the atomic label a bit further away from the atom by
            // monkey-patching 3Dmol
            $3Dmol.SpriteAlignment.bottomLeft = new $3Dmol.Vector2(1.5, 1.5);

            if (showLabels) {
                assert(this._current.atomLabels.length === 0);

                const structure = this._current.structure;
                for (let i = 0; i < structure.size; i++) {
                    const position = new $3Dmol.Vector3(
                        structure.x[i],
                        structure.y[i],
                        structure.z[i]
                    );
                    const name = structure.names[i];

                    let color = $3Dmol.elementColors.Jmol[name] || 0x000000;
                    if (color === 0xffffff || color === 'white') {
                        color = 0x000000;
                    }

                    const label = this._viewer.addLabel(name, {
                        position: position,
                        inFront: true,
                        fontColor: color,
                        fontSize: 14,
                        showBackground: false,
                        alignment: 'bottomLeft',
                    });
                    this._current.atomLabels.push(label);
                }
            } else {
                // remove all labels
                for (const label of this._current.atomLabels) {
                    this._viewer.removeLabel(label);
                }
                this._current.atomLabels = [];
            }
        };

        this._options.axes.onchange = (value) => {
            if (this._axes !== undefined) {
                this._viewer.removeShape(this._axes[0].arrow);
                this._viewer.removeLabel(this._axes[0].label);
                this._viewer.removeShape(this._axes[1].arrow);
                this._viewer.removeLabel(this._axes[1].label);
                this._viewer.removeShape(this._axes[2].arrow);
                this._viewer.removeLabel(this._axes[2].label);
                this._axes = undefined;
            }

            if (value === 'off') {
                // nothing to do
            } else if (value === 'xyz') {
                this._axes = [
                    this._addLabeledArrow([2, 0, 0], 'red', 'X'),
                    this._addLabeledArrow([0, 2, 0], 'green', 'Y'),
                    this._addLabeledArrow([0, 0, 2], 'blue', 'Z'),
                ];
            } else if (value === 'abc') {
                if (this._current === undefined || this._current.structure.cell === undefined) {
                    return;
                }
                const cell = this._current.structure.cell;
                const a: [number, number, number] = [cell[0], cell[1], cell[2]];
                const b: [number, number, number] = [cell[3], cell[4], cell[5]];
                const c: [number, number, number] = [cell[6], cell[7], cell[8]];
                this._axes = [
                    this._addLabeledArrow(a, 'red', 'A'),
                    this._addLabeledArrow(b, 'green', 'B'),
                    this._addLabeledArrow(c, 'blue', 'C'),
                ];
            }

            this._viewer.render();
        };

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
            assignBonds(this._current.model.selectedAtoms({}) as $3Dmol.AtomSpec[]);

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
        this._resetEnvCutoff = getByID<HTMLButtonElement>(`${this.guid}-env-reset`);
        this._resetEnvCutoff.onclick = () => {
            this._options.environments.cutoff.value = this._currentDefaultCutoff();
            restyleAndRender();
        };

        const alignX = getByID<HTMLButtonElement>(`${this.guid}-align-x`);
        alignX.onclick = () => this._viewAlong([1, 0, 0]);

        const alignY = getByID<HTMLButtonElement>(`${this.guid}-align-y`);
        alignY.onclick = () => this._viewAlong([0, 1, 0]);

        const alignZ = getByID<HTMLButtonElement>(`${this.guid}-align-z`);
        alignZ.onclick = () => this._viewAlong([0, 0, 1]);

        const alignA = getByID<HTMLButtonElement>(`${this.guid}-align-a`);
        alignA.onclick = () => {
            if (this._current === undefined || this._current.structure.cell === undefined) {
                return;
            }
            const cell = this._current.structure.cell;
            const a: [number, number, number] = [cell[0], cell[1], cell[2]];
            this._viewAlong(a);
        };

        const alignB = getByID<HTMLButtonElement>(`${this.guid}-align-b`);
        alignB.onclick = () => {
            if (this._current === undefined || this._current.structure.cell === undefined) {
                return;
            }
            const cell = this._current.structure.cell;
            const b: [number, number, number] = [cell[3], cell[4], cell[5]];
            this._viewAlong(b);
        };

        const alignC = getByID<HTMLButtonElement>(`${this.guid}-align-c`);
        alignC.onclick = () => {
            if (this._current === undefined || this._current.structure.cell === undefined) {
                return;
            }
            const cell = this._current.structure.cell;
            const c: [number, number, number] = [cell[6], cell[7], cell[8]];
            this._viewAlong(c);
        };

        this._resetSupercell = getByID<HTMLButtonElement>(`${this.guid}-reset-supercell`);
        this._resetSupercell.onclick = () => {
            assert(this._initialSupercell !== undefined);
            this._options.supercell[0].value = this._initialSupercell[0];
            this._options.supercell[1].value = this._initialSupercell[1];
            this._options.supercell[2].value = this._initialSupercell[2];
        };

        // Reset zoom level when double clicked
        this._root.ondblclick = () => {
            this._resetView(this._options.environments.center.value);
            this._viewer.render();
        };
    }

    /**
     * Function called whenever the user click on an atom in the viewer
     */
    private _selectAtom(atom: Partial<$3Dmol.AtomSpec>): void {
        // use atom.serial instead of atom.index to ensure we are getting the
        // id of the atom inside the central cell when using a supercell
        assert(atom.serial !== undefined);

        if (this._environments !== undefined) {
            this.highlight(atom.serial);
        }

        this.onselect(atom.serial);
    }

    /**
     * Update the styles of all atoms as required
     */
    private _updateStyle(): void {
        if (this._current === undefined) {
            return;
        }

        // if there is no environment to highlight, render all atoms with the
        // main style
        if (!this._environmentsEnabled()) {
            this._current.model.setStyle({}, this._mainStyle());
            this._viewer.render();
            return;
        }

        assert(this._highlighted !== undefined);
        // otherwise, render all atoms with hidden/background style
        this._current.model.setStyle({}, this._hiddenStyle());
        // the central atom with central/highlighted style
        this._current.model.setStyle(
            { index: this._highlighted.center },
            this._centralStyle(),
            /* add: */ true
        );

        // and the environment around the central atom with main style
        this._highlighted.model.setStyle({}, this._mainStyle());
    }

    /**
     * Get the main style used for all atoms/atoms inside the environment when
     * highlighting a specific environment
     */
    private _mainStyle(): Partial<$3Dmol.AtomStyleSpec> {
        const style: Partial<$3Dmol.AtomStyleSpec> = {
            sphere: {
                scale: this._options.spaceFilling.value ? 1.0 : 0.22,
            },
        };
        if (this._options.bonds.value) {
            style.stick = {
                radius: 0.15,
            };
        }

        return style;
    }

    /**
     * Get the style specification for the hidden/background atoms when
     * highlighting a specific environment
     */
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
                    scale: 0.22,
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

    /**
     * Get the style specification for the central atom when
     * highlighting a specific environment
     */
    private _centralStyle(): Partial<$3Dmol.AtomStyleSpec> {
        return {
            sphere: {
                scale: 0.4,
                color: 'green',
                opacity: 0.7,
            },
        };
    }

    /**
     * Show the information related to supercell in a small box on the bottom
     * right corner.
     */
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

    /**
     * Check whether environments are enabled or not
     */
    private _environmentsEnabled(): boolean {
        return this._highlighted !== undefined && !this._resetEnvCutoff.disabled;
    }

    /**
     * Enable (if `show` is true) or disable (if `show` is false) the settings
     * related to environments
     */
    private _enableEnvironmentSettings(show: boolean): void {
        const toggleGroup = getByID(`${this.guid}-env-activated`);
        assert(toggleGroup.parentElement !== null);
        const toggle = toggleGroup.parentElement.lastChild;
        assert(toggle !== null);
        const reset = this._resetEnvCutoff;

        if (show) {
            if (this._environmentsEnabled()) {
                // nothing to do
                return;
            }

            reset.disabled = false;
            toggle.nodeValue = 'Disable';

            this._options.environments.cutoff.enable();
            this._options.environments.bgStyle.enable();
            this._options.environments.bgColor.enable();
        } else {
            if (!this._environmentsEnabled()) {
                // nothing to do
                return;
            }

            reset.disabled = true;
            toggle.nodeValue = 'Enable';

            this._options.environments.cutoff.disable();
            this._options.environments.bgStyle.disable();
            this._options.environments.bgColor.disable();
        }
    }

    /**
     * Get the default cutoff for the currently displayed environment
     */
    private _currentDefaultCutoff(): number {
        if (this._highlighted === undefined) {
            throw Error('no central environments defined when calling _currentCutoff');
        } else {
            assert(this._environments !== undefined);
            return this._environments[this._highlighted.center].cutoff;
        }
    }

    /**
     * Change which central atom is highlighted in the system to `center`. If
     * `center` is undefined, this disable highlighting.
     *
     * @param center index of the atom to highlight
     */
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
            // for the background & highlighted atoms
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

    /**
     * Reset the view by re-centering it and zooming to fit the model as much as
     * possible inside the views.
     *
     * @param center Should we center on the highlighted central atom? If false,
     *               the camera center/look-at is the center of the whole
     *               structure.
     */
    private _resetView(center: boolean): void {
        // HACK: orthographic camera is broken in 3Dmol
        // (https://github.com/3dmol/3Dmol.js/issues/434) so we fake on by using
        // a perspective camera set far away, with a very narrow field of view.
        this._viewer.setCameraParameters({ fov: 0.01, z: 2000 });

        if (center && this._highlighted !== undefined) {
            this._viewer.zoomTo({ serial: this._highlighted.center });
            this._viewer.zoom(4.0);
        } else {
            this._viewer.zoomTo();
        }
        this._viewer.setSlab(-1000, 1000);
    }

    /**
     * Rotate the viewed group so that the given direction (in group
     * coordinates) is aligned with the z axis (in camera space)
     *
     * @param direction axis to align with the camera view
     */
    private _viewAlong(direction: [number, number, number]): void {
        const norm = Math.sqrt(
            direction[0] * direction[0] + direction[1] * direction[1] + direction[2] * direction[2]
        );

        // angle between Oz and the axis
        const angle = Math.acos(direction[2] / norm);

        const quaternion = [
            // rotation axis is direction ^ Oz
            // rotation[0] * sin(angle / 2)
            (direction[1] / norm) * Math.sin(angle / 2),
            // rotation[1] * sin(angle / 2)
            (-direction[0] / norm) * Math.sin(angle / 2),
            // rotation[2] * sin(angle / 2)
            0,
            Math.cos(angle / 2),
        ];

        const viewpoint = this._viewer.getView();
        viewpoint[4] = quaternion[0];
        viewpoint[5] = quaternion[1];
        viewpoint[6] = quaternion[2];
        viewpoint[7] = quaternion[3];
        this._viewer.setView(viewpoint);
    }

    /**
     * Add a labeled arrow of the given color from 0 to the given position, with
     * the given label.
     */
    private _addLabeledArrow(
        position: [number, number, number],
        color: string,
        label: string
    ): LabeledArrow {
        const pos = new $3Dmol.Vector3(position[0], position[1], position[2]);

        return {
            arrow: this._viewer.addArrow({
                start: new $3Dmol.Vector3(0, 0, 0),
                end: pos,
                radius: 0.1,
                color: color,
                midpos: -1,
            }),
            label: this._viewer.addLabel(label, {
                position: pos,
                inFront: true,
                fontColor: 'black',
                fontSize: 14,
                showBackground: false,
            }),
        };
    }
}
