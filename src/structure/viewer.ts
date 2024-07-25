/**
 * @packageDocumentation
 * @module structure
 */

import assert from 'assert';

import * as $3Dmol from '3dmol';
import { assignBonds } from './assignBonds';

import { arrayMaxMin, getElement, sendWarning, unreachable } from '../utils';
import { PositioningCallback } from '../utils';
import { Environment, Settings, Structure } from '../dataset';

import { Arrow, CustomShape, Cylinder, Ellipsoid, ShapeData, Sphere, add_shapes } from './shapes';

import { StructureOptions } from './options';

import { COLOR_MAPS } from '../map/colorscales';

const IS_SAFARI =
    navigator.vendor !== undefined &&
    navigator.vendor.indexOf('Apple') > -1 &&
    navigator.userAgent.indexOf('CriOS') === -1 &&
    navigator.userAgent.indexOf('FxiOS') === -1;

/** We need to use the same opacity everywhere since 3Dmol uses 1 opacity per model */
function defaultOpacity(): number {
    // Safari seems to have an issue rendering transparent things with the same
    // settings as FF/Chrome
    return IS_SAFARI ? 0.87 : 0.7;
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

interface ColorBar {
    /** The numbers scaling the color bar for the minimum, the middle
     * and the maximum values respectively displayed from left to right
     * under the color bar */
    min: $3Dmol.Label;
    mid: $3Dmol.Label;
    max: $3Dmol.Label;
    /** The property used in the color bar,
     * displayed under the number labels */
    property: $3Dmol.Label;
    // The color gradient (i.e. the main element) of the color bar
    gradient: $3Dmol.Label;
}

/** Possible options passed to `MoleculeViewer.load` */
export interface LoadOptions {
    /** Supercell to display (default: [1, 1, 1]) */
    supercell: [number, number, number];
    /** Should preserve we the current camera orientation (default: false) */
    keepOrientation: boolean;
    /** Are we loading a file part of a trajectory (default: false) */
    trajectory: boolean;
    /** List of atom-centered environments in the current structure, potentially
     * undefined if the environnement is not part of the dataset.
     * The `structure` field of `Environment` is ignored */
    environments: (Environment | undefined)[];
    /**
     * Index of the environment to highlight, this is only considered if
     * `environments` is set.
     */
    highlight: number;
}

/** */
export class MoleculeViewer {
    /** callback called when a new atom is clicked on */
    public onselect: (atom: number) => void;
    /**
     * Callback to get the initial positioning of the settings modal.
     *
     * The callback is called once, the first time the settings are opened.
     */
    public positionSettingsModal: PositioningCallback;

    /// Shadow root for isolation
    private _shadow: ShadowRoot;
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
    /// Information about the last highlighted environment
    private _lastHighlighted?: {
        /// index of the central atom
        center: number;
        /// user-specified cutoff
        cutoff: number;
    };
    /// Axes shown, if any
    private _axes?: [LabeledArrow, LabeledArrow, LabeledArrow];
    /// Color bar shown, if any
    private _colorBar?: ColorBar;

    /// Representation options from the HTML side
    public _options: StructureOptions;
    /// The supercell used to initialize the viewer
    private _initialSupercell?: [number, number, number];
    // button to reset the environment cutoff to its original value
    private _resetEnvCutoff!: HTMLButtonElement;
    // button to reset reset the supercell
    private _resetSupercell!: HTMLButtonElement;
    /// Show some information on the currently displayed cell to the user
    private _cellInfo: HTMLElement;
    /// Options related to the trajectory
    private _trajectoryOptions: HTMLElement;

    /// Dynamic CSS used to hide options as needed
    private _styles: {
        /// hide options related to unit cell if there is no unit cell in the
        /// current structure
        noCell: CSSStyleSheet;
        /// hide options related to environments if there are no environments
        /// for the current structure
        noEnvs: CSSStyleSheet;
        /// hide options related to shapes if there is no shapes in the
        /// current structure
        noShape: CSSStyleSheet;
        /// hide options related to atom coloring if there are no atom properties
        /// in the current structure
        noProperties: CSSStyleSheet;
    };
    /// List of atom-centered environments for the current structure
    private _environments?: (Environment | undefined)[];
    // List of properties for the current structure
    private _properties?: Record<string, (number | undefined)[]>;
    // Button used to reset the range of color axis
    private _colorReset: HTMLButtonElement;
    // Button used to see more color options
    private _colorMoreOptions: HTMLButtonElement;

    /**
     * Create a new `MoleculeViewer` inside the HTML DOM element with the given `id`.
     *
     * @param element HTML element or HTML id of the DOM element
     *                where the viewer will be created
     * @param propertiesName list of properties keys to be used as options
     */
    constructor(element: string | HTMLElement, propertiesName: string[]) {
        const containerElement = getElement(element);
        const hostElement = document.createElement('div');
        containerElement.appendChild(hostElement);

        hostElement.style.setProperty('height', '100%');
        this._shadow = hostElement.attachShadow({ mode: 'open' });

        this._root = document.createElement('div');
        this._shadow.appendChild(this._root);

        this._root.style.position = 'relative';
        this._root.style.width = '100%';
        this._root.style.height = '100%';

        const viewer = $3Dmol.createViewer(this._root, {
            antialias: true,
            defaultcolors: $3Dmol.elementColors.Jmol,
            backgroundColor: '0xffffff',
            backgroundAlpha: 0,
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
            'bg-light',
            'text-dark'
        );
        this._root.appendChild(this._cellInfo);

        const noCellStyle = new CSSStyleSheet();
        const noEnvsStyle = new CSSStyleSheet();
        const noShapeStyle = new CSSStyleSheet();
        const noPropsStyle = new CSSStyleSheet();

        this._styles = {
            noCell: noCellStyle,
            noEnvs: noEnvsStyle,
            noShape: noShapeStyle,
            noProperties: noPropsStyle,
        };

        this._shadow.adoptedStyleSheets = [
            ...(containerElement.getRootNode() as ShadowRoot).adoptedStyleSheets,
            noCellStyle,
            noEnvsStyle,
            noShapeStyle,
            noPropsStyle,
        ];

        this._options = new StructureOptions(
            this._root,
            (rect) => this.positionSettingsModal(rect),
            propertiesName
        );

        this._options.modal.shadow.adoptedStyleSheets = [
            ...this._options.modal.shadow.adoptedStyleSheets,
            noCellStyle,
            noEnvsStyle,
            noShapeStyle,
            noPropsStyle,
        ];
        this._colorReset = this._options.getModalElement<HTMLButtonElement>('atom-color-reset');
        this._colorMoreOptions =
            this._options.getModalElement<HTMLButtonElement>('atom-color-more-options');

        this._connectOptions();
        this._trajectoryOptions = this._options.getModalElement('trajectory-settings-group');

        // By default, position the modal for settings on top of the viewer,
        // centered horizontally
        this.positionSettingsModal = (rect: DOMRect) => {
            const rootRect = this._root.getBoundingClientRect();
            return {
                left: rootRect.left + rootRect.width / 2 - rect.width / 2,
                top: rootRect.top + 20,
            };
        };

        // Hack to reverse the scroll direction of 3dmol to match that of Plotly
        // The wheel event is captured on the parent of the canvas, modified
        // and then dispatched on the canvas.
        this._root.addEventListener(
            'wheel',
            (event) => {
                // Avoid an infinite loop by only intercepting the original
                // (= trusted) event and not synthetic ones.
                if (event.isTrusted) {
                    event.preventDefault();
                    event.stopImmediatePropagation();

                    const copy = new Event('wheel');

                    Object.assign(copy, {
                        ctrlKey: true, // This emulates a zoom action in the right direction
                        detail: event.detail,
                        wheelDelta: (event as unknown as { wheelDelta: number }).wheelDelta, // Deprecated but used by 3dmol
                        pageX: event.pageX,
                        pageY: event.pageY,
                    });

                    // Dispatch the synthetic event on event.target which
                    // represents the canvas
                    assert(event.target);
                    event.target.dispatchEvent(copy);
                }
            },
            { capture: true }
        );

        window.addEventListener('resize', () => this.resize());
        // waits for loading of widget, then triggers a redraw. fixes some glitches on the Jupyter side
        window.requestAnimationFrame(() => {
            window.dispatchEvent(new Event('resize'));
        });
    }

    /**
     * Recreates structure options
     * @param envView display mode, true if per environments
     * @param propertiesName property names used as the options in the modal
     */
    public refreshOptions(envView: boolean, propertiesName?: string[]): void {
        // Save the current adopted styles
        const adoptedStyleSheets = this._options.modal.shadow.adoptedStyleSheets;

        // Remove the current modal
        this._options.remove();

        // Recreate the StructureOptions
        this._options = new StructureOptions(
            this._root,
            (rect) => this.positionSettingsModal(rect),
            propertiesName
        );

        // Restore the previously saved styles
        this._options.modal.shadow.adoptedStyleSheets = adoptedStyleSheets;

        // Cache the button elements
        this._colorReset = this._options.getModalElement<HTMLButtonElement>('atom-color-reset');
        this._colorMoreOptions =
            this._options.getModalElement<HTMLButtonElement>('atom-color-more-options');
        this._trajectoryOptions = this._options.getModalElement('trajectory-settings-group');

        // Disable as by default, color property is element
        this._colorReset.disabled = true;
        this._colorMoreOptions.disabled = true;

        // Connect event handlers for the new options
        this._connectOptions();

        // Adapt settings to the display mode
        if (envView) {
            this._options.environments.activated.enable();
        } else {
            this._options.environments.activated.disable();
        }
    }

    /**
     * Remove all HTML added by this {@link MoleculeViewer} in the current document
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
        this._updateColorBar();
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
     * @param value list of atom-centered environments for the current structure
     */
    public set environments(value: (Environment | undefined)[] | undefined) {
        this._environments = value;
    }

    /**
     * Load the given `structure` in this viewer.
     *
     * @param structure structure to load
     * @param properties properties to used to load the structure
     * @param options options for the new structure
     * @param onLoadingDone fired when viewer finished rendering
     */
    public load(
        structure: Structure,
        properties?: Record<string, (number | undefined)[]> | undefined,
        options: Partial<LoadOptions> = {},
        onLoadingDone?: () => void
    ): void {
        // if the canvas size changed since last structure, make sure we update
        // everything
        this.resize();

        this._properties = properties;

        let previousDefaultCutoff = undefined;
        if (this._highlighted !== undefined) {
            previousDefaultCutoff = this._defaultCutoff(this._highlighted.center);
        }

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
            this._styles.noCell.replaceSync('.chsp-hide-if-no-cell { display: none; }');
        } else {
            this._styles.noCell.replaceSync('');
        }
        this._showSupercellInfo();

        if (options.trajectory === undefined || !options.trajectory) {
            this._trajectoryOptions.style.display = 'none';
        } else {
            this._trajectoryOptions.style.display = 'block';
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

        if (this._highlighted !== undefined) {
            this._viewer.removeModel(this._highlighted.model);
            this._highlighted = undefined;
        }

        // load new structure
        this._current = {
            model: this._viewer.addModel(),
            structure: structure,
            atomLabels: [],
        };

        setup3DmolStructure(this._current.model, structure);
        const supercell_a = this._options.supercell[0].value;
        const supercell_b = this._options.supercell[1].value;
        const supercell_c = this._options.supercell[2].value;
        if (!(supercell_a === 1 && supercell_b === 1 && supercell_c === 1)) {
            this._viewer.replicateUnitCell(
                this._options.supercell[0].value,
                this._options.supercell[1].value,
                this._options.supercell[2].value,
                this._current.model
            );
        }

        if (this._options.unitCell.value) {
            this._viewer.addUnitCell(this._current.model, {
                box: { color: 'black' },
                astyle: { hidden: true },
                bstyle: { hidden: true },
                cstyle: { hidden: true },
            });
        }
        assignBonds(this._current.model.selectedAtoms({}));

        if (this._environments === undefined) {
            this._styles.noEnvs.replaceSync('.chsp-hide-if-no-environments { display: none; }');
        } else {
            this._styles.noEnvs.replaceSync('');
            assert(this._environments.length === structure.size);
            this._setEnvironmentInteractions();
            const newCenter = options.highlight === undefined ? 0 : options.highlight;
            this._changeHighlighted(newCenter, previousDefaultCutoff);
        }

        if (this._properties === undefined) {
            this._styles.noProperties.replaceSync(
                '.chsp-hide-if-no-atom-properties { display: none; }'
            );
        } else {
            this._styles.noProperties.replaceSync('');
        }
        // update the options for shape visualization
        if (structure.shapes === undefined) {
            this._styles.noShape.replaceSync('.chsp-hide-if-no-shapes { display: none; }');
        } else {
            this._styles.noShape.replaceSync('');
            const selectShape = this._options.getModalElement<HTMLSelectElement>('shapes');
            selectShape.options.length = 0;
            for (const key of Object.keys(structure['shapes'])) {
                selectShape.options.add(new Option(key, key));
            }

            // leave space for up to 3 shapes in the settings, the other one
            // will be accessible with a scrollbar
            selectShape.size = Math.min(Object.keys(structure['shapes']).length, 3);

            // set the highlighting to reflect the options

            this._options.shape.bind(selectShape, 'multival');
        }

        this._updateStyle();

        if (!keepOrientation) {
            this._resetView();
        }

        if (this._options.environments.center.value) {
            this._centerView();
        }

        // make sure to reset axes/labels when the structure changes
        this._options.axes.changed('JS');
        this._options.atomLabels.changed('JS');

        this._viewer.render(() => {
            if (onLoadingDone) {
                onLoadingDone();
            }
        });
    }

    /** Setup interaction (click & hover) for environments highlighting */
    private _setEnvironmentInteractions() {
        if (this._environments === undefined) {
            return;
        }

        // eslint-disable-next-line @typescript-eslint/no-non-null-assertion
        const active = this._environments.filter((e) => e !== undefined).map((e) => e!.center);

        // resets clickable and hoverable state interactions
        for (const atom of this._viewer.selectedAtoms({})) {
            atom.clickable = false;
            atom.hoverable = false;
        }

        assert(this._current !== undefined);
        this._viewer.setClickable(
            { model: this._current.model, index: active },
            true,
            (atom: $3Dmol.AtomSelectionSpec) => this._selectAtom(atom)
        );

        this._setHoverable(active);
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
        let previousDefaultCutoff = undefined;
        if (this._highlighted !== undefined) {
            previousDefaultCutoff = this._defaultCutoff(this._highlighted.center);
        }

        this._changeHighlighted(center, previousDefaultCutoff);
        this._updateStyle();

        const centerView = this._options.environments.center.value;

        if (this._highlighted !== undefined && centerView) {
            if (!this._options.keepOrientation.value) {
                this._resetView();
            }

            this._centerView();
        }

        this._viewer.render();
    }

    /**
     * Applies saved settings, possibly filling in with default values
     */
    public applySettings(settings: Settings): void {
        this._options.applySettings(settings);
    }

    /**
     * Save the values of the current settings in a way that an be used with
     * {@link applySettings} or saved to JSON.
     */
    public saveSettings(): Settings {
        return this._options.saveSettings();
    }

    /**
     * Add the given `callback` to be called whenever a setting changes. The
     * callback will be given the path to the settings as a list of keys; and
     * the new value of the setting.
     *
     * There is currently no way to remove a callback.
     */
    public onSettingChange(callback: (keys: string[], value: unknown) => void): void {
        this._options.onSettingChange(callback);
    }

    /**
     * Returns a PNG screenshot of the viewer as a URI string
     */
    public exportPNG(): string {
        return this._viewer.pngURI();
    }

    /* Set the given list of active atoms as hoverable **/
    public _setHoverable(active: number[]): void {
        assert(this._current !== undefined);
        this._viewer.setHoverDuration(0);

        const n_atoms = this.natoms();
        assert(n_atoms !== undefined);
        const hovered: Map<number, $3Dmol.GLModel> = new Map();

        const hoverStart = (atom: $3Dmol.AtomSpec) => {
            assert(atom.index !== undefined);
            if (hovered.get(atom.index) !== undefined) {
                // the 'hover' model for this atom already exists
                return;
            }

            const model = this._viewer.addModel();
            model.addAtoms([atom]);
            model.setStyle({}, this._centralStyle(0.3));

            // the atom in the new is still marked as clickable, but we only
            // want the atom in the main model to react to interactions
            const sel = model.selectedAtoms({});
            sel[0].clickable = false;
            sel[0].hoverable = false;

            hovered.set(atom.index, model);
            this._viewer.render();
        };

        const hoverStop = (atom: $3Dmol.AtomSpec) => {
            assert(atom.index !== undefined);
            const model = hovered.get(atom.index);
            if (model === undefined) {
                return;
            }

            this._viewer.removeModel(model);
            hovered.delete(atom.index);
            this._viewer.render();
        };

        this._viewer.setHoverable(
            {
                model: this._current.model,
                index: active,
            },
            true,
            hoverStart,
            hoverStop
        );
    }

    private _connectOptions(): void {
        const restyleAndRender = () => {
            this._updateStyle();
            this._viewer.render();
        };

        this._options.spaceFilling.onchange.push(restyleAndRender);
        this._options.bonds.onchange.push(restyleAndRender);
        this._options.atoms.onchange.push(restyleAndRender);

        this._options.atomLabels.onchange.push((showLabels) => {
            if (this._current === undefined) {
                return;
            }

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
                        alignment: new $3Dmol.Vector2(2.0, 2.0),
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
        });

        this._options.axes.onchange.push((value) => this._addAxes(value));

        this._options.rotation.onchange.push((rotate) => {
            if (rotate) {
                this._viewer.spin('vy');
            } else {
                this._viewer.spin(false);
            }
        });

        this._options.unitCell.onchange.push((add) => {
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
        });

        this._options.shape.onchange.push(() => {
            if (this._current === undefined) {
                return;
            }

            this._viewer.removeAllShapes();
            this._updateStyle();
            this._viewer.render();
        });
        const changedSuperCell = () => {
            this._showSupercellInfo();
            if (this._current === undefined) {
                return;
            }

            // remove the atoms previously added by `replicateUnitCell`
            const atoms = this._current.model.selectedAtoms({});
            this._current.model.removeAtoms(atoms.splice(this._current.structure.size));

            const supercell_a = this._options.supercell[0].value;
            const supercell_b = this._options.supercell[1].value;
            const supercell_c = this._options.supercell[2].value;
            if (!(supercell_a === 1 && supercell_b === 1 && supercell_c === 1)) {
                this._viewer.replicateUnitCell(
                    this._options.supercell[0].value,
                    this._options.supercell[1].value,
                    this._options.supercell[2].value,
                    this._current.model
                );
            }
            assignBonds(this._current.model.selectedAtoms({}));

            if (this._highlighted !== undefined) {
                const previousDefaultCutoff = this._defaultCutoff(this._highlighted.center);
                this._changeHighlighted(this._highlighted.center, previousDefaultCutoff);
            }
            if (this._environmentsEnabled()) {
                this._setEnvironmentInteractions();
            }
            this._updateStyle();
            this._viewer.render();
        };
        this._options.supercell[0].onchange.push(changedSuperCell);
        this._options.supercell[1].onchange.push(changedSuperCell);
        this._options.supercell[2].onchange.push(changedSuperCell);

        this._options.environments.bgColor.onchange.push(restyleAndRender);
        this._options.environments.bgStyle.onchange.push(restyleAndRender);
        this._options.environments.cutoff.onchange.push(() => {
            if (this._highlighted !== undefined) {
                const previousDefaultCutoff = this._defaultCutoff(this._highlighted.center);
                this._changeHighlighted(this._highlighted.center, previousDefaultCutoff);
            }
            restyleAndRender();
        });
        this._options.environments.center.onchange.push((center) => {
            if (!this._options.keepOrientation.value) {
                this._resetView();
            }

            if (center) {
                this._centerView();
            }

            this._viewer.render();
        });

        // Deal with activation/de-activation of environments
        this._options.environments.activated.onchange.push((value) => {
            this._enableEnvironmentSettings(value);

            if (value) {
                assert(this._lastHighlighted !== undefined);
                // add the 3DMol model for highlighted environment
                this._options.environments.cutoff.value = this._lastHighlighted.cutoff;
                const center = this._lastHighlighted.center;
                const previousDefaultCutoff = this._defaultCutoff(center);
                this._changeHighlighted(center, previousDefaultCutoff);
            } else {
                assert(this._highlighted !== undefined);
                // remove the 3DMol model for highlighted environment
                this._viewer.removeModel(this._highlighted.model);
                // keep information about the last highlighted point around
                this._lastHighlighted = {
                    center: this._highlighted.center,
                    cutoff: this._options.environments.cutoff.value,
                };
            }

            restyleAndRender();
        });

        // ======= color settings
        // setup state when the property changes
        const colorPropertyChanged = () => {
            const property = this._options.color.property.value;

            if (property !== 'element') {
                this._options.color.transform.enable();
                this._options.color.transform.value = 'linear';
                this._options.color.min.enable();
                this._options.color.max.enable();

                this._colorReset.disabled = false;
                this._colorMoreOptions.disabled = false;
                this._options.color.palette.enable();

                const values = this._colorValues(property, 'linear');

                if (values.some((v) => v === null)) {
                    sendWarning(
                        'The selected structure has undefined properties for some atoms, these atoms will be colored in light gray.'
                    );
                }

                // To change min and max values when the transform has been changed
                const { max, min } = arrayMaxMin(values);

                // We have to set max first and min second here to avoid sending
                // a spurious warning in `colorRangeChange` below in case the
                // new min is bigger than the old max.
                this._options.color.min.value = Number.NEGATIVE_INFINITY;
                this._options.color.max.value = max;
                this._options.color.min.value = min;
                this._setScaleStep([min, max]);
            } else {
                this._options.color.transform.disable();
                this._options.color.min.disable();
                this._options.color.max.disable();

                this._colorReset.disabled = true;
                this._colorMoreOptions.disabled = true;
                this._options.color.palette.disable();

                this._viewer.setColorByElement({}, $3Dmol.elementColors.Jmol);
            }

            this._updateColorBar();
            restyleAndRender();
        };
        this._options.color.property.onchange.push(colorPropertyChanged);

        const colorRangeChange = (minOrMax: 'min' | 'max') => {
            const min = this._options.color.min.value;
            const max = this._options.color.max.value;
            if (min > max) {
                sendWarning(
                    `The inserted min and max values in color are such that min > max! The last inserted value was reset.`
                );
                if (minOrMax === 'min') {
                    this._options.color.min.reset();
                } else {
                    this._options.color.max.reset();
                }
                return;
            }
            this._setScaleStep([min, max]);
            this._updateColorBar();
            restyleAndRender();
        };

        // ======= color transform
        this._options.color.transform.onchange.push(() => {
            const property = this._options.color.property.value;
            assert(property !== 'element');
            const transform = this._options.color.transform.value;

            const values = this._colorValues(property, transform);
            // To change min and max values when the transform has been changed
            const { min, max } = arrayMaxMin(values);

            // to avoid sending a spurious warning in `colorRangeChange` below
            // in case the new min is bigger than the old max.
            this._options.color.min.value = Number.NEGATIVE_INFINITY;
            this._options.color.max.value = max;
            this._options.color.min.value = min;
            this._setScaleStep([min, max]);
            this._updateColorBar();
            restyleAndRender();
        });

        // ======= color min
        this._options.color.min.onchange.push(() => {
            colorRangeChange('min');
            restyleAndRender();
        });

        // ======= color max
        this._options.color.max.onchange.push(() => {
            colorRangeChange('max');
            restyleAndRender();
        });

        // ======= color reset
        this._colorReset.addEventListener('click', () => {
            const properties = JSON.parse(JSON.stringify(this._properties)) as Record<
                string,
                (number | undefined)[]
            >;
            const property: string = this._options.color.property.value;
            // Use map to extract the specified property values into an array
            const values: number[] = properties[property].filter(
                (value) => !isNaN(Number(value))
            ) as number[];
            // To change min and max values when the transform has been changed
            const [min, max]: [number, number] = [Math.min(...values), Math.max(...values)];
            this._options.color.min.value = min;
            this._options.color.max.value = max;
            this._setScaleStep([min, max]);
            this._updateColorBar();
            restyleAndRender();
        });

        // ======= color palette
        this._options.color.palette.onchange.push(() => {
            this._updateColorBar();
            restyleAndRender();
        });

        // Setup various buttons
        this._resetEnvCutoff = this._options.getModalElement<HTMLButtonElement>('env-reset');
        this._resetEnvCutoff.onclick = () => {
            assert(this._highlighted !== undefined);
            this._options.environments.cutoff.value = this._defaultCutoff(this._highlighted.center);
            restyleAndRender();
        };

        const alignX = this._options.getModalElement<HTMLButtonElement>('align-x');
        alignX.onclick = () => this._viewAlong([1, 0, 0]);

        const alignY = this._options.getModalElement<HTMLButtonElement>('align-y');
        alignY.onclick = () => this._viewAlong([0, 1, 0]);

        const alignZ = this._options.getModalElement<HTMLButtonElement>('align-z');
        alignZ.onclick = () => this._viewAlong([0, 0, 1]);

        const alignA = this._options.getModalElement<HTMLButtonElement>('align-a');
        alignA.onclick = () => {
            if (this._current === undefined || this._current.structure.cell === undefined) {
                return;
            }
            const cell = this._current.structure.cell;
            const a: [number, number, number] = [cell[0], cell[1], cell[2]];
            this._viewAlong(a);
        };

        const alignB = this._options.getModalElement<HTMLButtonElement>('align-b');
        alignB.onclick = () => {
            if (this._current === undefined || this._current.structure.cell === undefined) {
                return;
            }
            const cell = this._current.structure.cell;
            const b: [number, number, number] = [cell[3], cell[4], cell[5]];
            this._viewAlong(b);
        };

        const alignC = this._options.getModalElement<HTMLButtonElement>('align-c');
        alignC.onclick = () => {
            if (this._current === undefined || this._current.structure.cell === undefined) {
                return;
            }
            const cell = this._current.structure.cell;
            const c: [number, number, number] = [cell[6], cell[7], cell[8]];
            this._viewAlong(c);
        };

        this._resetSupercell = this._options.getModalElement<HTMLButtonElement>('reset-supercell');
        this._resetSupercell.onclick = () => {
            assert(this._initialSupercell !== undefined);
            this._options.supercell[0].value = this._initialSupercell[0];
            this._options.supercell[1].value = this._initialSupercell[1];
            this._options.supercell[2].value = this._initialSupercell[2];
        };

        // Reset zoom level and centering when double clicked
        this._root.ondblclick = () => {
            this._resetView();

            if (this._options.environments.center.value) {
                this._centerView();
            }

            this._viewer.render();
        };
    }

    /**
     * Function called whenever the user click on an atom in the viewer
     */
    private _selectAtom(atom: $3Dmol.AtomSelectionSpec): void {
        // use atom.serial instead of atom.index to ensure we are getting the
        // id of the atom inside the central cell when using a supercell
        assert(atom.serial !== undefined);

        this.onselect(atom.serial);
    }

    /**
     * Update the styles of all atoms as required
     */
    private _updateStyle(): void {
        if (this._current === undefined) {
            return;
        }

        if (this._options.shape.value !== '') {
            this._addShapes();
        }

        // if there is no environment to highlight, render all atoms with the
        // main style
        if (!this._environmentsEnabled()) {
            this._current.model.setStyle({}, this._mainStyle());
            return;
        }

        assert(this._highlighted !== undefined);
        // otherwise, render all atoms with hidden/background style
        this._current.model.setStyle({}, this._hiddenStyle());
        // and the central atom with central/highlighted style
        this._current.model.setStyle(
            { index: this._highlighted.center },
            this._centralStyle(0.4),
            /* add */ true
        );

        // and the environment around the central atom with main style
        this._highlighted.model.setStyle({}, this._mainStyle());
    }

    private _addAxes(value: string): void {
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
    }

    private _addShapes(): void {
        if (this._current === undefined) {
            return;
        }
        if (this._options.shape.value === '') {
            return;
        }

        this._viewer.removeAllShapes();

        // removeAllShapes also removes the unit cell, so let's add it back
        if (this._options.unitCell.value) {
            this._viewer.addUnitCell(this._current.model, {
                box: { color: 'black' },
                astyle: { hidden: true },
                bstyle: { hidden: true },
                cstyle: { hidden: true },
            });
        }

        // removeAllShapes also removes the axes, so let's add them back
        this._addAxes(this._options.axes.value);

        assert(this._current.atomLabels.length === 0);

        const structure = this._current.structure;
        assert(!(structure.shapes === undefined));

        const active_shapes = this._options.shape.value.split(',');

        // consolidates all shapes in a single CustomShapeSpec
        const all_shapes: $3Dmol.CustomShapeSpec = {
            vertexArr: [],
            normalArr: [],
            faceArr: [],
            color: [],
        };

        // color function for atoms
        const colorfunc = this._colorFunction();
        const atomspec: $3Dmol.AtomSpec = {};

        for (const shape of active_shapes) {
            if (shape === '') {
                continue;
            }
            assert(shape in structure.shapes);
            const current_shape = structure.shapes[shape];
            const supercell_a = this._options.supercell[0].value;
            const supercell_b = this._options.supercell[1].value;
            const supercell_c = this._options.supercell[2].value;
            let cell = this._current.structure.cell;

            if ((supercell_a > 1 || supercell_b > 1 || supercell_c > 1) && cell === undefined) {
                return;
            } else if (cell === undefined) {
                cell = [1, 0, 0, 0, 1, 0, 0, 0, 1];
            }

            for (let a = 0; a < supercell_a; a++) {
                for (let b = 0; b < supercell_b; b++) {
                    for (let c = 0; c < supercell_c; c++) {
                        let shape_data: Partial<ShapeData> = { ...current_shape.parameters.global };

                        shape_data.position = [0, 0, 0]; // defaults
                        if (current_shape.parameters.structure) {
                            shape_data = {
                                ...shape_data,
                                ...current_shape.parameters.structure[0],
                            };
                        }
                        if (current_shape.parameters.atom) {
                            for (let i = 0; i < structure.size; i++) {
                                const name = structure.names[i];
                                assert(i < current_shape.parameters.atom.length);
                                const atom_pars = {
                                    ...shape_data,
                                    ...current_shape.parameters.atom[i],
                                };

                                let position: [number, number, number] = [
                                    structure.x[i],
                                    structure.y[i],
                                    structure.z[i],
                                ];

                                // only overrides the atom position if it's given explicitly
                                const atom_position = current_shape.parameters.atom[i].position;
                                if (atom_position !== undefined) {
                                    position = atom_position;
                                }

                                // adds supercell offset
                                position[0] += a * cell[0] + b * cell[3] + c * cell[6];
                                position[1] += a * cell[1] + b * cell[4] + c * cell[7];
                                position[2] += a * cell[2] + b * cell[5] + c * cell[8];

                                atom_pars.position = position;
                                // obey explicit color specification if given,
                                // otherwise color as the corresponding atom
                                if (!atom_pars.color) {
                                    if (colorfunc) {
                                        atomspec.serial = i;
                                        atom_pars.color = colorfunc(atomspec);
                                    } else {
                                        atom_pars.color = $3Dmol.elementColors.Jmol[name];
                                    }
                                }

                                if (current_shape.kind === 'sphere') {
                                    const shape = new Sphere(atom_pars);
                                    //this._viewer.addCustom(
                                    add_shapes(
                                        all_shapes,
                                        shape.outputTo3Dmol(atom_pars.color || 0xffffff),
                                        this._viewer,
                                        32768
                                    );
                                } else if (current_shape.kind === 'ellipsoid') {
                                    const shape = new Ellipsoid(atom_pars);
                                    add_shapes(
                                        all_shapes,
                                        shape.outputTo3Dmol(atom_pars.color || 0xffffff),
                                        this._viewer,
                                        32768
                                    );
                                } else if (current_shape.kind === 'cylinder') {
                                    const shape = new Cylinder(atom_pars);
                                    add_shapes(
                                        all_shapes,
                                        shape.outputTo3Dmol(atom_pars.color || 0xffffff),
                                        this._viewer,
                                        32768
                                    );
                                } else if (current_shape.kind === 'arrow') {
                                    const shape = new Arrow(atom_pars);
                                    add_shapes(
                                        all_shapes,
                                        shape.outputTo3Dmol(atom_pars.color || 0xffffff),
                                        this._viewer,
                                        32768
                                    );
                                } else {
                                    assert(current_shape.kind === 'custom');
                                    const shape = new CustomShape(atom_pars);
                                    add_shapes(
                                        all_shapes,
                                        shape.outputTo3Dmol(atom_pars.color || 0xffffff),
                                        this._viewer,
                                        32768
                                    );
                                }
                            }
                        } else {
                            // the shape is defined only at the structure level
                            if (current_shape.kind === 'sphere') {
                                const shape = new Sphere(shape_data);
                                add_shapes(
                                    all_shapes,
                                    shape.outputTo3Dmol(shape_data.color || 0xffffff),
                                    this._viewer,
                                    32768
                                );
                            } else if (current_shape.kind === 'ellipsoid') {
                                const shape = new Ellipsoid(shape_data);
                                add_shapes(
                                    all_shapes,
                                    shape.outputTo3Dmol(shape_data.color || 0xffffff),
                                    this._viewer,
                                    32768
                                );
                            } else if (current_shape.kind === 'arrow') {
                                const shape = new Arrow(shape_data);
                                add_shapes(
                                    all_shapes,
                                    shape.outputTo3Dmol(shape_data.color || 0xffffff),
                                    this._viewer,
                                    32768
                                );
                            } else {
                                assert(current_shape.kind === 'custom');
                                const shape = new CustomShape(shape_data);
                                add_shapes(
                                    all_shapes,
                                    shape.outputTo3Dmol(shape_data.color || 0xffffff),
                                    this._viewer,
                                    32768
                                );
                            }
                        }
                    }
                }
            }
        }

        if (Array.isArray(all_shapes.faceArr) && all_shapes.faceArr.length > 0) {
            // adds all shapes that have been accumulated
            this._viewer.addCustom(all_shapes);
        }
        this._viewer.render();
    }
    /**
     * Get the main style used for all atoms/atoms inside the environment when
     * highlighting a specific environment
     */
    private _mainStyle(): Partial<$3Dmol.AtomStyleSpec> {
        const style: Partial<$3Dmol.AtomStyleSpec> = {};
        if (this._options.atoms.value) {
            style.sphere = {
                scale: this._options.spaceFilling.value ? 1.0 : 0.22,
                colorfunc: this._colorFunction(),
            } as unknown as $3Dmol.SphereStyleSpec;
        }
        if (this._options.bonds.value) {
            style.stick = {
                radius: 0.15,
                colorfunc: this._colorFunction(),
            } as unknown as $3Dmol.StickStyleSpec;
        }

        return style;
    }

    /**
     * Get the values that should be used to color atoms when coloring them
     * according to properties
     */
    private _colorValues(property: string, transform: string): Array<number | null> {
        assert(this._properties !== undefined);
        assert(Object.keys(this._properties).includes(property));

        // JSON.parse & JSON.stringify to make a deep copy of the properties to
        // avoid modifying the original ones.
        //
        // This also transforms all `undefined` values into `null`
        let currentProperty = JSON.parse(JSON.stringify(this._properties[property])) as Array<
            number | null
        >;

        let transformFunc = (x: number) => x;
        if (transform === 'log') {
            transformFunc = Math.log10;
        } else if (transform === 'sqrt') {
            transformFunc = Math.sqrt;
        } else if (transform === 'inverse') {
            transformFunc = (x) => 1 / x;
        }

        currentProperty = currentProperty.map((value) => {
            if (value !== null && !isNaN(value)) {
                return transformFunc(value);
            } else {
                return value;
            }
        });

        return currentProperty;
    }

    /**
     * Get a function computing the atom color, that can be used as 3Dmol
     * `colorfunc`
     */
    private _colorFunction(): ((atom: $3Dmol.AtomSpec) => number) | undefined {
        if (this._properties === undefined) {
            return undefined;
        }

        const property = this._options.color.property.value;
        if (property === 'element') {
            return undefined;
        }

        const transform = this._options.color.transform.value;

        const currentProperty = this._colorValues(property, transform);
        const min = this._options.color.min.value;
        const max = this._options.color.max.value;

        const palette = COLOR_MAPS[this._options.color.palette.value];
        const colors = [];
        // NB: 3Dmol does not consider midpoints. Palettes should be given as
        // uniformly spaced colors, otherwise the mapping will be incorrect
        for (let c = 0; c < palette.length; c++) {
            colors.push(palette[c][1]);
        }
        const grad: $3Dmol.Gradient = new $3Dmol.Gradient.CustomLinear(min, max, colors);

        return (atom: $3Dmol.AtomSpec) => {
            assert(atom.serial !== undefined);
            const value = currentProperty[atom.serial];
            if (value === null) {
                // missing values
                return 0xdddddd;
            } else if (isNaN(value)) {
                // NaN values
                return 0x222222;
            } else {
                if (Number.isFinite(min) && Number.isFinite(max)) {
                    return grad.valueToHex(value);
                } else {
                    return 0x222222;
                }
            }
        };
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
                // slightly smaller radius than the main style
                radius: 0.149,
                opacity: defaultOpacity(),
                hidden: !this._options.bonds.value,
            };

            if (bgStyle === 'ball-stick' && this._options.shape.value === '') {
                style.sphere = {
                    // slightly smaller scale than the main style
                    scale: this._options.spaceFilling.value ? 0.999 : 0.219,
                    opacity: defaultOpacity(),
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
                style.stick.color = 0x808080;
            }

            if (style.sphere !== undefined) {
                style.sphere.color = 0x808080;
            }
        } else if (bgColor === 'property') {
            if (style.stick !== undefined) {
                style.stick = {
                    // slightly smaller radius than the main style
                    radius: 0.149,
                    opacity: defaultOpacity(),
                    hidden: !this._options.bonds.value,
                    colorfunc: this._colorFunction(),
                } as unknown as $3Dmol.StickStyleSpec;
            }

            if (style.sphere !== undefined) {
                style.sphere = {
                    // slightly smaller scale than the main style
                    scale: this._options.spaceFilling.value ? 0.999 : 0.219,
                    opacity: defaultOpacity(),
                    colorfunc: this._colorFunction(),
                } as unknown as $3Dmol.SphereStyleSpec;
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
    private _centralStyle(scale: number): Partial<$3Dmol.AtomStyleSpec> {
        return {
            sphere: {
                scale: scale,
                color: 'green',
                opacity: defaultOpacity(),
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
        const toggle = this._options.getModalElement(`env-activated`)
            .nextElementSibling as HTMLLabelElement;
        assert(toggle !== null);
        const reset = this._resetEnvCutoff;

        if (show) {
            if (this._environmentsEnabled()) {
                // nothing to do
                return;
            }

            reset.disabled = false;
            toggle.innerText = 'Disable';

            this._options.environments.cutoff.enable();
            this._options.environments.bgStyle.enable();
            this._options.environments.bgColor.enable();
        } else {
            if (!this._environmentsEnabled()) {
                // nothing to do
                return;
            }

            reset.disabled = true;
            toggle.innerText = 'Enable';

            this._options.environments.cutoff.disable();
            this._options.environments.bgStyle.disable();
            this._options.environments.bgColor.disable();
        }
    }

    /**
     * Get the cutoff for the environment around the given `center`
     */
    private _defaultCutoff(center: number): number {
        assert(this._environments !== undefined);
        const environment = this._environments[center];
        assert(environment !== undefined);
        return environment.cutoff;
    }

    /**
     * Change which central atom is highlighted in the system to `center`. If
     * `center` is undefined, this disable highlighting.
     *
     * @param center index of the atom to highlight
     */
    private _changeHighlighted(center?: number, previousDefaultCutoff?: number): void {
        if (this._highlighted !== undefined) {
            this._viewer.removeModel(this._highlighted.model);
        }

        if (center === undefined) {
            this._options.environments.cutoff.value = 0;
            this._highlighted = undefined;
        } else {
            if (this._environments === undefined) {
                throw Error('can not highlight an atom without having a list of environments');
            }

            const environment = this._environments[center];
            if (environment === undefined) {
                throw Error(
                    `can not highlight atom ${center}: it is not part of the list of environments`
                );
            }

            // only change the cutoff if it was not changed manually by the user
            const htmlCutoff = this._options.environments.cutoff.value;
            if (previousDefaultCutoff === undefined || previousDefaultCutoff === htmlCutoff) {
                this._options.environments.cutoff.value = environment.cutoff;
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
            // initialize with main style
            this._highlighted.model.setStyle({}, this._mainStyle());
        }
    }

    /**
     * Reset the view by re-centering it and zooming to fit the model as much as
     * possible inside the views.
     */
    private _resetView(): void {
        this._viewer.zoomTo();
        this._viewer.zoom(2.0);
        this._viewer.setSlab(-1000, 1000);
    }

    /**
     * Centers the view around the selected atom (if there is one)
     */
    private _centerView(): void {
        if (this._highlighted !== undefined && this._current !== undefined) {
            // use index rather than serial to specify the selection, to avoid picking also the
            // periodic replicas. however we then have to specify the model id, otherwise it'd
            // pick the index from the highlighted selection, which does not match the serial ID
            this._viewer.center({ index: this._highlighted.center, model: this._current.model });
        }
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

    /** Changes the step of the arrow buttons in min/max input based on dataset range*/
    private _setScaleStep(axisBounds: number[]): void {
        if (axisBounds !== undefined) {
            // round to 10 decimal places so it does not break in Firefox
            const step = Math.round(((axisBounds[1] - axisBounds[0]) / 20) * 10 ** 10) / 10 ** 10;
            const minElement = this._options.getModalElement<HTMLInputElement>(`atom-color-min`);
            const maxElement = this._options.getModalElement<HTMLInputElement>(`atom-color-max`);
            minElement.step = `${step}`;
            maxElement.step = `${step}`;
        }
    }

    // this is duplicated from map.ts - should really move all this stuff in a utils module
    // and uniform the names and modes of operation
    private _colorTitle(): string {
        let title = this._options.color.property.value;
        switch (this._options.color.transform.value) {
            case 'inverse':
                title = `(${title})⁻¹`;
                break;
            case 'log':
                title = `log<sub>10</sub>(${title})`;
                break;
            case 'sqrt':
                title = `&#x221A;(${title})`;
                break;
            case 'linear':
                break;
            default:
                break;
        }
        return title;
    }

    /**
     * Add a color bar corresponding to the color palette and the atom
     * properties displaying the minimum and maximum values
     */
    private _addColorBar(): ColorBar {
        const palette = this._options.color.palette.value;
        const title = this._colorTitle();

        const viewerWidth = this._viewer.container?.clientWidth;
        const viewerHeight = this._viewer.container?.clientHeight;
        assert(viewerWidth !== undefined && viewerHeight !== undefined);
        const width = 20;
        const horizontalShift = 5;

        const titleHeight = 30;
        const height = viewerHeight - (titleHeight + 50);

        // Create the color bar gradient using canvas
        const gradientCanvas = document.createElement('canvas');

        gradientCanvas.width = width;
        gradientCanvas.height = height;
        const mapColor = COLOR_MAPS[palette];

        const ctx = gradientCanvas.getContext('2d');

        if (!ctx) {
            throw new Error('Could not get 2D context from canvas');
        }

        // Create a linear gradient from colorStart to colorEnd
        const gradient = ctx.createLinearGradient(0, height, 0, 0);
        for (let c = 0; c < mapColor.length; c++) {
            gradient.addColorStop(mapColor[c][0], mapColor[c][1]);
        }

        ctx.fillStyle = gradient;
        ctx.fillRect(0, 0, width, height);

        const min = this._options.color.min.value;
        const max = this._options.color.max.value;
        const mid = (min + max) / 2;

        return {
            min: this._viewer.addLabel(
                fixedWidthFloat(min) + '-',
                this._colorBarLabelsSpec({
                    x: viewerWidth - (horizontalShift + width - 5),
                    y: viewerHeight - (titleHeight + 8),
                    alignment: 'centerRight',
                    fontSize: 14,
                })
            ),
            mid: this._viewer.addLabel(
                fixedWidthFloat(mid) + '-',
                this._colorBarLabelsSpec({
                    x: viewerWidth - (horizontalShift + width - 5),
                    y: viewerHeight - (titleHeight + height / 2),
                    alignment: 'centerRight',
                    fontSize: 14,
                })
            ),
            max: this._viewer.addLabel(
                fixedWidthFloat(max) + '-',
                this._colorBarLabelsSpec({
                    x: viewerWidth - (horizontalShift + width - 5),
                    y: viewerHeight - (titleHeight + height - 8),
                    alignment: 'centerRight',
                    fontSize: 14,
                })
            ),
            property: this._viewer.addLabel(
                title,
                this._colorBarLabelsSpec({
                    x: viewerWidth - 6,
                    y: viewerHeight - 15,
                    alignment: 'centerRight',
                    fontSize: 16,
                })
            ),
            gradient: this._viewer.addLabel('.', {
                position: new $3Dmol.Vector3(
                    viewerWidth - (width + horizontalShift),
                    viewerHeight - (height + titleHeight),
                    0
                ),
                backgroundImage: gradientCanvas,
                borderColor: 'black',
                borderOpacity: 1,
                borderThickness: 1,
                fontColor: 'white',
                fontSize: 2,
                fontOpacity: 0,
                inFront: true,
                useScreen: true,
            }),
        };
    }

    /** Generate a LabelSpec for the key values in the color bar */
    private _colorBarLabelsSpec(options: {
        x: number;
        y: number;
        alignment: string;
        fontSize: number;
    }): $3Dmol.LabelSpec {
        return {
            alignment: options.alignment,
            position: new $3Dmol.Vector3(options.x, options.y, 0),
            font: 'sans',
            fontColor: 'black',
            fontSize: options.fontSize,
            showBackground: false,
            inFront: true,
            useScreen: true,
            backgroundColor: 'red',
        };
    }

    private _deleteColorBar() {
        if (this._colorBar !== undefined) {
            this._viewer.removeLabel(this._colorBar.min);
            this._viewer.removeLabel(this._colorBar.mid);
            this._viewer.removeLabel(this._colorBar.max);
            this._viewer.removeLabel(this._colorBar.property);
            this._viewer.removeLabel(this._colorBar.gradient);
            this._colorBar = undefined;
        }
    }

    private _updateColorBar() {
        this._deleteColorBar();
        if (this._options.color.property.value !== 'element') {
            this._colorBar = this._addColorBar();
        }
    }
}

/** Try to format a float to ensure it fits in a 4 characters limit */
function fixedWidthFloat(value: number): string {
    if (value === 0) {
        return '0';
    } else if (value > 0) {
        // we can use the full 4 characters for the number
        if (value < 10000 && value > 0.009) {
            if (Number.isInteger(value)) {
                return value.toString();
            } else {
                // convert to fixed format with 3 decimals, then truncate to fit
                // in 4 characters
                let result = value.toFixed(3).substring(0, 4);
                if (result[3] === '.') {
                    result = result.substring(0, 3);
                }
                return result;
            }
        } else {
            // convert to exponential format with no decimal part
            return value.toExponential(0).replace('+', '');
        }
    } else {
        // we only have 3 characters, we need one for the '-' sign
        if (value > -1000 && value < -0.09) {
            if (Number.isInteger(value)) {
                return value.toString();
            } else {
                // convert to fixed format with 2 decimals, then truncate to fit in
                // 4 characters
                let result = value.toFixed(3).substring(0, 4);
                if (result[3] === '.') {
                    result = result.substring(0, 3);
                }
                return result;
            }
        } else {
            // convert to exponential format with no decimal part
            return value.toExponential(0).replace('+', '');
        }
    }
}
