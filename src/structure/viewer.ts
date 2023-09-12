/**
 * @packageDocumentation
 * @module structure
 */

import assert from 'assert';

import * as $3Dmol from '3dmol';
import { assignBonds } from './assignBonds';

import { getElement, sendWarning, unreachable } from '../utils';
import { PositioningCallback } from '../utils';
import { Environment, Property, Settings, Structure } from '../dataset';

import {
    CustomShape,
    CustomShapeData,
    Ellipsoid,
    EllipsoidData,
    Sphere,
    SphereData,
} from './shapes';

import { MapData } from '../map/data';
import { EnvironmentIndexer } from '../indexer';
import { StructureOptions } from './options';

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
function setup3DmolStructure(model: $3Dmol.GLModel, structure: Structure, properties?: Record<string, number | undefined>[] | undefined): void {
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
        
        if (properties !== undefined) {
            atoms.push({
                serial: i,
                elem: structure.names[i],
                properties: properties[i],
                x: x,
                y: y,
                z: z,
            });
        } else {
            atoms.push({
                serial: i,
                elem: structure.names[i],
                x: x,
                y: y,
                z: z,
            });
        }
    }

    model.addAtoms(atoms);
}

interface LabeledArrow {
    label: $3Dmol.Label;
    arrow: $3Dmol.GLShape;
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
    /// Axes shown, if any
    private _axes?: [LabeledArrow, LabeledArrow, LabeledArrow];

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
    };
    /// List of atom-centered environments for the current structure
    private _environments?: (Environment | undefined)[];
    // List of properties for the current structure
    private _properties?: Record<string, number | undefined>[] | undefined;
    // All known properties
    private _data: MapData;
    // environment indexer
    private _indexer: EnvironmentIndexer;
    /**
     * Create a new `MoleculeViewer` inside the HTML DOM element with the given `id`.
     *
     * @param element HTML element or HTML id of the DOM element
     *                where the viewer will be created
     */
    constructor(
        element: string | HTMLElement,
        indexer: EnvironmentIndexer,
        properties: { [name: string]: Property }) {
        const containerElement = getElement(element);
        const hostElement = document.createElement('div');
        containerElement.appendChild(hostElement);

        hostElement.style.setProperty('height', '100%');
        this._shadow = hostElement.attachShadow({ mode: 'open' });

        this._root = document.createElement('div');
        this._shadow.appendChild(this._root);

        this._root.style.position = 'relative';
        // this._root.id = this.guid;
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

        this._styles = {
            noCell: noCellStyle,
            noEnvs: noEnvsStyle,
            noShape: noShapeStyle,
        };

        this._shadow.adoptedStyleSheets = [
            ...(containerElement.getRootNode() as ShadowRoot).adoptedStyleSheets,
            noCellStyle,
            noEnvsStyle,
            noShapeStyle,
        ];

        // Options reuse the same style sheets so they must be created after these.
        this._indexer = indexer;
        this._data = new MapData(properties);
        this._options = new StructureOptions(
            this._root, 
            this._data,
            (rect) => this.positionSettingsModal(rect)
        );

        this._options.modal.shadow.adoptedStyleSheets = [
            ...this._options.modal.shadow.adoptedStyleSheets,
            noCellStyle,
            noEnvsStyle,
            noShapeStyle,
        ];

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
    public load(structure: Structure, properties: Record<string, number | undefined>[] | undefined, options: Partial<LoadOptions> = {}): void {
        // if the canvas size changed since last structure, make sure we update
        // everything
        this.resize();

        let previousDefaultCutoff = undefined;
        if (this._highlighted !== undefined) {
            previousDefaultCutoff = this._defaultCutoff(this._highlighted.center);
        }

        // Deal with loading options
        this._environments = options.environments;
        this._properties = properties;

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
        setup3DmolStructure(this._current.model, structure, properties);
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

            this._options.shape.bind(selectShape, 'options');
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

        this._viewer.render();
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

        this._options.axes.onchange.push((value) => {
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
        });

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

            this._viewer.replicateUnitCell(
                this._options.supercell[0].value,
                this._options.supercell[1].value,
                this._options.supercell[2].value,
                this._current.model
            );
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
            restyleAndRender();
        });

        this._options.color.property.onchange.push(() => {
            if (this._options.color.property.value !== 'element') {
                this._options.color.map.enable();
                if (this._properties !== undefined) {
                    if (this._properties.some((record) => Object.values(record).some((v) => v === undefined))) {
                        sendWarning("The selected structure has undefined properties for some atoms, these atoms will be colored in grey.");
                    }
                };
            } else {
                this._options.color.map.disable();
                // const sel: Partial<$3Dmol.AtomStyleSpec> = {};
                this._viewer.setColorByElement({}, $3Dmol.elementColors.Jmol);
            }
            restyleAndRender();
        });
        this._options.color.map.onchange.push(restyleAndRender);

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
            this._viewer.render();
            return;
        }

        assert(this._highlighted !== undefined);
        // otherwise, render all atoms with hidden/background style
        this._current.model.setStyle({}, this._hiddenStyle());
        // the central atom with central/highlighted style
        this._current.model.setStyle(
            { index: this._highlighted.center },
            this._centralStyle(0.4),
            /* add: */ true
        );

        // and the environment around the central atom with main style
        this._highlighted.model.setStyle({}, this._mainStyle());
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

        assert(this._current.atomLabels.length === 0);

        const structure = this._current.structure;
        assert(!(structure.shapes === undefined));

        const active_shapes = this._options.shape.value.split(',');

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
                        for (let i = 0; i < structure.size; i++) {
                            const name = structure.names[i];
                            const position: [number, number, number] = [
                                structure.x[i] + a * cell[0] + b * cell[3] + c * cell[6],
                                structure.y[i] + a * cell[1] + b * cell[4] + c * cell[7],
                                structure.z[i] + a * cell[2] + b * cell[5] + c * cell[8],
                            ];

                            if (current_shape[i].kind === 'ellipsoid') {
                                const data = current_shape[i] as unknown as EllipsoidData;
                                const shape = new Ellipsoid(position, data);
                                this._viewer.addCustom(
                                    shape.outputTo3Dmol($3Dmol.elementColors.Jmol[name] || 0x000000)
                                );
                            } else if (current_shape[i].kind === 'custom') {
                                const data = current_shape[i] as unknown as CustomShapeData;
                                const shape = new CustomShape(position, data);
                                this._viewer.addCustom(
                                    shape.outputTo3Dmol($3Dmol.elementColors.Jmol[name] || 0x000000)
                                );
                            } else {
                                assert(current_shape[i].kind === 'sphere');
                                const data = current_shape[i] as unknown as SphereData;
                                const shape = new Sphere(position, data);
                                this._viewer.addCustom(
                                    shape.outputTo3Dmol($3Dmol.elementColors.Jmol[name] || 0x000000)
                                );
                            }
                        }
                    }
                }
            }
        }

        this._viewer.render();
    }
    /**
     * Get the main style used for all atoms/atoms inside the environment when
     * highlighting a specific environment
     */
    private _mainStyle(): Partial<$3Dmol.AtomStyleSpec> {
        const [min, max] = $3Dmol.getPropertyRange(this._current?.model.selectedAtoms({}), this._options.color.property.value) as [number, number];
        let colorScheme: { prop: string; gradient: $3Dmol.Gradient } | undefined;
        if (this._options.color.map.value === 'rwb') {
            // min and max are swapped to ensure red is used for high values, blue for low values
            colorScheme = { prop: this._options.color.property.value, gradient: new $3Dmol.Gradient.RWB(max, min) };
        } else if (this._options.color.map.value === 'sinebow') {
            colorScheme = { prop: this._options.color.property.value, gradient: new $3Dmol.Gradient.Sinebow(min, max) };
        }

        const style: Partial<$3Dmol.AtomStyleSpec> = {};

        if (this._options.atoms.value) {
            style.sphere = {
                scale: this._options.spaceFilling.value ? 1.0 : 0.22,
                colorscheme: colorScheme,
            };
        }
        if (this._options.bonds.value) {
            style.stick = {
                radius: 0.15,
                colorscheme: colorScheme,
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
                // slightly smaller radius than the main style
                radius: 0.149,
                opacity: defaultOpacity(),
                hidden: !this._options.bonds.value,
            };

            if (bgStyle === 'ball-stick' && this._options.shape.value === '') {
                style.sphere = {
                    // slightly smaller scale than the main style
                    scale: 0.219,
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
}
