/**
 * @packageDocumentation
 * @module settings
 */

import assert from 'assert';

import Collapse from '../collapse';
import Modal from '../modal';
import { Settings } from '../dataset';
import { HTMLOption, JSOption, OptionsGroup } from '../options';
import { optionValidator } from '../options';
import { PositioningCallback, Warnings, arrayMaxMin, getByID, makeDraggable } from '../utils';
import { CameraState, validateCamera } from '../utils/camera';
import { NumericProperties, NumericProperty } from './data';
import * as styles from '../styles';

import { COLOR_MAPS } from './colorscales';

import BARS_SVG from '../static/bars.svg';
import HTML_OPTIONS from './options.html.in';

// in 3D mode, only strings are supported for 'marker.symbol'.
/* eslint-disable */
const markers3d = require('./plotly/markers3d');
const POSSIBLE_SYMBOLS_IN_3D: string[] = Object.keys(markers3d.default);
/* eslint-enable */
export function get3DSymbol(i: number): string {
    return POSSIBLE_SYMBOLS_IN_3D[i % POSSIBLE_SYMBOLS_IN_3D.length];
}

/** HTML element holding settings for a given axis (x, y, z, color) */
export class AxisOptions extends OptionsGroup {
    public property: HTMLOption<'string'>;
    public scale: HTMLOption<'string'>;
    public min: HTMLOption<'number'>;
    public max: HTMLOption<'number'>;

    constructor(validProperties: string[], warnings: Warnings) {
        assert(validProperties.length > 0);
        super(warnings);
        this.max = new HTMLOption('number', NaN);
        this.min = new HTMLOption('number', NaN);
        this.property = new HTMLOption('string', validProperties[0]);
        this.scale = new HTMLOption('string', 'linear');

        this.property.validate = optionValidator(validProperties, 'axis');
        this.scale.validate = optionValidator(['linear', 'log'], 'axis scale');
    }

    /** Disable auxiliary settings (min/max/scale) related to this axis */
    public disable(): void {
        this.max.disable();
        this.min.disable();
        this.scale.disable();
    }

    /** Enable auxiliary settings (min/max/scale) related to this axis */
    public enable(): void {
        this.max.enable();
        this.min.enable();
        this.scale.enable();
    }
}

export class MapOptions extends OptionsGroup {
    public x: AxisOptions;
    public y: AxisOptions;
    public z: AxisOptions;
    public symbol: HTMLOption<'string'>;
    public color: {
        mode: HTMLOption<'string'>;
        property: HTMLOption<'string'>;
        min: HTMLOption<'number'>;
        max: HTMLOption<'number'>;
        palette: HTMLOption<'string'>;
        opacity: HTMLOption<'number'>;
        select: {
            mode: HTMLOption<'string'>;
            category: HTMLOption<'string'>;
            min: HTMLOption<'number'>;
            max: HTMLOption<'number'>;
        };
    };
    public size: {
        factor: HTMLOption<'number'>;
        mode: HTMLOption<'string'>;
        property: HTMLOption<'string'>;
    };
    public markerOutline: HTMLOption<'boolean'>;
    public joinPoints: HTMLOption<'boolean'>;
    public useLOD: HTMLOption<'boolean'>;
    public camera: JSOption<CameraState | undefined>;

    /// The HTML button to open the settings modal
    private _openModal: HTMLElement;
    /// The Modal instance
    private _modal: Modal;
    // Callback to get the initial positioning of the settings modal.
    private _positionSettingsModal: PositioningCallback;

    constructor(
        root: Node,
        properties: NumericProperties,
        positionSettings: PositioningCallback,
        settings: Settings = {},
        warnings: Warnings
    ) {
        super(warnings);

        // Setup axes
        const propertiesNames = Object.keys(properties);
        if (propertiesNames.length < 2) {
            throw new Error(
                'Cannot show a map because the dataset contains fewer than two properties.'
            );
        }
        this.x = new AxisOptions(propertiesNames, this.warnings);
        this.y = new AxisOptions(propertiesNames, this.warnings);
        // For z and color '' is a valid value
        this.z = new AxisOptions(propertiesNames.concat(['']), this.warnings);
        // Initialise symbol
        this.symbol = new HTMLOption('string', '');
        const validSymbols = [''];
        for (const key in properties) {
            if (properties[key].string !== undefined) {
                validSymbols.push(key);
            }
        }
        this.symbol.validate = optionValidator(validSymbols, 'symbol');

        // Initialise color
        this.color = {
            mode: new HTMLOption('string', 'linear'),
            property: new HTMLOption('string', ''),
            min: new HTMLOption('number', NaN),
            max: new HTMLOption('number', NaN),
            palette: new HTMLOption('string', 'inferno'),
            opacity: new HTMLOption('number', 100),
            select: {
                mode: new HTMLOption('string', 'all'),
                category: new HTMLOption('string', ''),
                min: new HTMLOption('number', NaN),
                max: new HTMLOption('number', NaN),
            },
        };
        this.color.property.validate = optionValidator(propertiesNames.concat(['']), 'color');
        this.color.mode.validate = optionValidator(['linear', 'log', 'sqrt', 'inverse'], 'mode');
        this.color.palette.validate = optionValidator(Object.keys(COLOR_MAPS), 'palette');
        this.color.opacity.validate = (value) => {
            if (value < 1 || value > 100) {
                throw Error(`opacity must be between 1 and 100, got ${value}`);
            }
        };
        this.color.select.mode.validate = optionValidator(
            ['all', 'range-gray', 'category-gray', 'range-hide', 'category-hide'],
            'selection mode'
        );
        // Initialise size
        this.size = {
            factor: new HTMLOption('number', 50),
            mode: new HTMLOption('string', 'linear'),
            property: new HTMLOption('string', ''),
        };
        this.size.property.validate = optionValidator(propertiesNames.concat(['']), 'size');
        this.size.factor.validate = (value) => {
            if (value < 1 || value > 100) {
                throw Error(`size factor must be between 0 and 100, got ${value}`);
            }
        };
        this.size.mode.validate = optionValidator(
            ['linear', 'log', 'sqrt', 'inverse', 'flip-linear', 'proportional'],
            'mode'
        );

        // outline and line
        this.markerOutline = new HTMLOption('boolean', true);
        this.joinPoints = new HTMLOption('boolean', false);
        // adaptive resolution for large datasets
        this.useLOD = new HTMLOption('boolean', true);
        this.camera = new JSOption<CameraState | undefined>(undefined);
        // custom equality check to avoid triggering unnecessary updates
        this.camera.equals = (a, b) => {
            if (a === undefined && b === undefined) {
                return true;
            }

            if (a === undefined || b === undefined) {
                return false;
            }

            const vecEquals = (
                v1: { x: number; y: number; z: number },
                v2: { x: number; y: number; z: number }
            ) => {
                return v1.x === v2.x && v1.y === v2.y && v1.z === v2.z;
            };

            return (
                vecEquals(a.eye, b.eye) &&
                vecEquals(a.center, b.center) &&
                vecEquals(a.up, b.up) &&
                a.zoom === b.zoom
            );
        };
        this.camera.validate = (value: CameraState | undefined) => {
            if (value !== undefined) {
                validateCamera(value);
                if (value.zoom <= 0) {
                    throw Error(`zoom factor must be greater than zero, got ${value.zoom}`);
                }
            }
        };

        // Setup default values
        this.x.property.value = propertiesNames[0];
        this.y.property.value = propertiesNames[1];
        this.z.property.value = '';
        if (propertiesNames.length > 2) {
            this.color.property.value = propertiesNames[2];
        } else {
            this.color.property.value = '';
        }

        // Setup modal
        this._positionSettingsModal = positionSettings;
        const { openModal, modal } = this._createSettingsHTML();
        this._modal = modal;
        this._openModal = openModal;
        root.appendChild(this._openModal);

        // Attach callbacks
        this._bind(properties);

        this.color.select.mode.onchange.push(() => {
            if (this.color.select.mode.value.startsWith('range')) {
                if (this.color.select.min.value > this.color.select.max.value) {
                    this.warnings.sendMessage(
                        'The selection min value is greater than the max value! Resetting the range.'
                    );
                    this.color.select.min.reset();
                    this.color.select.max.reset();
                }
            }
        });

        const validateSelectRange = (minOrMax: 'min' | 'max') => {
            return () => {
                if (this.color.select.mode.value.startsWith('range')) {
                    const min = this.color.select.min.value;
                    const max = this.color.select.max.value;
                    if (!isNaN(min) && !isNaN(max) && min > max) {
                        this.warnings.sendMessage(
                            `The selection ${minOrMax} value makes min > max! Resetting it.`
                        );
                        if (minOrMax === 'min') {
                            this.color.select.min.reset();
                        } else {
                            this.color.select.max.reset();
                        }
                    }
                }
            };
        };
        this.color.select.min.onchange.push(validateSelectRange('min'));
        this.color.select.max.onchange.push(validateSelectRange('max'));

        this.color.property.onchange.push(() => {
            // reset selection range when the property changes
            this.color.select.min.value = NaN;
            this.color.select.max.value = NaN;

            // disable range selection if color is fixed
            const selectMode = this.getModalElement<HTMLSelectElement>('map-color-select-mode');
            const rangeOptions = selectMode.querySelectorAll('option[value^="range"]');
            if (this.color.property.value === '') {
                if (this.color.select.mode.value.startsWith('range')) {
                    this.color.select.mode.value = 'all';
                }
                rangeOptions.forEach((opt) => ((opt as HTMLOptionElement).disabled = true));
            } else {
                rangeOptions.forEach((opt) => ((opt as HTMLOptionElement).disabled = false));
            }
        });

        // Apply new settings to the modal options
        this.applySettings(settings);
    }

    /** Get in a element in the modal from its id */
    // eslint-disable-next-line @typescript-eslint/no-unnecessary-type-parameters
    public getModalElement<T extends HTMLElement = HTMLElement>(id: string): T {
        return getByID(id, this._modal.shadow);
    }

    /**
     * Apply saved settings to all the map options
     *
     * @param settings settings for all panels
     */
    public applySettings(settings: Settings): void {
        // deal with backward compatibility: size.mode === 'constant' should be
        // the same as `size.property === ''`
        if ('size' in settings) {
            const size = settings.size as Settings;
            if ('mode' in size && size.mode === 'constant') {
                delete size.mode;
                size.property = '';
            }
        }
        if ('color' in settings) {
            const color = settings.color as Settings;
            if ('scale' in color) {
                color.mode = color.scale;
                delete color.scale;
            }
            if ('select' in color) {
                const select = color.select as Settings;
                if (select.mode === 'range') {
                    select.mode = 'range-gray';
                } else if (select.mode === 'category') {
                    select.mode = 'category-gray';
                }
            }
        }
        if ('palette' in settings) {
            if ('color' in settings) {
                const color = settings.color as Settings;
                color.palette = settings.palette;
            } else {
                settings.color = { palette: settings.palete };
            }
            delete settings.palette;
        }
        super.applySettings(settings);
    }

    /**
     * Remove all HTML added by this {@link MapSettings} in the current document
     */
    public remove(): void {
        this._modal.remove();
        this._openModal.remove();
    }

    /** Is the current plot in 3D mode? */
    public is3D(): boolean {
        return this.z.property.value !== '';
    }

    /** Does the current plot use color values? */
    public hasColors(): boolean {
        return this.color.property.value !== '';
    }

    /** Get the plotly hovertemplate depending on `this._current.color` */
    public hovertemplate(): string {
        if (this.hasColors()) {
            let property = this.color.property.value;
            switch (this.color.mode.value) {
                case 'inverse':
                    property = `(${property})<sup>-1</sup>`;
                    break;
                case 'log':
                    property = `log(${property})`;
                    break;
                case 'sqrt':
                    property = `sqrt(${property})`;
                    break;
                case 'linear':
                    break;
                default:
                    break;
            }

            return property + ': %{customdata:.2f}<extra></extra>';
        } else {
            return '%{x:.2f}, %{y:.2f}<extra></extra>';
        }
    }

    /**
     * Get the values to use as colors according to the user-selected property
     * and scale
     */
    public calculateColors(rawColors: number[]): Array<number | string> {
        let scaleMode = this.color.mode.value;
        const { min, max } = arrayMaxMin(rawColors);
        if (max === min) {
            scaleMode = 'fixed';
        }

        const values = rawColors.map((v: number) => {
            let transformed = 0.5; // default
            switch (scaleMode) {
                case 'inverse':
                    transformed = 1.0 / v;
                    break;
                case 'log':
                    transformed = Math.log10(v);
                    break;
                case 'sqrt':
                    transformed = Math.sqrt(v);
                    break;
                case 'linear':
                    transformed = 1.0 * v;
                    break;
                default:
                    // corresponds to 'fixed'
                    transformed = 0.5;
                    break;
            }

            return isNaN(transformed) ? '#aaaaaa' : transformed;
        });

        return values;
    }

    /**
     * Get the values to use as marker size with the given plotly `trace`, or
     * all of them if `trace === undefined`.
     */
    public calculateSizes(rawSizes: number[]): number[] {
        const logSlider = (value: number) => {
            const min_slider = 1;
            const max_slider = 100;

            // go from 1/6th of the size to 6 time the size
            const min_value = Math.log(1.0 / 6.0);
            const max_value = Math.log(6.0);

            const tmp = (max_value - min_value) / (max_slider - min_slider);
            return Math.exp(min_value + tmp * (value - min_slider));
        };

        const userFactor = logSlider(this.size.factor.value);

        let scaleMode = this.size.mode.value;
        const { min, max } = arrayMaxMin(rawSizes);
        const defaultSize = this.is3D() ? 800 : 300;
        const bottomLimit = 0.1; // lower limit to prevent size of 0
        const defaultScaled = 0.3;
        const nonzeromin = min > 0 ? min : 1e-6 * (max - min); // non-zero minimum value for scales needing it
        const values = rawSizes.map((v: number) => {
            // normalize between 0 and 1, then scale by the user provided value
            let scaled = defaultScaled; // default
            if (max === min) {
                scaleMode = 'fixed';
            } else {
                scaled = (v - min) / (max - min);
            }
            switch (scaleMode) {
                case 'proportional':
                    // absolude proportionality - zero to max
                    // nb this will break for negative values!
                    scaled = v / Math.abs(max);
                    break;
                case 'inverse':
                    // inverse mapping
                    // nb this will break for negative values!
                    scaled = nonzeromin / v;
                    break;
                case 'log':
                    // log scale magnitude
                    // nb this will break for negative values!
                    scaled = Math.log(v / nonzeromin) / Math.log(max / nonzeromin);
                    break;
                case 'sqrt':
                    // sqrt mapping
                    // nb this will break for negative values!
                    scaled = Math.sqrt(v / Math.abs(max));
                    break;
                case 'linear':
                    scaled = 1.0 * scaled;
                    break;
                case 'flip-linear':
                    scaled = 1.0 - scaled;
                    break;
                default:
                    // corresponds to 'constant'
                    scaled = defaultScaled - bottomLimit;
                    break;
            }
            scaled = scaled + bottomLimit; // minimum size is enforced
            // nb: we use scalemode=area so the property value is linked to the
            // _area_ of the points (which is the perceptually correct thing to do)
            return defaultSize * scaled * userFactor;
        });

        const someValuesNaN = values.some((value) => isNaN(value) || value < 0);

        if (someValuesNaN) {
            this.warnings.sendMessage(
                `After applying the selected scaling mode ${scaleMode}, some point sizes ` +
                    `evaluated to invalid values. These points will be displayed at the minimum size.`
            );
            return values.map((v: number) => {
                return isNaN(v) ? bottomLimit * defaultSize * userFactor : v;
            });
        } else {
            return values;
        }
    }

    /** Given the property, return the symbols */
    public getSymbols(property: NumericProperty): number[] | string[] {
        /** How many different symbols are being displayed */
        assert(property.string !== undefined);
        const symbolsCount = property.string.strings().length;

        if (this.is3D()) {
            // If we need more symbols than available, we'll send a warning
            // and repeat existing ones
            if (symbolsCount > POSSIBLE_SYMBOLS_IN_3D.length) {
                this.warnings.sendMessage(
                    `${symbolsCount} symbols are required, but we only have ${POSSIBLE_SYMBOLS_IN_3D.length}. Some symbols will be repeated`
                );
            }
            const symbols = property.values.map(get3DSymbol);
            return symbols;
        } else {
            return property.values;
        }
    }

    /**
     * Create the settings modal by adding HTML to the page
     * @return the modal HTML element, not yet inserted in the document
     */
    private _createSettingsHTML(): { openModal: HTMLElement; modal: Modal } {
        const template = document.createElement('template');
        template.innerHTML = `<button
            class="btn btn-light btn-sm chsp-viewer-button"
            style="top: 4px; left: 5px; opacity: 1;">
                <div>${BARS_SVG}</div>
            </button>`;
        const openModal = template.content.firstChild as HTMLElement;

        // replace id to ensure they are unique even if we have multiple viewers
        // on a single page
        // prettier-ignore
        template.innerHTML = HTML_OPTIONS;

        const modalElement = template.content.querySelector('.modal');
        assert(modalElement !== null && modalElement instanceof HTMLElement);
        const modalDialog = modalElement.querySelector('.modal-dialog');
        assert(modalDialog !== null && modalDialog instanceof HTMLElement);

        // make the settings modal draggable
        makeDraggable(modalDialog, '.modal-header');

        const modal = new Modal(modalElement);
        modal.shadow.adoptedStyleSheets = [styles.bootstrap, styles.chemiscope];

        Collapse.initialize(modalElement);

        // Position modal according to this._positionSettingsModal
        openModal.addEventListener('click', () => {
            // only set style once, on first open, and keep previous position
            // on next open to keep the 'dragged-to' position
            if (modalDialog.getAttribute('data-initial-modal-positions-set') === null) {
                modalDialog.setAttribute('data-initial-modal-positions-set', 'true');

                // display: block to ensure modalDialog.offsetWidth is non-zero
                (modalDialog.parentNode as HTMLElement).style.display = 'block';

                // set width first, since setting position can influence it
                // scale width of the larger modal-lg class
                modalDialog.style.width = `${modalDialog.offsetWidth / 1.5}px`;
                // minimum width so that text in rows remains on a single line
                modalDialog.style.minWidth = `500px`;
                // unset margins when using position: fixed
                modalDialog.style.margin = '0';
                modalDialog.style.position = 'fixed';

                const { top, left } = this._positionSettingsModal(
                    modalDialog.getBoundingClientRect()
                );

                modalDialog.style.top = `${top}px`;
                modalDialog.style.left = `${left}px`;
            }

            modal.open();
        });

        // Stop propagation of keydown events. This is required for the Jupyter integration,
        // otherwise jupyter tries to interpret key press in the modal as its own input
        modalElement.addEventListener('keydown', (event) => {
            event.stopPropagation();
        });

        return { openModal, modal };
    }

    /** Bind all options to the corresponding HTML elements */
    private _bind(properties: NumericProperties): void {
        // ======= data used as x values
        const selectXProperty = this.getModalElement<HTMLSelectElement>('map-x-property');
        selectXProperty.options.length = 0;
        for (const key in properties) {
            selectXProperty.options.add(new Option(key, key));
        }
        this.x.property.bind(selectXProperty, 'value');
        this.x.min.bind(this.getModalElement('map-x-min'), 'value');
        this.x.max.bind(this.getModalElement('map-x-max'), 'value');
        this.x.scale.bind(this.getModalElement('map-x-scale'), 'value');

        // ======= data used as y values
        const selectYProperty = this.getModalElement<HTMLSelectElement>('map-y-property');

        selectYProperty.options.length = 0;
        for (const key in properties) {
            selectYProperty.options.add(new Option(key, key));
        }
        this.y.property.bind(selectYProperty, 'value');
        this.y.min.bind(this.getModalElement('map-y-min'), 'value');
        this.y.max.bind(this.getModalElement('map-y-max'), 'value');
        this.y.scale.bind(this.getModalElement('map-y-scale'), 'value');

        // ======= data used as z values
        const selectZProperty = this.getModalElement<HTMLSelectElement>('map-z-property');
        // first option is 'none'
        selectZProperty.options.length = 0;
        selectZProperty.options.add(new Option('none', ''));
        for (const key in properties) {
            selectZProperty.options.add(new Option(key, key));
        }
        this.z.property.bind(selectZProperty, 'value');
        this.z.min.bind(this.getModalElement('map-z-min'), 'value');
        this.z.max.bind(this.getModalElement('map-z-max'), 'value');
        this.z.scale.bind(this.getModalElement('map-z-scale'), 'value');

        // ======= data used as color values
        const selectColorProperty = this.getModalElement<HTMLSelectElement>('map-color-property');
        // first option is 'fixed'
        selectColorProperty.options.length = 0;
        selectColorProperty.options.add(new Option('fixed', ''));
        for (const key in properties) {
            selectColorProperty.options.add(new Option(key, key));
        }
        this.color.property.bind(selectColorProperty, 'value');
        this.color.mode.bind(this.getModalElement('map-color-transform'), 'value');
        this.color.min.bind(this.getModalElement('map-color-min'), 'value');
        this.color.max.bind(this.getModalElement('map-color-max'), 'value');
        this.color.opacity.bind(this.getModalElement('map-opacity'), 'value');

        // ======= color palette
        const selectPalette = this.getModalElement<HTMLSelectElement>('map-color-palette');
        selectPalette.length = 0;
        for (const key in COLOR_MAPS) {
            selectPalette.options.add(new Option(key, key));
        }
        this.color.palette.bind(selectPalette, 'value');

        // ======= color select
        const selectSelectMode = this.getModalElement<HTMLSelectElement>('map-color-select-mode');
        this.color.select.mode.bind(selectSelectMode, 'value');

        const selectSelectCategory = this.getModalElement<HTMLSelectElement>(
            'map-color-select-category'
        );
        selectSelectCategory.length = 0;
        let hasCategorical = false;
        for (const key in properties) {
            const prop = properties[key];
            if (prop.string !== undefined) {
                hasCategorical = true;
                const values = prop.string.strings();
                for (const val of values) {
                    const optionValue = `${key}/${val}`;
                    selectSelectCategory.add(new Option(optionValue, optionValue));
                }
            }
        }

        // disable category selection if there are no categorical properties
        const categoryOptions = selectSelectMode.querySelectorAll('option[value^="category"]');
        if (!hasCategorical) {
            categoryOptions.forEach((opt) => ((opt as HTMLOptionElement).disabled = true));
        }

        // disable range selection if color is fixed
        const rangeOptions = selectSelectMode.querySelectorAll('option[value^="range"]');
        if (this.color.property.value === '') {
            rangeOptions.forEach((opt) => ((opt as HTMLOptionElement).disabled = true));
        }

        if (selectSelectCategory.options.length > 0) {
            this.color.select.category.value = selectSelectCategory.options[0].value;
        }
        this.color.select.category.bind(selectSelectCategory, 'value');

        this.color.select.min.bind(this.getModalElement('map-color-select-min'), 'value');
        this.color.select.max.bind(this.getModalElement('map-color-select-max'), 'value');

        const updateSelectVisibility = () => {
            const mode = this.color.select.mode.value;
            const containerMode = this.getModalElement('map-color-select-container');
            const containerCategory = this.getModalElement('map-color-select-category-container');
            const containerMin = this.getModalElement('map-color-select-min-container');
            const containerMax = this.getModalElement('map-color-select-max-container');

            if (mode.startsWith('range')) {
                containerMode.style.gridColumn = 'auto / span 2';
                containerCategory.style.display = 'none';
                containerMin.style.display = 'flex';
                containerMin.style.gridColumn = 'span 2';
                containerMax.style.display = 'flex';
                containerMax.style.gridColumn = 'span 2';
            } else if (mode.startsWith('category')) {
                containerMode.style.gridColumn = 'span 3';
                containerCategory.style.display = 'flex';
                containerCategory.style.gridColumn = 'span 3';
                containerMin.style.display = 'none';
                containerMax.style.display = 'none';
            } else {
                containerMode.style.gridColumn = 'span 3';
                containerCategory.style.display = 'none';
                containerMin.style.display = 'none';
                containerMax.style.display = 'none';
            }
        };
        this.color.select.mode.onchange.push(updateSelectVisibility);
        // Call it once to set initial state (e.g. if settings loaded 'range')
        setTimeout(updateSelectVisibility, 0);

        // ======= marker symbols
        const selectSymbolProperty = this.getModalElement<HTMLSelectElement>('map-symbol-property');
        // first option is 'fixed'
        selectSymbolProperty.options.length = 0;
        selectSymbolProperty.options.add(new Option('fixed', ''));
        for (const key in properties) {
            if (properties[key].string !== undefined) {
                selectSymbolProperty.options.add(new Option(key, key));
            }
        }
        this.symbol.bind(selectSymbolProperty, 'value');

        // ======= marker size
        const selectSizeProperty = this.getModalElement<HTMLSelectElement>('map-size-property');
        // first option is 'fixed'
        selectSizeProperty.options.length = 0;
        selectSizeProperty.options.add(new Option('fixed', ''));
        for (const key in properties) {
            selectSizeProperty.options.add(new Option(key, key));
        }
        this.size.property.bind(selectSizeProperty, 'value');
        this.size.factor.bind(this.getModalElement('map-size-factor'), 'value');
        this.size.mode.bind(this.getModalElement('map-size-transform'), 'value');
        // ====== marker outline and line trace
        this.markerOutline.bind(this.getModalElement('map-marker-outline'), 'checked');
        this.joinPoints.bind(this.getModalElement('map-join-points'), 'checked');
        // ====== use LOD for large datasets
        this.useLOD.bind(this.getModalElement('map-lod'), 'checked');
    }

    /** Get the colorscale to use for markers in the main plotly trace */
    public colorScale(): Plotly.ColorScale {
        return COLOR_MAPS[this.color.palette.value];
    }

    /** Changes the min/max range label between linear and log appropriately */
    public setLogLabel(axis: AxisOptions, axisName: string): void {
        const minInputLabel = this.getModalElement(`map-${axisName}-min-label`);
        const maxInputLabel = this.getModalElement(`map-${axisName}-max-label`);

        if (axis.scale.value === 'log') {
            minInputLabel.innerHTML = 'min: 10^';
            maxInputLabel.innerHTML = 'max: 10^';
        } else {
            minInputLabel.innerHTML = 'min:';
            maxInputLabel.innerHTML = 'max:';
        }
    }

    /** Show or hide the LOD toggle based on dataset size */
    public showLODOption(show: boolean): void {
        const el = this.getModalElement('map-lod-div');
        el.style.display = show ? 'flex' : 'none';
    }

    public getCategorySelectionProperty(): string | null {
        if (!this.color.select.mode.value.startsWith('category')) return null;
        const val = this.color.select.category.value;
        const sep = val.indexOf('/');
        if (sep === -1) return null;
        return val.substring(0, sep);
    }

    /**
     * Get a mask of selected points based on the current selection mode
     *
     * @param values values of the property used for coloring (for range mode)
     * @param categoryValues values of the property used for categorization (property/value string)
     */
    public getFilterMask(values: number[], categoryValues?: string[]): boolean[] {
        const mode = this.color.select.mode.value;
        if (mode === 'all') {
            return new Array(values.length).fill(true) as boolean[];
        }

        if (mode.startsWith('range')) {
            const min = this.color.select.min.value;
            const max = this.color.select.max.value;
            // if min/max are NaN, they are ignored (open range)
            return values.map((v) => {
                if (!isNaN(min) && v < min) return false;
                if (!isNaN(max) && v > max) return false;
                return true;
            });
        }

        if (mode.startsWith('category')) {
            const target = this.color.select.category.value; // "property/value"
            if (!categoryValues) {
                return new Array(values.length).fill(true) as boolean[];
            }

            // let's split target
            const sepIndex = target.indexOf('/');
            if (sepIndex === -1) return new Array(values.length).fill(true) as boolean[];
            const targetVal = target.substring(sepIndex + 1);

            return categoryValues.map((v) => v === targetVal);
        }

        return new Array(values.length).fill(true) as boolean[];
    }

    /**
     * Convert a numeric value to a color string based on the current palette
     */
    public valueToColor(value: number, min: number, max: number): string {
        const palette = COLOR_MAPS[this.color.palette.value];
        // palette is [stop, colorString][]
        // normalize value
        let t = (value - min) / (max - min);
        if (isNaN(t)) t = 0.5;
        if (t < 0) t = 0;
        if (t > 1) t = 1;

        // find segment
        // palette is sorted by stop 0..1
        let i = 0;
        while (i < palette.length - 1 && palette[i + 1][0] < t) {
            i++;
        }
        // segment is i to i+1
        const start = palette[i];
        const end = palette[i + 1] || start; // handle t=1 or end

        const t0 = start[0];
        const t1 = end[0];
        const c0 = start[1]; // "rgb(r, g, b)"
        const c1 = end[1];

        if (t1 === t0) return c0;

        const f = (t - t0) / (t1 - t0);

        const parseRgb = (s: string) => {
            const m = s.match(/rgb\((\d+),\s*(\d+),\s*(\d+)\)/);
            if (m) return [parseInt(m[1], 10), parseInt(m[2], 10), parseInt(m[3], 10)];
            return [0, 0, 0];
        };

        const rgb0 = parseRgb(c0);
        const rgb1 = parseRgb(c1);

        const r = Math.round(rgb0[0] + f * (rgb1[0] - rgb0[0]));
        const g = Math.round(rgb0[1] + f * (rgb1[1] - rgb0[1]));
        const b = Math.round(rgb0[2] + f * (rgb1[2] - rgb0[2]));

        return `rgb(${r}, ${g}, ${b})`;
    }
}
