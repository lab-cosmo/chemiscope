/**
 * @packageDocumentation
 * @module settings
 */

import assert from 'assert';

import Collapse from '../collapse';
import Modal from '../modal';
import { Settings } from '../dataset';
import { HTMLOption, OptionsGroup } from '../options';
import { optionValidator } from '../options';
import { PositioningCallback, Warnings, arrayMaxMin, getByID, makeDraggable } from '../utils';
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
    };
    public size: {
        factor: HTMLOption<'number'>;
        mode: HTMLOption<'string'>;
        property: HTMLOption<'string'>;
    };
    public markerOutline: HTMLOption<'boolean'>;
    public joinPoints: HTMLOption<'boolean'>;

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
        };
        this.color.property.validate = optionValidator(propertiesNames.concat(['']), 'color');
        this.color.mode.validate = optionValidator(['linear', 'log', 'sqrt', 'inverse'], 'mode');
        this.color.palette.validate = optionValidator(Object.keys(COLOR_MAPS), 'palette');
        this.color.opacity.validate = (value) => {
            if (value < 1 || value > 100) {
                throw Error(`opacity must be between 1 and 100, got ${value}`);
            }
        };
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

            return property + ': %{marker.color:.2f}<extra></extra>';
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
                `After applying the selected scaling mode ${scaleMode}, some point sizes` +
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

        // ======= color palette
        const selectPalette = this.getModalElement<HTMLSelectElement>('map-color-palette');
        selectPalette.length = 0;
        for (const key in COLOR_MAPS) {
            selectPalette.options.add(new Option(key, key));
        }
        this.color.palette.bind(selectPalette, 'value');
        this.color.opacity.bind(this.getModalElement('map-opacity'), 'value');

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
}
