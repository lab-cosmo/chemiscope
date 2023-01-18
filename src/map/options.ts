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
import { PositioningCallback } from '../utils';
import { arrayMaxMin, getByID, makeDraggable, sendWarning } from '../utils';
import { NumericProperties, NumericProperty } from './data';
import * as styles from '../styles';

import { COLOR_MAPS } from './colorscales';

import BARS_SVG from '../static/bars.svg';
import HTML_OPTIONS from './options.html.in';

// in 3D mode, only strings are supported for 'marker.symbol', and only very few
// of them. See https://github.com/plotly/plotly.js/issues/4205 as the plotly
// issue tracking more symbols in 3D mode.
const POSSIBLE_SYMBOLS_IN_3D = ['circle', 'square', 'diamond', 'cross', 'x'];

export function get3DSymbol(i: number): string {
    return POSSIBLE_SYMBOLS_IN_3D[i % POSSIBLE_SYMBOLS_IN_3D.length];
}

/** HTML element holding settings for a given axis (x, y, z, color) */
export class AxisOptions extends OptionsGroup {
    public property: HTMLOption<'string'>;
    public scale: HTMLOption<'string'>;
    public min: HTMLOption<'number'>;
    public max: HTMLOption<'number'>;

    constructor(validProperties: string[]) {
        assert(validProperties.length > 0);
        super();
        this.max = new HTMLOption('number', 0);
        this.min = new HTMLOption('number', 0);
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
    public color: AxisOptions;
    public palette: HTMLOption<'string'>;
    public symbol: HTMLOption<'string'>;
    public size: {
        factor: HTMLOption<'number'>;
        mode: HTMLOption<'string'>;
        property: HTMLOption<'string'>;
        reverse: HTMLOption<'boolean'>;
    };

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
        settings: Settings = {}
    ) {
        super();
        const propertiesName = Object.keys(properties);
        assert(propertiesName.length >= 2);

        this.x = new AxisOptions(propertiesName);
        this.y = new AxisOptions(propertiesName);
        // For z and color '' is a valid value
        this.z = new AxisOptions(propertiesName.concat(['']));
        this.color = new AxisOptions(propertiesName.concat(['']));

        this.symbol = new HTMLOption('string', '');
        const validSymbols = [''];
        for (const key in properties) {
            if (properties[key].string !== undefined) {
                validSymbols.push(key);
            }
        }
        this.symbol.validate = optionValidator(validSymbols, 'symbol');

        this.palette = new HTMLOption('string', 'inferno');
        this.palette.validate = optionValidator(Object.keys(COLOR_MAPS), 'palette');

        this.size = {
            factor: new HTMLOption('number', 50),
            mode: new HTMLOption('string', 'linear'),
            property: new HTMLOption('string', ''),
            reverse: new HTMLOption('boolean', false),
        };

        this.size.property.validate = optionValidator(propertiesName.concat(['']), 'size');
        this.size.factor.validate = (value) => {
            if (value < 1 || value > 100) {
                throw Error(`size factor must be between 0 and 100, got ${value}`);
            }
        };
        this.size.mode.validate = optionValidator(['linear', 'log', 'sqrt', 'inverse'], 'mode');

        this.x.property.value = propertiesName[0];
        this.y.property.value = propertiesName[1];
        this.z.property.value = '';

        if (propertiesName.length > 2) {
            this.color.property.value = propertiesName[2];
        } else {
            this.color.property.value = '';
        }

        this._positionSettingsModal = positionSettings;

        const { openModal, modal } = this._createSettingsHTML();
        this._modal = modal;
        this._openModal = openModal;
        root.appendChild(this._openModal);

        this._bind(properties);
        this.applySettings(settings);
    }

    /** Get in a element in the modal from its id */
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
            return this.color.property.value + ': %{marker.color:.2f}<extra></extra>';
        } else {
            return '%{x:.2f}, %{y:.2f}<extra></extra>';
        }
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
        const reversed = this.size.reverse.value;
        const { min, max } = arrayMaxMin(rawSizes);
        const defaultSize = this.is3D() ? 800 : 300;
        const bottomLimit = 0.1; // lower limit to prevent size of 0

        const values = rawSizes.map((v: number) => {
            // normalize between 0 and 1, then scale by the user provided value
            let scaled = 0.55; // default
            if (max === min) {
                scaleMode = 'fixed';
            } else {
                scaled = (v + bottomLimit - min) / (max - min);
                if (reversed) {
                    scaled = 1.0 + bottomLimit - scaled;
                }
            }
            switch (scaleMode) {
                case 'inverse':
                    scaled = 1.0 / scaled;
                    break;
                case 'log':
                    scaled = Math.log(scaled);
                    break;
                case 'sqrt':
                    scaled = Math.sqrt(scaled);
                    break;
                case 'linear':
                    scaled = 1.0 * scaled;
                    break;
                default:
                    // corresponds to 'constant'
                    scaled = 0.55;
                    break;
            }
            // since we are using scalemode: 'area', square the scaled value
            return defaultSize * scaled * scaled * userFactor * userFactor;
        });
        return values;
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
                sendWarning(
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
     * @param  guid unique identifier of this map, used as prefix for all
     *              elements ID
     * @return      the modal HTML element, not yet inserted in the document
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
                modalDialog.style.minWidth = `400px`;
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
        this.color.min.bind(this.getModalElement('map-color-min'), 'value');
        this.color.max.bind(this.getModalElement('map-color-max'), 'value');

        // ======= color palette
        const selectPalette = this.getModalElement<HTMLSelectElement>('map-color-palette');
        selectPalette.length = 0;
        for (const key in COLOR_MAPS) {
            selectPalette.options.add(new Option(key, key));
        }
        this.palette.bind(selectPalette, 'value');

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
        this.size.mode.bind(this.getModalElement('map-size-scale'), 'value');
        this.size.reverse.bind(this.getModalElement('map-size-reverse'), 'checked');
    }

    /** Get the colorscale to use for markers in the main plotly trace */
    public colorScale(): Plotly.ColorScale {
        return COLOR_MAPS[this.palette.value];
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
