/**
 * @packageDocumentation
 * @module settings
 */

import assert from 'assert';

import { HTMLOption, OptionsGroup, SavedSettings } from '../options';
import { optionValidator } from '../options';
import { PositioningCallback, arrayMaxMin, getByID, makeDraggable, sendWarning } from '../utils';
import { NumericProperties, NumericProperty } from './data';

import { COLOR_MAPS } from './colorscales';

import BARS_SVG from '../static/bars.svg';
import HTML_OPTIONS from './options.html';

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

    /// The HTML element containing the settings modal
    private _modal: HTMLElement;
    // Callback to get the initial positioning of the settings modal.
    private _positionSettingsModal: PositioningCallback;

    constructor(
        root: HTMLElement,
        properties: NumericProperties,
        positionSettings: PositioningCallback,
        settings: SavedSettings = {}
    ) {
        super();
        const propertiesName = Object.keys(properties);
        if (propertiesName.length < 2) {
            throw Error('we need at least two properties to plot in the map');
        }

        this.x = new AxisOptions(propertiesName);
        this.y = new AxisOptions(propertiesName);
        // For z and color, '' is a valid value
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
            mode: new HTMLOption('string', 'constant'),
            property: new HTMLOption('string', propertiesName[0]),
            reverse: new HTMLOption('boolean', false),
        };
        this.size.property.validate = optionValidator(propertiesName, 'size');
        this.size.factor.validate = (value) => {
            if (value < 1 || value > 100) {
                throw Error(`size factor must be between 0 and 100, got ${value}`);
            }
        };
        this.size.mode.validate = optionValidator(
            ['constant', 'linear', 'log', 'sqrt', 'inverse'],
            'size'
        );

        this.x.property.value = propertiesName[0];
        this.y.property.value = propertiesName[1];
        this.z.property.value = '';

        if (propertiesName.length >= 2) {
            this.color.property.value = propertiesName[2];
        } else {
            this.color.property.value = '';
        }

        this._positionSettingsModal = positionSettings;
        this._modal = this._insertSettingsHTML(root);
        document.body.appendChild(this._modal);

        this._bind(properties);
        this.applySettings(settings);
    }

    /**
     * Apply saved settings to all the map options
     *
     * @param settings settings for all panels
     */
    public applySettings(settings: SavedSettings): void {
        // deal with backward compatibility: size.property === '' used to mean
        // "use a constant scaling"
        if ('size' in settings) {
            const size = settings.size as SavedSettings;
            if ('property' in size && size.property === '') {
                delete size.property;
                size.mode = 'constant';
            }
        }
        super.applySettings(settings);
    }

    /**
     * Remove all HTML added by this [[MapSettings]] in the current document
     */
    public remove(): void {
        if (this._modal.classList.contains('show')) {
            const close = this._modal.querySelector('.close');
            assert(close !== null);
            (close as HTMLElement).click();
        }
        this._modal.remove();
    }

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

        let values;
        if (this.size.mode.value !== 'constant') {
            const scaleMode = this.size.mode.value;
            const reversed = this.size.reverse.value;
            const { min, max } = arrayMaxMin(rawSizes);
            const defaultSize = this.is3D() ? 2000 : 150;
            const bottomLimit = 0.1; // lower limit to prevent size of 0

            values = rawSizes.map((v: number) => {
                // normalize between 0 and 1, then scale by the user provided value
                let scaled = (v + bottomLimit - min) / (max - min);
                if (reversed) {
                    scaled = 1.0 + bottomLimit - scaled;
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
                        scaled = 1.0;
                        break;
                }
                // since we are using scalemode: 'area', square the scaled value
                return defaultSize * scaled * scaled * userFactor * userFactor;
            });
        } else {
            // we need to use an array instead of a single value because of
            // https://github.com/plotly/plotly.js/issues/2735
            const defaultSize = this.is3D() ? 500 : 50;
            values = Array(rawSizes.length).fill(defaultSize * userFactor) as number[];
        }
        return values;
    }

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
     * @param  root root element in which the 'open settings' button will be placed
     * @return      the modal HTML element, not yet inserted in the document
     */
    private _insertSettingsHTML(root: HTMLElement): HTMLElement {
        const template = document.createElement('template');
        template.innerHTML = `<button
            class="btn btn-light btn-sm chsp-viewer-button"
            data-target='#chsp-settings'
            data-toggle="modal"
            style="top: 4px; left: 5px; opacity: 1;">
                <div>${BARS_SVG}</div>
            </button>`;
        const openSettings = template.content.firstChild as HTMLElement;
        root.append(openSettings);

        // TODO: set unique HTML id in the settings to allow multiple map in
        // the same page
        template.innerHTML = HTML_OPTIONS;
        const modal = template.content.firstChild as HTMLElement;

        const modalDialog = modal.childNodes[1] as HTMLElement;
        assert(modalDialog !== undefined);
        assert(modalDialog.classList.contains('modal-dialog'));
        // make the settings modal draggable
        makeDraggable(modalDialog, '.modal-header');

        // Position modal according to this._positionSettingsModal
        openSettings.addEventListener('click', () => {
            // only set style once, on first open, and keep previous position
            // on next open to keep the 'draged-to' position
            if (modalDialog.getAttribute('data-initial-modal-positions-set') === null) {
                modalDialog.setAttribute('data-initial-modal-positions-set', 'true');

                // display: block to ensure modalDialog.offsetWidth is non-zero
                (modalDialog.parentNode as HTMLElement).style.display = 'block';

                const { top, left } = this._positionSettingsModal(
                    modalDialog.getBoundingClientRect()
                );

                // set width first, since setting position can influence it
                modalDialog.style.width = `${modalDialog.offsetWidth}px`;
                // unset margins when using position: fixed
                modalDialog.style.margin = '0';
                modalDialog.style.position = 'fixed';
                modalDialog.style.top = `${top}px`;
                modalDialog.style.left = `${left}px`;
            }
        });

        return modal;
    }

    /** Bind all options to the corresponding HTML elements */
    private _bind(properties: NumericProperties): void {
        // ======= data used as x values
        const selectXProperty = getByID<HTMLSelectElement>('chsp-x');
        selectXProperty.options.length = 0;
        for (const key in properties) {
            selectXProperty.options.add(new Option(key, key));
        }
        this.x.property.bind(selectXProperty, 'value');
        this.x.min.bind('chsp-x-min', 'value');
        this.x.max.bind('chsp-x-max', 'value');
        this.x.scale.bind('chsp-x-scale', 'value');

        // ======= data used as y values
        const selectYProperty = getByID<HTMLSelectElement>('chsp-y');
        selectYProperty.options.length = 0;
        for (const key in properties) {
            selectYProperty.options.add(new Option(key, key));
        }
        this.y.property.bind(selectYProperty, 'value');
        this.y.min.bind('chsp-y-min', 'value');
        this.y.max.bind('chsp-y-max', 'value');
        this.y.scale.bind('chsp-y-scale', 'value');

        // ======= data used as z values
        const selectZProperty = getByID<HTMLSelectElement>('chsp-z');
        // first option is 'none'
        selectZProperty.options.length = 0;
        selectZProperty.options.add(new Option('none', ''));
        for (const key in properties) {
            selectZProperty.options.add(new Option(key, key));
        }
        this.z.property.bind(selectZProperty, 'value');
        this.z.min.bind('chsp-z-min', 'value');
        this.z.max.bind('chsp-z-max', 'value');
        this.z.scale.bind('chsp-z-scale', 'value');

        // ======= data used as color values
        const selectColorProperty = getByID<HTMLSelectElement>('chsp-color');
        // first option is 'none'
        selectColorProperty.options.length = 0;
        selectColorProperty.options.add(new Option('none', ''));
        for (const key in properties) {
            selectColorProperty.options.add(new Option(key, key));
        }
        this.color.property.bind(selectColorProperty, 'value');
        this.color.min.bind('chsp-color-min', 'value');
        this.color.max.bind('chsp-color-max', 'value');

        // ======= color palette
        const selectPalette = getByID<HTMLSelectElement>('chsp-palette');
        selectPalette.length = 0;
        for (const key in COLOR_MAPS) {
            selectPalette.options.add(new Option(key, key));
        }
        this.palette.bind(selectPalette, 'value');

        // ======= marker symbols
        const selectSymbolProperty = getByID<HTMLSelectElement>('chsp-symbol');
        // first option is 'default'
        selectSymbolProperty.options.length = 0;
        selectSymbolProperty.options.add(new Option('default', ''));
        for (const key in properties) {
            if (properties[key].string !== undefined) {
                selectSymbolProperty.options.add(new Option(key, key));
            }
        }
        this.symbol.bind(selectSymbolProperty, 'value');

        // ======= marker size
        const selectSizeProperty = getByID<HTMLSelectElement>('chsp-size');
        selectSizeProperty.options.length = 0;
        for (const key in properties) {
            selectSizeProperty.options.add(new Option(key, key));
        }
        this.size.property.bind(selectSizeProperty, 'value');
        this.size.factor.bind('chsp-size-factor', 'value');
        this.size.mode.bind('chsp-size-mode', 'value');
        this.size.reverse.bind('chsp-size-reverse', 'checked');
    }

    /** Get the colorscale to use for markers in the main plotly trace */
    public colorScale(): Plotly.ColorScale {
        return COLOR_MAPS[this.palette.value];
    }
}
