/**
 * @packageDocumentation
 * @module settings
 */

import assert from 'assert';

import {HTMLSetting, PositioningCallback} from '../utils';
import {getByID, makeDraggable, settingsValidator} from '../utils';
import {NumericProperties} from './data';

import {COLOR_MAPS} from './colorscales';

import BARS_SVG from '../static/bars.svg';
import HTML_SETTINGS from './settings.html';

export interface AxisPresets {
    /// Which property should we use for this axis
    property: string;
    /// Which scale (linear/log) should we use
    scale: string;
    /// The minimal value for this axis
    min: number;
    /// The maximal value for this axis
    max: number;
}

const AXIS_DEFAULTS: AxisPresets = {
     max: 0,
     min : 0,
     property : '',
     scale : 'linear',
};

/** HTML element holding settings for a given axis (x, y, z, color) */
export class AxisSettings {
    public property: HTMLSetting<'string'>;
    public scale: HTMLSetting<'string'>;
    public min: HTMLSetting<'number'>;
    public max: HTMLSetting<'number'>;

    constructor(validProperties: string[]) {
        this.max = new HTMLSetting('number', AXIS_DEFAULTS.max);
        this.min = new HTMLSetting('number', AXIS_DEFAULTS.min);
        this.property = new HTMLSetting('string', AXIS_DEFAULTS.property);
        this.scale = new HTMLSetting('string', AXIS_DEFAULTS.scale);

        this.property.validate = settingsValidator(validProperties, 'axis');
        this.scale.validate = settingsValidator(['linear', 'log'], 'axis scale');
    }

    /**
     * Applies presets, possibly filling in with default values
     */
    public applyPresets(presets: Partial<AxisPresets> = {}): void {
        const initial: AxisPresets = {
            ...AXIS_DEFAULTS,
            ...presets,
        };

        this.max.value = initial.max;
        this.min.value = initial.min;
        this.property.value = initial.property;
        this.scale.value = initial.scale;
    }

    /**
     * Dumps presets, in a way that can e.g. be serialized to json
     */
    public dumpPresets(): AxisPresets {
        return {
            max: this.max.value,
            min: this.min.value,
            property: this.property.value,
            scale: this.scale.value,
        };
    }

    /** Disable auxiliary settings (min/max/scale) related to this axis */
    public disable() {
        this.max.disable();
        this.min.disable();
        this.scale.disable();
    }

    /** Enable auxiliary settings (min/max/scale) related to this axis */
    public enable() {
        this.max.enable();
        this.min.enable();
        this.scale.enable();
    }
}

export interface MapPresets {
    x: Partial<AxisPresets>;
    y: Partial<AxisPresets>;
    z: Partial<AxisPresets>;
    color: Partial<AxisPresets>;
    palette: string;
    symbol: string;
    size: {
        property?: string,
        factor?: number,
        mode: string,
    };
}

const MAP_DEFAULTS: MapPresets = {
    color: AXIS_DEFAULTS,
    x: AXIS_DEFAULTS,
    y: AXIS_DEFAULTS,
    z: AXIS_DEFAULTS,

    palette: 'inferno',
    size: {
        factor: 50.0,
        mode: '',
        property: '',
    },
    symbol: '',
};

export class MapSettings {
    public x: AxisSettings;
    public y: AxisSettings;
    public z: AxisSettings;
    public color: AxisSettings;
    public palette: HTMLSetting<'string'>;
    public symbol: HTMLSetting<'string'>;
    public size: {
        property: HTMLSetting<'string'>;
        factor: HTMLSetting<'number'>;
        mode: HTMLSetting<'string'>;
    };

    /// Default properties used
    private readonly _defaultProperties: Partial<MapPresets>;
    /// The HTML element containing the settings modal
    private _settingsModal: HTMLElement;
    // Callback to get the initial positioning of the settings modal.
    private _positionSettingsModal: PositioningCallback;

    constructor(
        root: HTMLElement,
        properties: NumericProperties,
        positionSettings: PositioningCallback,
        presets: Partial<MapPresets> = {},
    ) {
        const propertiesName = Object.keys(properties);
        if (propertiesName.length < 2) {
            throw Error('we need at least two properties to plot in the map');
        }

        this._defaultProperties = {
            x: {property: propertiesName[0]},
            y: {property: propertiesName[1]},
        };

        if (propertiesName.length >= 2) {
            this._defaultProperties.color = {property: propertiesName[2]};
        }

        this.x = new AxisSettings(propertiesName);
        this.y = new AxisSettings(propertiesName);
        // For z and color, '' is a valid value
        this.z = new AxisSettings(propertiesName.concat(['']));
        this.color = new AxisSettings(propertiesName.concat(['']));

        this.symbol = new HTMLSetting('string', '');
        this.symbol.validate = settingsValidator(propertiesName.concat(['']), 'symbol');

        this.palette = new HTMLSetting('string', '');
        this.palette.validate = settingsValidator(Object.keys(COLOR_MAPS), 'palette');

        this.size = {
            factor : new HTMLSetting('number', 0),
            property : new HTMLSetting('string', ''),
            mode : new HTMLSetting('string', ''),
        };
        this.size.property.validate = settingsValidator(propertiesName.concat(['']), 'size');
        this.size.factor.validate = (value) => {
            if (value < 0 || value > 100) {
                throw Error(`size factor must be between 0 and 100, got ${value}`);
            }
        };

        this._positionSettingsModal = positionSettings;
        this._settingsModal = this._insertSettingsHTML(root);
        document.body.appendChild(this._settingsModal);

        this._bind(properties);
        this.applyPresets(presets);
    }

    /**
     * Applies presets, possibly filling in with default values
     */
    public applyPresets(presets: Partial<MapPresets> = {}): void {
        const initial: MapPresets = {
            ...MAP_DEFAULTS,
            ...this._defaultProperties,
            ...presets,
        };

        this.x.applyPresets(initial.x);
        this.y.applyPresets(initial.y);
        this.z.applyPresets(initial.z);
        this.color.applyPresets(initial.color);

        this.symbol.value = initial.symbol;
        this.palette.value = initial.palette;

        assert(initial.size.factor !== undefined);
        assert(initial.size.property !== undefined);
        assert(initial.size.mode !== undefined);
        this.size.factor.value = initial.size.factor;
        this.size.property.value = initial.size.property;
        this.size.mode.value = initial.size.mode;
    }

    /**
     * Dumps presets, in a way that can e.g. be serialized to json
     */
    public dumpPresets(): MapPresets {
        return {
            color: this.color.dumpPresets(),
            x: this.x.dumpPresets(),
            y: this.y.dumpPresets(),
            z: this.z.dumpPresets(),

            palette: this.palette.value,
            size: {
                factor: this.size.factor.value,
                mode: this.size.mode.value,
                property: this.size.property.value,
            },
            symbol: this.symbol.value,
        };
    }

    /**
     * Remove all HTML added by this [[MapSettings]] in the current document
     */
    public remove(): void {
        this._settingsModal.remove();
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
        const openSettings = template.content.firstChild!;
        root.append(openSettings);

        // TODO: set unique HTML id in the settings to allow multiple map in
        // the same page
        template.innerHTML = HTML_SETTINGS;
        const modal = template.content.firstChild! as HTMLElement;

        const modalDialog = modal.childNodes[1]! as HTMLElement;
        if (!modalDialog.classList.contains('modal-dialog')) {
            throw Error('internal error: missing modal-dialog class');
        }
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

                const {top, left} = this._positionSettingsModal(modalDialog.getBoundingClientRect());

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

    /**
     * Bind all settings to the corresponding HTML elements
     */
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
        // first option is 'default'
        selectSizeProperty.options.length = 0;
        selectSizeProperty.options.add(new Option('default', ''));
        for (const key in properties) {
            selectSizeProperty.options.add(new Option(key, key));
        }
        this.size.property.bind(selectSizeProperty, 'value');
        this.size.factor.bind('chsp-size-factor', 'value');
        this.size.mode.bind('chsp-size-mode', 'value');
    }
}
