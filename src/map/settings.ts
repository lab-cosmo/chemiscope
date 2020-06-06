/**
 * @packageDocumentation
 * @module settings
 */

import assert from 'assert';

import {HTMLSetting, settingsValidator} from '../utils';

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
    };

    private readonly _defaultProperties: Partial<MapPresets>;

    constructor(properties: string[], presets: Partial<MapPresets> = {}) {
        if (properties.length < 2) {
            throw Error('we need at least two properties to plot in the map');
        }
        this._defaultProperties = {
            x: {property: properties[0]},
            y: {property: properties[1]},
        };

        if (properties.length >= 2) {
            this._defaultProperties.color = {property: properties[2]};
        }

        this.x = new AxisSettings(properties);
        this.y = new AxisSettings(properties);
        // For z and color, '' is a valid value
        this.z = new AxisSettings(properties.concat(['']));
        this.color = new AxisSettings(properties.concat(['']));

        this.symbol = new HTMLSetting('string', '');
        this.symbol.validate = settingsValidator(properties.concat(['']), 'symbol');

        this.palette = new HTMLSetting('string', '');
        this.symbol.validate = settingsValidator(properties.concat(['']), 'symbol');

        this.size = {
            factor : new HTMLSetting('number', 0),
            property : new HTMLSetting('string', ''),
        };
        this.size.property.validate = settingsValidator(properties.concat(['']), 'size');
        this.size.factor.validate = (value) => {
            if (value < 0 || value > 100) {
                throw Error(`size factor must be between 0 and 100, got ${value}`);
            }
        };

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
        this.size.factor.value = initial.size.factor;
        this.size.property.value = initial.size.property;
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
                property: this.size.property.value,
            },
            symbol: this.symbol.value,
        };
    }
}
