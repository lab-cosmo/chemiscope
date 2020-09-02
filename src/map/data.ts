/**
 * @packageDocumentation
 * @module map
 */

import { Property } from '../dataset';
import { sendWarning } from '../utils';
import { DisplayMode } from '../indexer';

/** @hidden
 * Properties turned into numeric values to be displayed on the map.
 */
export interface NumericProperty {
    /** values for the property */
    values: number[];
    /** string interner if the property was a string-valued property */
    string?: StringInterner;
    /** unit of the property */
    units?: string;
}

/** @hidden
 * A simple string => number association, giving a single numeric id to each
 * string.
 */
export class StringInterner {
    private _values: string[];

    /** Create a new empty [[StringInterner]] */
    constructor() {
        this._values = [];
    }

    /** Get the index associated with the given string `value` */
    public get(value: string): number {
        let index = this._values.findIndex((x) => x === value);
        if (index === -1) {
            index = this._values.length;
            this._values.push(value);
        }

        return index;
    }

    /** Get the string associated with the given `index` */
    public string(index: number): string {
        if (index >= this._values.length) {
            throw Error('requested unknown string from interner');
        }
        return this._values[index];
    }

    /** Get all the strings known by this interner */
    public strings(): string[] {
        return this._values;
    }
}

/** Transform a property to a numeric property */
function propertyToNumeric(name: string, property: Property): NumericProperty {
    const prop_type = typeof property.values[0];
    if (prop_type === 'number') {
        return {
            values: property.values as number[],
            units: property.units,
        };
    } else if (prop_type === 'string') {
        const interner = new StringInterner();
        const values = [];
        for (const value of property.values as string[]) {
            const numeric = interner.get(value);
            values.push(numeric);
            // string properties are assumed to be categories, and each one of
            // them have an associated trace to be able to display the symbols
            // legend. Having more than a handfull of traces drastically slows
            // down the whole page.
            if (numeric > 20) {
                throw Error(
                    `the '${name}' property contains more than 50 different values, \
                    it can not be interpreted as categories`
                );
            }
        }

        return {
            string: interner,
            values: values,
            units: property.units,
        };
    } else {
        throw Error(`unexpected property type '${prop_type}'`);
    }
}

/** Sanity check that all properties have the same size */
function checkSize(name: string, properties: { [key: string]: NumericProperty }) {
    let size;
    let initial;
    for (const key in properties) {
        if (size === undefined) {
            size = properties[key].values.length;
            initial = key;
            continue;
        }

        if (properties[key].values.length !== size) {
            throw Error(
                `${name} property '${key}' do not have the same size as the first property '${
                    initial as string
                }': expected ${size}, got ${properties[key].values.length}`
            );
        }
    }
}

export interface NumericProperties {
    [name: string]: NumericProperty;
}

/** @hidden
 * Data storage for maps, differenciating between numeric/string and
 * structure/atom properties
 *
 * String properties can be used to assing classes to each data point. They are
 * mapped to a numeric value using the `StringInterner` for easier
 * visualization.
 */
export class MapData {
    /** Structure properties */
    public structure: NumericProperties;
    /** Atomic properties */
    public atom: NumericProperties;
    /** Maximal number of symbols (i.e. different values in string properties) in this dataset */
    public maxSymbols: number;

    /** Create a new [[MapData]] containing values from the given properties */
    constructor(properties: { [name: string]: Property }) {
        this.structure = {};
        this.atom = {};
        this.maxSymbols = -1;

        for (const name in properties) {
            let property;
            try {
                property = propertyToNumeric(name, properties[name]);
            } catch (e) {
                sendWarning(`warning: ${(e as Error).message}`);
                continue;
            }
            this[properties[name].target][name] = property;
            if (property.string !== undefined) {
                this.maxSymbols = Math.max(this.maxSymbols, property.string.strings().length);
            }
        }

        // sanity checks
        checkSize('structure', this.structure);
        checkSize('atom', this.atom);
    }

    public length(mode: DisplayMode): number {
        return this[mode][Object.keys(this[mode])[0]].values.length;
    }
}
