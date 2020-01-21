/**
 * @packageDocumentation
 * @module map
 */

import {Property} from '../dataset'

/** @hidden
 * Properties turned into numeric values to be displayed on the map.
 */
export interface NumericProperty {
    /** values for the property */
    values: number[]
    /** string interner if the property was a string-valued property */
    string?: StringInterner;
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
        let index = this._values.findIndex(x => x === value);
        if (index === -1) {
            index = this._values.length;
            this._values.push(value);
        }

        return index;
    }

    /** Get the string associated with the given `index` */
    public string(index: number): string {
        if (index >= this._values.length) {
            throw Error("requested unknown string from interner")
        }
        return this._values[index]
    }

    /** Get all the strings known by this interner */
    public strings(): string[] {
        return this._values;
    }
}

/** Transform a property to a numeric property */
function propertyToNumeric(property: number[] | string[]): NumericProperty {
    const prop_type = typeof property[0];
    if (prop_type === "number") {
        return {
            values: property as number[],
        }
    } else if (prop_type === "string") {
        const interner = new StringInterner();
        const values = [];
        for (const value of (property as string[])) {
            values.push(interner.get(value));
        }

        return {
            values: values,
            string: interner,
        }
    } else {
        throw Error(`unexpected property type '${prop_type}'`);
    }
}

/** Sanity check that all properties have the same size */
function checkSize(name: string, properties: { [key: string]: NumericProperty }) {
    let size = undefined;
    for (const key in properties) {
        if (size === undefined) {
            size = properties[key].values.length;
            continue
        }

        if (properties[key].values.length !== size) {
            throw Error(`${name} properties do not all have the same size`);
        }
    }
}

/** @hidden
 * Data storage for maps, differenciating between numeric/string and
 * structure/atom properties
 *
 * String properties can be used to assing classes to each data point. They are
 * mapped to a numeric value using the `StringInterner` for easier
 * vizualization.
 */
export class MapData {
    /** Structure properties */
    public structure: {
        [name: string]: NumericProperty
    }

    /** Atomic properties */
    public atom: {
        [name: string]: NumericProperty
    }

    /** Maximal number of symbols (i.e. different values in string properties) in this dataset */
    public maxSymbols: number;

    /** Create a new [[MapData]] containing values from the given properties */
    constructor(properties: {[name: string]: Property}) {
        this.structure = {}
        this.atom = {}
        this.maxSymbols = -1;

        for (const name in properties) {
            if (properties[name].target == 'structure') {
                const property = propertyToNumeric(properties[name].values);
                this.structure[name] = property;
                if (property.string !== undefined) {
                    this.maxSymbols = Math.max(this.maxSymbols, property.string.strings().length)
                }
            } else if (properties[name].target == 'atom') {
                const property = propertyToNumeric(properties[name].values);
                this.atom[name] = property;
                if (property.string !== undefined) {
                    this.maxSymbols = Math.max(this.maxSymbols, property.string.strings().length)
                }
            }
        }

        // sanity checks
        checkSize('structure', this.structure);
        checkSize('atom', this.atom);
    }
}
