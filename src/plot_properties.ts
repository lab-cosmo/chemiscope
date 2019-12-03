import {Property} from './dataset'

export interface NumericProperty {
    /// All properties have numeric values. String properties also have a string
    /// interner
    values: number[]
    string?: StringInterner;
}

/// A simple string => number association, giving a single id to each string
export class StringInterner {
    private _values: string[];

    constructor() {
        this._values = [];
    }

    /// Get the index associated with the given string `value`
    public get(value: string): number {
        let index = this._values.findIndex(x => x === value);
        if (index === -1) {
            index = this._values.length;
            this._values.push(value);
        }

        return index;
    }

    /// Get the string associated with the given `index`
    public string(index: number): string {
        if (index >= this._values.length) {
            throw Error("requested unknown string from interner")
        }
        return this._values[index]
    }
}

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

/// Data storage for maps, differenciating between numeric/string and
/// structure/atom properties
///
/// String properties can be used to assing classes to each data point. They are
/// mapped to a numeric value using the `StringInterner` for easier
/// vizualization.
export class PlotProperties {
    /// Structure properties
    public structure: {
        [name: string]: NumericProperty
    }

    /// Atomic properties
    public atom: {
        [name: string]: NumericProperty
    }

    constructor(properties: {[name: string]: Property}) {
        this.structure = {}
        this.atom = {}

        for (const name in properties) {
            if (properties[name].target == 'structure') {
                this.structure[name] = propertyToNumeric(properties[name].values);
            } else if (properties[name].target == 'atom') {
                this.atom[name] = propertyToNumeric(properties[name].values);
            }
        }

        // sanity checks
        checkSize('structure', this.structure);
        checkSize('atom', this.atom);
    }
}
