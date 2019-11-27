export interface MapInputData {
    [key: string]: number[] | string[];
}

export interface MapInput {
    name: string;
    data: MapInputData;
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

/// Data storage for maps, differenciating between numeric and string properties
///
/// String properties can be used to assing classes to each data point. They are
/// mapped to a numeric value using the `StringInterner` for easier vizualization.
export class MapData {
    /// All properties have a name, and numeric values. String properties also
    /// have a string interner
    [name: string]: {
        values: number[]
        string?: StringInterner;
    };

    constructor(data: MapInputData) {
        // check that all properties have the same size
        let size = -1;
        for (const key in data) {
            if (size === -1) {
                size = data[key].length;
            }

            if (data[key].length !== size) {
                throw Error("not all properties have the same size")
            }
        }

        if (size === 0) {
            return;
        }

        for (const name in data) {
            const prop_type = typeof(data[name][0]);
            if (prop_type === "number") {
                this[name] = {
                    values: data[name] as number[],
                }
            } else if (prop_type === "string") {
                const interner = new StringInterner();
                const values = [];
                for (const value of (data[name] as string[])) {
                    values.push(interner.get(value));
                }

                this[name] = {
                    values: values,
                    string: interner,
                }
            } else {
                throw Error(`unexpected property type '${prop_type}'`);
            }
        }
    }
}
