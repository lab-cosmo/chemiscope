/**
 * @packageDocumentation
 * @module utils
 */

/**
 * Possible HTML attributes to attach to a setting
 */
// this is mostly to catch typo early. Feel free to add more!
type Attribute = 'value' | 'checked' | 'innerText';

/// Type mapping for settings
interface SettingsTypeMap {
    // I would really like to use `'string': string` here as the mapping, but
    // this is forbidden by the typescript compiler (for reasons unclear to me,
    // it thinks that `SettingsTypeMap[T extends keyof SettingsTypeMap]` is
    // `never`). https://github.com/microsoft/TypeScript/issues/33014 seems to
    // be related.
    //
    // Using a tuple/array type (`[string]`) seems to trick the compiler/improve
    // the way it performs type inference.
    'string': [string];
    'int': [number];
    'number': [number];
    'boolean': [boolean];
}
type SettingsType = keyof SettingsTypeMap;
type SettingsValue<T extends SettingsType> = SettingsTypeMap[T][0];

function parse<T extends SettingsType>(type: T, value: string): SettingsValue<T> {
    if (type === 'string') {
        return value;
    } else if (type === 'int') {
        // FIXME: this parses `134hdvao` as `134`
        return parseInt(value, 10);
    } else if (type === 'number') {
        // FIXME: this parses `1.34hdvao` as `1.34`
        return parseFloat(value);
    } else if (type === 'boolean') {
        if (value === 'false') {
            return false;
        } else if (value === 'true') {
            return true;
        } else {
            throw Error(`invalid value for boolean: ${value}`);
        }
    }

    throw Error(`unknown type '${type}' passed to parse`);
}

/**
 * Simple two-way data binding implementation to store settings in chemiscope.
 *
 * This class manages a setting single value, and bind it to HTML elements.
 * Whenever the value is changed from javascript code, or updated in any of the
 * linked HTML element, it is automatically updated everywhere.
 */
export class HTMLSetting<T extends SettingsType> {
    /** The value of the setting, accessible to JS/TS code */
    declare public value: SettingsValue<T>;
    /** Additional callback to run whenever the setting value changes */
    public onchange: (value: SettingsValue<T>) => void;

    // Store the expected type of the value
    private _type: T;
    // Store the value itself
    private _value: SettingsValue<T>;
    // Store linked HTML elements: any update to the setting `value` or any of
    // the element will change both the settting value and all of the linked
    // elements
    private _elements: Array<[HTMLElement, Attribute]>;

    /**
     * Create a new [[HTMLSetting]] containing a value of the given type.
     * Possible type/values combinations are described in the [[SettingsTypeMap]]
     * interface.
     *
     * @param type  type of the setting
     * @param value initial value of the setting
     */
    constructor(type: T, value: SettingsValue<T>) {
        this._type = type;
        this._value = value;
        this._elements = [];
        this.onchange = () => {};

        // uses a Proxy on a frozen object to catch get/set events
        Object.preventExtensions(this);
        return new Proxy(this, {
            get: (_, key) => {
                if (key === 'value') {
                    return this._value;
                } else {
                    return Reflect.get(this, key);
                }
            },
            set: (_, key, v) => {
                if (key === 'value') {
                    this._update(v.toString());
                    return true;
                } else {
                    return Reflect.set(this, key, v);
                }
            },
        });
    }

    /**
     * Add a new HTML element to the list of synchronized element.
     * `element.attribute` will be set to the setting value, and the setting
     * value will be updated everytime the 'change' event is emmited.
     *
     * @param  element   HTML DOM element to watch for updates
     * @param  attribute attribute of the HTML element to use as value
     */
    public bind(element: HTMLElement, attribute: Attribute): void {
        this._elements.push([element, attribute]);
        (element as any)[attribute] = this._value;
        element.addEventListener('change', () => {
            this._update((element as any)[attribute].toString());
        });
    }

    /**
     * Disable all HTML elements linked to this by setting the `disabled`
     * attribute to `true` if it exists on the element. This should only
     * disable HTMLInputElement.
     */
    public disable(): void {
        for (const pair of this._elements) {
            const element = pair[0];
            if ('disabled' in element) {
                (element as any).disabled = true;
            }
        }
    }

    /**
     * Enable all HTML elements linked to this by setting the `disabled`
     * attribute to `false` if it exists on the element. This should only
     * enable HTMLInputElement.
     */
    public enable(): void {
        for (const pair of this._elements) {
            const element = pair[0];
            if ('disabled' in element) {
                (element as any).disabled = false;
            }
        }
    }

    // Actually perform the update of the value and all linked elements. The
    // value is passed around as a string, since that's the easiest way to merge
    // different data types coming from the DOM.
    private _update(value: string) {
        this._value = parse<T>(this._type, value);
        for (const [element, attribute] of this._elements) {
            (element as any)[attribute] = this._value;
        }
        this.onchange(this._value);
    }
}
