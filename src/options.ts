/**
 * @packageDocumentation
 * @module utils
 */

import assert from 'assert';

import { getByID, sendWarning } from './utils';

/**
 * Possible HTML attributes to attach to a setting
 */
// this is mostly to catch typo early. Feel free to add more!
type Attribute = 'value' | 'checked' | 'innerText';

/// Type mapping for options
interface OptionsTypeMap {
    // I would really like to use `'string': string` here as the mapping, but
    // this is forbidden by the typescript compiler (for reasons unclear to me,
    // it thinks that `OptionsTypeMap[T extends keyof OptionsTypeMap]` is
    // `never`). https://github.com/microsoft/TypeScript/issues/33014 seems to
    // be related.
    //
    // Using a tuple/array type (`[string]`) seems to trick the compiler/improve
    // the way it performs type inference.
    string: [string];
    int: [number];
    number: [number];
    boolean: [boolean];
}
type OptionsType = keyof OptionsTypeMap;
type OptionsValue<T extends OptionsType> = OptionsTypeMap[T][0];

function parse<T extends OptionsType>(type: T, value: string): OptionsValue<T> {
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

/** Possible origins of a option change: from JS code, or from DOM interaction */
export type OptionModificationOrigin = 'JS' | 'DOM';

/** A single DOM element that will be synchronized with the javascript value of a setting */
interface HTMLOptionElement {
    /** The DOM element that we should check for changes */
    element: HTMLElement;
    /** Which atrribute of the element should we check */
    attribute: Attribute;
    /** Store the function used to listen on the element to be able to remove it */
    listener: (e: Event) => void;
}

/**
 * Creates a validating function that checks value against a list of valid
 * entries. This is intended to be used as a `validate` callback on an
 * [[HTMLOption]].
 *
 * @param  valid list of valid values for the setting
 * @param  name  name of the setting for better error messages
 * @return       a function that can be used to validate a new setting value
 */
export function optionValidator(valid: string[], name = ''): (value: string) => void {
    return (value: string) => {
        if (valid.includes(value)) {
            return;
        }
        throw Error(
            `invalid property '${value}' for ${name}, are you sure the settings correspond to the current dataset?`
        );
    };
}

/**
 * Simple two-way data binding implementation to store settings in chemiscope.
 *
 * This class manages a setting single value, and bind it to HTML elements.
 * Whenever the value is changed from javascript code, or updated in any of the
 * linked HTML element, it is automatically updated everywhere.
 */
export class HTMLOption<T extends OptionsType> {
    /** the type of the value stored by this setting */
    public readonly type: T;
    /** Callback to validate the new value before propagating changes. */
    public validate: (value: OptionsValue<T>) => void;
    /** Additional callback to run whenever the setting value changes */
    public onchange: (value: OptionsValue<T>, origin: OptionModificationOrigin) => void;

    // storage of the actual value, this is separated from this.value to allow
    // intercepting set
    private _value: OptionsValue<T>;
    // list of bound HTML elements: any update to the setting `value` or any of
    // the element will change both the settting value and all of the linked
    // elements
    private _boundList: HTMLOptionElement[];

    /**
     * Create a new [[HTMLOption]] containing a value of the given type.
     * Possible type/values combinations are described in the [[OptionsTypeMap]]
     * interface.
     *
     * @param type  type of the setting
     * @param value initial value of the setting
     */
    constructor(type: T, value: OptionsValue<T>) {
        this.type = type;
        this._value = value;
        this._boundList = [];
        this.validate = () => {};
        this.onchange = () => {};

        Object.preventExtensions(this);
    }

    /** Get the value of this setting */
    public get value(): OptionsValue<T> {
        return this._value;
    }

    /** Set a new value for this setting */
    public set value(v: OptionsValue<T>) {
        this._update(v.toString(), 'JS');
    }

    /**
     * Add a new HTML element to the list of synchronized element.
     * `element.attribute` will be set to the setting value, and the setting
     * value will be updated everytime the 'change' event is emmited.
     *
     * @param  element   HTML DOM element to watch for updates, or string id
     *                   of such element
     * @param  attribute attribute of the HTML element to use as value
     */
    public bind(element: HTMLElement | string, attribute: Attribute): void {
        if (typeof element === 'string') {
            element = getByID(element);
        }
        element = element as HTMLElement;

        const listener = (event: Event) => {
            assert(event.target !== null);
            // eslint-disable-next-line @typescript-eslint/no-unsafe-member-access, @typescript-eslint/no-unsafe-call, @typescript-eslint/no-explicit-any
            this._update((event.target as any)[attribute].toString(), 'DOM');
        };

        // eslint-disable-next-line @typescript-eslint/no-explicit-any, @typescript-eslint/no-unsafe-member-access
        (element as any)[attribute] = this._value;
        element.addEventListener('change', listener);

        this._boundList.push({ element, attribute, listener });
    }

    /**
     * Disable all HTML elements linked to this by setting the `disabled`
     * attribute to `true` if it exists on the element. This should only
     * disable HTMLInputElement.
     */
    public disable(): void {
        for (const bound of this._boundList) {
            const element = (bound.element as unknown) as Record<string, unknown>;
            if ('disabled' in element) {
                element.disabled = true;
            }
        }
    }

    /**
     * Enable all HTML elements linked to this by setting the `disabled`
     * attribute to `false` if it exists on the element. This should only
     * enable HTMLInputElement.
     */
    public enable(): void {
        for (const bound of this._boundList) {
            const element = (bound.element as unknown) as Record<string, unknown>;
            if ('disabled' in element) {
                element.disabled = false;
            }
        }
    }

    /**
     * Remove all HTML elements bound to this setting, and the associated event
     * listeners.
     */
    public unbindAll(): void {
        for (const bound of this._boundList) {
            bound.element.removeEventListener('change', bound.listener);
        }
        this._boundList = [];
    }

    // Actually perform the update of the value and all linked elements. The
    // value is passed around as a string, since that's the easiest way to merge
    // different data types coming from the DOM.
    private _update(value: string, origin: OptionModificationOrigin) {
        const updated = parse<T>(this.type, value);
        this.validate(updated);

        this._value = updated;
        for (const bound of this._boundList) {
            // eslint-disable-next-line @typescript-eslint/no-unsafe-member-access, @typescript-eslint/no-explicit-any
            (bound.element as any)[bound.attribute] = updated;
        }
        this.onchange(updated, origin);
    }
}

/**
 * Type defintion for saveSettings output: should be a simple object containing
 * either value, or a nested SavedSettings.
 */
// eslint-disable-next-line @typescript-eslint/no-empty-interface
export interface SavedSettings extends Record<string, string | number | boolean | SavedSettings> {}

/**
 * Callback function to use with [[OptionsGroup.foreachSetting]]
 */
// eslint-disable-next-line @typescript-eslint/no-explicit-any
export type OptionsCallback = (keys: string[], setting: HTMLOption<any>) => void;

/**
 * Abstract base class to use for a group of settings.
 *
 * This class implement saving current settings as [[SavedSettings]]; and
 * applying saved settings to the setting group.
 *
 * # Example
 *
 * ```typescript
 * class MyOptions extends OptionsGroup {
 *     public cats: HTMLOption<'int'>;
 *     public dogs: {
 *          husky: HTMLOption<'string'>;
 *          labrador: HTMLOption<'string'>;
 *     };
 *
 *     constructor() {
 *         super();
 *         this.cats = new HTMLOption('int', 3);
 *         this.dogs = {
 *             husky: new HTMLOption('string', 'a cold dog'),
 *             labrador: new HTMLOption('string', 'a long dog'),
 *         };
 *     }
 * }
 *
 * const settings = new MyOptions();
 *
 * settings.saveSettings();
 * // {cats: 3, dogs: {husky: 'a cold dog', labrador: 'a long dog'}}
 *
 * settings.applySettings({cats: 66, dogs: {husky: 'a good dog'}});
 * settings.saveSettings();
 * // {cats: 66, dogs: {husky: 'a good dog', labrador: 'a long dog'}}
 * ```
 */
export abstract class OptionsGroup {
    /**
     * Save the current values of all HTMLOption properties of the class,
     * including nested ones.
     *
     * Properties which name starts with an underscore are ignored.
     *
     * @return An object with the same structure as this class containing the
     *         values of all settings.
     */
    public saveSettings(): SavedSettings {
        const settings = {};
        /* eslint-disable */
        this.foreachOption((keys, option) => {
            assert(keys.length >= 1);
            const value = option.value;

            let root = settings as any;
            for (const key of keys.slice(0, keys.length - 1)) {
                if (!(key in root)) {
                    root[key] = {};
                }
                root = root[key];
            }
            const lastkey = keys[keys.length - 1];
            root[lastkey] = value;
        });
        /* eslint-enable */
        return settings;
    }

    /**
     * Set values from `settings` to the [[HTMLOption]] properties of this class,
     * matching the properties names.
     *
     * Properties starting with an underscore are ignored.
     */
    public applySettings(settings: SavedSettings): void {
        // make a copy of the settings since we will be changing it below
        const copy = JSON.parse(JSON.stringify(settings)) as SavedSettings;

        this.foreachOption((keys, option) => {
            /* eslint-disable */
            assert(keys.length >= 1);
            let root = copy as any;
            let parent;
            for (const key of keys.slice(0, keys.length - 1)) {
                if (!(key in root)) {
                    // this key is missing from the settings
                    return;
                }
                parent = root;
                root = root[key];
            }
            const lastkey = keys[keys.length - 1];

            if (lastkey in root) {
                option.value = root[lastkey];

                // remove used keys from the settings to be able to warn on
                // unused keys
                delete root[lastkey];
                if (parent !== undefined && Object.keys(root).length === 0) {
                    // if we removed all keys from a sub-object, remove the sub-object
                    assert(keys.length >= 2);
                    delete parent[keys[keys.length - 2]];
                }
            }
            /* eslint-enable */
        });

        if (Object.keys(copy).length !== 0) {
            sendWarning(`ignored unkown settings '${JSON.stringify(copy)}'`);
        }
    }

    /**
     * Call the given callback on each setting inside the given SettingGroup.
     *
     * Keys starting with an underscore character are ignored.
     *
     * @param settings group of settings
     * @param callback callback operating on a single setting
     */
    protected foreachOption(callback: OptionsCallback): void {
        foreachOptionImpl(this as Record<string, unknown>, callback, []);
    }
}

/**
 * Recursive implementation of foreachOption
 *
 * This function looks through all properties of `options`, ignoring the ones
 * starting with an underscore. If the property is an instance of
 * [[HTMLOption]], the `callback` is called on the option, together with the
 * keys used to access the property from the root object.
 *
 * @param options  object containing the options
 * @param callback function to call on each HTMLOption
 * @param keys     current keys to access the `settings` object from the root
 */
function foreachOptionImpl(
    options: Record<string, unknown>,
    callback: OptionsCallback,
    keys: string[] = []
): void {
    if (keys.length > 10) {
        throw Error('setting object is too deep');
    }

    for (const key in options) {
        if (key.startsWith('_')) {
            continue;
        }

        const currentKeys = keys.concat([key]);
        const element = options[key];
        if (element instanceof HTMLOption) {
            callback(currentKeys, element);
        } else if (typeof element === 'object' && element !== null) {
            foreachOptionImpl(element as Record<string, unknown>, callback, currentKeys);
        }
    }
}
