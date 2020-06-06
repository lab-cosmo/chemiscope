/**
 * @packageDocumentation
 * @module utils
 */

import assert from 'assert';

import {getByID} from './index';

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

/** Possible origins of a setting change: from JS code, or from DOM interaction */
export type SettingModificationOrigin = 'JS' | 'DOM';

/** A single DOM element that will be synchronized with the javascript value of a setting */
interface HTMLSettingElement {
    /** The DOM element that we should check for changes */
    element: HTMLElement;
    /** Which atrribute of the element should we check */
    attribute: Attribute;
    /** Store the function used to listen on the element to be able to remove it */
    listener: (e: Event) => void;
}

/** Multiple settings stored together as (potentially nested) javascript object */
export interface SettingGroup {
    [key: string]: HTMLSetting<any> | SettingGroup;
}

/**
 * Call the given callback on each setting inside the given SettingGroup
 *
 * @param settings group of settings
 * @param callback callback operating on a single setting
 */
export function foreachSetting(settings: SettingGroup, callback: (s: HTMLSetting<any>) => void): void {
    for (const key in settings) {
        const element = settings[key];
        if (element instanceof HTMLSetting) {
            callback(element);
        } else {
            foreachSetting(element as SettingGroup, callback);
        }
    }
}

/**
 * Creates a valitating function that checks value against a list of valid
 * entries. This is intended to be used as a `validate` callback on an
 * [[HTMLSetting]].
 *
 * @param  valid list of valid values for the setting
 * @param  name  name of the setting for better error messages
 * @return       a function that can be used to validate a new setting value
 */
export function settingsValidator<T>(valid: T[], name = ''): (value: T) => void {
    return (value: T) => {
        if (valid.includes(value)) {
            return;
        }
        throw Error(`invalid property '${value}' for ${name}`);
    };
}

/**
 * Simple two-way data binding implementation to store settings in chemiscope.
 *
 * This class manages a setting single value, and bind it to HTML elements.
 * Whenever the value is changed from javascript code, or updated in any of the
 * linked HTML element, it is automatically updated everywhere.
 */
export class HTMLSetting<T extends SettingsType> {
    /** the type of the value stored by this setting */
    public readonly type: T;
    /** Callback to validate the new value before propagating changes. */
    public validate: (value: SettingsValue<T>) => void;
    /** Additional callback to run whenever the setting value changes */
    public onchange: (value: SettingsValue<T>, origin: SettingModificationOrigin) => void;

    // storage of the actual value, this is separated from this.value to allow
    // intercepting set
    private _value: SettingsValue<T>;
    // list of bound HTML elements: any update to the setting `value` or any of
    // the element will change both the settting value and all of the linked
    // elements
    private _boundList: HTMLSettingElement[];

    /**
     * Create a new [[HTMLSetting]] containing a value of the given type.
     * Possible type/values combinations are described in the [[SettingsTypeMap]]
     * interface.
     *
     * @param type  type of the setting
     * @param value initial value of the setting
     */
    constructor(type: T, value: SettingsValue<T>) {
        this.type = type;
        this._value = value;
        this._boundList = [];
        this.validate = () => {};
        this.onchange = () => {};

        Object.preventExtensions(this);
    }

    /** Get the value of this setting */
    public get value(): SettingsValue<T> {
        return this._value;
    }

    /** Set a new value for this setting */
    public set value(v: SettingsValue<T>) {
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
            this._update((event.target as any)[attribute].toString(), 'DOM');
        };

        (element as any)[attribute] = this._value;
        element.addEventListener('change', listener);

        this._boundList.push({element, attribute, listener});
    }

    /**
     * Disable all HTML elements linked to this by setting the `disabled`
     * attribute to `true` if it exists on the element. This should only
     * disable HTMLInputElement.
     */
    public disable(): void {
        for (const bound of this._boundList) {
            if ('disabled' in bound.element) {
                (bound.element as any).disabled = true;
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
            if ('disabled' in bound.element) {
                (bound.element as any).disabled = false;
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
    private _update(value: string, origin: SettingModificationOrigin) {
        const updated = parse<T>(this.type, value);
        this.validate(updated);

        this._value = updated;
        for (const bound of this._boundList) {
            (bound.element as any)[bound.attribute] = updated;
        }
        this.onchange(updated, origin);
    }
}
