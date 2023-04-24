/**
 * @packageDocumentation
 * @module utils
 */

import assert from 'assert';

export { makeDraggable } from './draggable';
export { addWarningHandler, sendWarning, WarningHandler } from './warnings';

/** Callback type to position an HTML element.
 *
 * The callback gets the current placement of the settings as a
 * [DOMRect](https://developer.mozilla.org/en-US/docs/Web/API/DOMRect), and
 * should return top and left positions in pixels, used with `position:
 * fixed`.
 */
export type PositioningCallback = (rect: DOMRect) => { top: number; left: number };

const STANDARD_COLORS = [
    'red',
    'yellow',
    'green',
    'blue',
    'orange',
    'aqua',
    'purple',
    'teal',
    'silver',
];

export function getNextColor(alreadyUsedColors: string[]): string {
    for (const color of STANDARD_COLORS) {
        if (!alreadyUsedColors.includes(color)) {
            return color;
        }
    }
    throw Error('required more colors than available');
}

/** Type to separate generic strings from GUID */
declare const tag: unique symbol;
export type GUID = string & { readonly [tag]: 'guid' };

/** Generate a new GUID */
export function generateGUID(): GUID {
    return 'xxxxxxxx-xxxx-4xxx-yxxx-xxxxxxxxxxxx'.replace(/[xy]/g, (c) => {
        const r = (Math.random() * 16) | 0;
        const v = c === 'x' ? r : (r & 0x3) | 0x8;
        return v.toString(16);
    }) as GUID;
}

/**
 * Get an HTML element by id, looking inside the `root` (by default the whole `document`).
 *
 * The generic parameter `HTMLType` can be used to cast the element to a given type (e.g. `getById<HTMLInputElement>("foo")`) . The element type is not checked by this function.
 *
 * @throws if there is not element with the given id.
 */
export function getByID<HTMLType = HTMLElement>(
    id: string,
    root: Document | HTMLElement | ShadowRoot = document
): HTMLType {
    let e;
    if (root instanceof HTMLElement) {
        e = root.querySelector(`#${id}`);
    } else {
        e = root.getElementById(id);
    }
    if (e === null) {
        throw Error(`unable to get element with id ${id}`);
    }
    return e as unknown as HTMLType;
}

export function enumerate<T>(iterable: Iterable<T>): Iterable<[number, T]> {
    const iterator = function* () {
        let index = 0;
        for (const element of iterable) {
            yield [index, element];
            index++;
        }
    };

    return {
        [Symbol.iterator]: iterator,
    } as Iterable<[number, T]>;
}

/**
 * Get the first key of a Map, potentially excluding a specific key value.
 *
 * If `excluding` is specified and the first key is equal to `excluding`; the
 * second key will be returned.
 *
 * @param map       the map where to look for keys
 * @param excluding do not use this specific value if it is the first key
 * @return the first or second key depending on `excluding`
 */
export function getFirstKey<K, V>(map: Map<K, V>, excluding?: K): K {
    assert(map.size >= 1);
    const keys = map.keys();
    const first = keys.next().value as K;
    if (excluding !== undefined) {
        if (first === excluding) {
            assert(map.size >= 2);
            return keys.next().value as K;
        }
    }
    return first;
}

// get the max/min of an array. Math.min(...array) fails with very large arrays
export function arrayMaxMin(values: number[]): { max: number; min: number } {
    let max = Number.NEGATIVE_INFINITY;
    let min = Number.POSITIVE_INFINITY;
    for (const value of values) {
        if (value > max && isFinite(value)) {
            max = value;
        }
        if (value < min && isFinite(value)) {
            min = value;
        }
    }
    assert(isFinite(min) && isFinite(max));
    return { max, min };
}

/** unreachable marker for cases that REALLY should not happen */
export function unreachable(): never {
    throw Error('unreachable code entered, this is a bug');
}

/** Returns the element based on ID or the element itself*/
export function getElement<HTMLElement>(element: string | HTMLElement): HTMLElement {
    if (typeof element !== 'string') {
        assert(element instanceof HTMLElement);
        return element;
    } else {
        return getByID<HTMLElement>(element);
    }
}

/**
 * check for the position of `value` in `array` using binary search.
 * The array must be sorted for this to work
 */
export function binarySearch(array: number[], value: number): number {
    let min = 0;
    let max = array.length - 1;
    let mid;
    while (min <= max) {
        mid = (min + max) >>> 1;
        if (array[mid] === value) {
            return mid;
        } else if (array[mid] < value) {
            min = mid + 1;
        } else {
            max = mid - 1;
        }
    }

    return -1;
}
