/**
 * @packageDocumentation
 * @module utils
 */

import {makeDraggable} from './draggable';
import {Indexes, EnvironmentIndexer, DisplayMode} from './indexer';
import {sendWarning, addWarningHandler} from './warnings';

function generateGUID() {
    return "xxxxxxxx-xxxx-4xxx-yxxx-xxxxxxxxxxxx".replace(/[xy]/g, (c) => {
        var r = Math.random() * 16 | 0, v = c === "x" ? r : (r & 0x3 | 0x8);
        return v.toString(16);
    });
}

/** Get an HTML element by id */
function getByID<HTMLType>(id: string): HTMLType {
    const e = document.getElementById(id);
    if (e === null) {
        throw Error(`unable to get element with id ${id}`);
    }
    return e as unknown as HTMLType;
}

export {
    makeDraggable,
    sendWarning,
    addWarningHandler,
    generateGUID,
    getByID,
    Indexes,
    DisplayMode,
    EnvironmentIndexer,
}
