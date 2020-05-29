/**
 * @packageDocumentation
 * @module utils
 */

import {makeDraggable} from './draggable';
import {DisplayMode, EnvironmentIndexer, Indexes} from './indexer';
import {foreachSetting, HTMLSetting, SettingGroup, SettingModificationOrigin} from './settings';
import {addWarningHandler, sendWarning} from './warnings';

function generateGUID() {
    return 'xxxxxxxx-xxxx-4xxx-yxxx-xxxxxxxxxxxx'.replace(/[xy]/g, (c) => {
        // tslint:disable-next-line:no-bitwise
        const r = Math.random() * 16 | 0;
        // tslint:disable-next-line:no-bitwise
        const v = c === 'x' ? r : (r & 0x3 | 0x8);
        return v.toString(16);
    });
}

/** Get an HTML element by id */
function getByID<HTMLType = HTMLElement>(id: string): HTMLType {
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
    foreachSetting,
    generateGUID,
    getByID,
    Indexes,
    DisplayMode,
    EnvironmentIndexer,
    HTMLSetting,
    SettingGroup,
    SettingModificationOrigin,
};
