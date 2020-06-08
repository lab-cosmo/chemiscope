/**
 * @packageDocumentation
 * @module utils
 */

import {makeDraggable} from './draggable';
import {DisplayMode, EnvironmentIndexer, Indexes} from './indexer';
import {foreachSetting, HTMLSetting, SettingGroup, SettingModificationOrigin} from './settings';
import {addWarningHandler, sendWarning} from './warnings';

const LIGHT_COLORS = ['White', 'Snow', 'Honeydew', 'MintCream',
                     'Azure', 'AliceBlue', 'GhostWhite', 'WhiteSmoke',
                     'Seashell', 'Beige', 'OldLace', 'FloralWhite',
                     'Ivory', 'AntiqueWhite', 'Linen', 'LavenderBlush',
                     'MistyRose'];
const STANDARD_COLORS = ['Red', 'Yellow', 'Green', 'Blue',
                   'Orange', 'Aqua', 'Purple', 'Teal',
                   'White', 'Gray', 'Black',
                   'Maroon', 'Olive', 'Silver', 'Lime',
                   'Navy', 'Fuchsia'];

// to be replaced later
export function getRandomColor() {
  const letters = '0123456789ABCDEF';
  let color = '#';
  for (let i = 0; i < 6; i++) {
    color += letters[Math.floor(Math.random() * 16)];
  }
  return color;
}

export function getNextColor(colors: string[]) {
    for (const color of STANDARD_COLORS) {
      if ( colors.indexOf(color) < 0) {
        return color;
      }
    }
    return STANDARD_COLORS[1];
}

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
