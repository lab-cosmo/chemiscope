/* eslint-disable */
const { JSDOM } = require('jsdom');

/**
 * Copy all JS properties from src to target
 */
function copyProps(src: any, target: any) {
    Object.defineProperties(target, {
        ...Object.getOwnPropertyDescriptors(src),
        ...Object.getOwnPropertyDescriptors(target),
    });
}

/**
 * Setup JSDom and the corresponding globals on `window`
 */
export default function () {
    const jsdom = new JSDOM('<!doctype html><html><body></body></html>', {
        pretendToBeVisual: true,
    });
    const { window } = jsdom;

    global.window = window;
    global.document = window.document;
    global.navigator = {
        userAgent: 'node.js',
    } as any;

    copyProps(window, global);
}
