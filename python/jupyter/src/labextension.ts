/* eslint-disable */
import { IJupyterWidgetRegistry } from '@jupyter-widgets/base';
import * as jquery from 'jquery';

// Force export jquery as a global for 3Dmol to find and use.
//
// When loading code in jupyter lab, jquery does not export itself as a global:
// https://github.com/jquery/jquery/blob/5d5ea015114092c157311c4948f7cc3d8c8e7f8a/src/wrapper.js#L25
// calls the factory with `noGlobal = true`, preventing the creation of `window.$`:
// https://github.com/jquery/jquery/blob/5d5ea015114092c157311c4948f7cc3d8c8e7f8a/src/exports/global.js#L26-L28
(window as any).$ = jquery;

// Use require instead of import since otherwise webpack moves the import before
// installing jquery above
const widgetExports = require('./widget');

const PACKAGE_NAME = 'chemiscope';
const PACKAGE_VERSION = require('../../../package.json').version;

const PLUGIN = {
    id: 'chemiscope:plugin',
    requires: [IJupyterWidgetRegistry],
    activate: (app: unknown, registry: IJupyterWidgetRegistry) => {
        registry.registerWidget({
            name: PACKAGE_NAME,
            version: PACKAGE_VERSION,
            exports: widgetExports,
        });
    },
    autoStart: true,
};

export default PLUGIN;
