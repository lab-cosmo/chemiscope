/* eslint-disable @typescript-eslint/no-unsafe-member-access */
import { Config, ConfigOptions } from 'karma';
import webpack from 'webpack';

import { WEBPACK_CONFIG } from './webpack.config';

/**
 * Small webpack plugin to remove typescript declaration files (.d.ts) from the
 * list of assets.
 *
 * karma-webpack is using this list to know which file to load for tests, which
 * contains both the path to test files compiled to javascript, and the path to
 * typescript declaration files. Webpack output path is set to a temporary
 * directory, and the javascript files are inside this directory. The typescript
 * declaration files are not in this temporary directory (they are emitted in
 * `chemiscope/dist/` by tsc), and for a still unknown reason that can cause
 * karma not to found these file in some very specific cases (macOS and/or ZSH
 * seems to be required to hit this issue).
 *
 * Since we do not care about the `.d.ts` files in karma, this plugin works
 * around the issue by removing said files from the assets list. This can break
 * if webpack plugin API changes, or if karma-webpack stops using the assets.
 */
class RemoveDeclarationsFromAssets implements webpack.WebpackPluginInstance {
    apply(compiler: webpack.Compiler) {
        compiler.hooks.done.tap('RemoveDeclarationsFromAssets', (stats) => {
            const toRemove = [];
            for (const path in stats.compilation.assets) {
                if (path.endsWith('.d.ts')) {
                    toRemove.push(path);
                }
            }

            for (const key of toRemove) {
                delete stats.compilation.assets[key];
            }
        });
    }
}

WEBPACK_CONFIG.plugins?.push(new RemoveDeclarationsFromAssets());

module.exports = (config: Config) => {
    config.set({
        browserNoActivityTimeout: 8000,
        client: {
            mocha: {
                timeout: 8000,
            },
        },

        files: [
            // FIXME: we should not have to manually load jquery, but we
            // currently don't include it in the main bundle
            'node_modules/jquery/dist/jquery.min.js',
            'tests/**/*.test.ts',
        ],
        frameworks: ['webpack', 'mocha', 'detectBrowsers'],

        preprocessors: {
            'tests/**/*.test.ts': 'webpack',
        },
        reporters: ['progress'],
        singleRun: true,

        webpack: WEBPACK_CONFIG,
        webpackMiddleware: {
            stats: 'errors-only',
        },

        detectBrowsers: {
            postDetection: function (availableBrowsers: string[]) {
                // Remove IE
                const IEindex = availableBrowsers.indexOf('IE');
                if (IEindex !== -1) {
                    availableBrowsers.splice(IEindex);
                }

                // Rename Safari to use SafariNative karma launcher
                const SafariIndex = availableBrowsers.indexOf('Safari');
                if (SafariIndex !== -1) {
                    availableBrowsers[SafariIndex] = 'SafariNative';
                }

                return availableBrowsers;
            },
            // we can not enable headless mode since firefox does not support
            // WebGL in this case (https://bugzilla.mozilla.org/show_bug.cgi?id=1375585)
            preferHeadless: false,
            usePhantomJS: false,
        },
    } as ConfigOptions);
};
