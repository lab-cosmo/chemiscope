import { WEBPACK_CONFIG } from '../webpack.config';
import { Config, ConfigOptions } from 'karma';

module.exports = (config: Config) => {
    config.set({
        exclude: [],
        files: ['./**/*.test.ts'],
        frameworks: ['webpack', 'mocha', 'detectBrowsers'],

        preprocessors: {
            './**/*.test.ts': ['webpack'],
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
                return availableBrowsers;
            },
            preferHeadless: false,
            usePhantomJS: false,
        },
    } as ConfigOptions);
};
