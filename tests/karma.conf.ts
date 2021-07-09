import { BASE_CONFIG } from '../webpack.config';
import path from 'path';

module.exports = (config: any) => {
    config.set({
        autoWatch: false,
        browserNoActivityTimeout: 10000,
        concurrency: 1,
        client: {
            mocha: {
                timeout: 8000,
            },
        },
        exclude: [],
        files: ['./**/*.test.ts'],
        frameworks: ['webpack', 'mocha', 'detectBrowsers'],

        preprocessors: {
            './**/*.test.ts': ['webpack'],
        },
        reporters: ['progress'],
        singleRun: true,

        webpack: BASE_CONFIG,

        webpackMiddleware: {
            stats: 'errors-only',
        },

        detectBrowsers: {
            postDetection: function (availableBrowsers: any) {
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
    });
};
