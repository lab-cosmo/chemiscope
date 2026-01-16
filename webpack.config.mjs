import path from 'path';
import webpack from 'webpack';

import { execSync } from 'child_process';

const GIT_VERSION = execSync('git describe --tags --dirty').toString().trim();

export const WEBPACK_CONFIG = {
    plugins: [
        new webpack.DefinePlugin({
            CHEMISCOPE_GIT_VERSION: `"${GIT_VERSION}"`,
        }),
        new webpack.ProvidePlugin({
            process: 'process/browser',
        }),
    ],
    resolve: {
        extensions: ['.js', '.ts'],
        modules: ['./node_modules'],
    },
    module: {
        rules: [
            { test: /\.ts$/, use: ['ts-loader', './utils/webpack-assert-message.js'] },
            { test: /\.css$/, use: ['style-loader', 'css-loader'] },
            { test: /\.less$/, use: ['style-loader', 'css-loader', 'less-loader'] },
            { test: /\.html\.in$/, loader: 'raw-loader' },
            { test: /\.svg$/, loader: 'raw-loader' },
            // this is required by plotly, since we are building our own bundle
            { test: /\.js$/, use: ['ify-loader'] },
        ],
    },
};

const config = (env, argv) => {
    const mode = argv.mode || process.env.CHEMISCOPE_WEBPACK_MODE || 'production';
    return {
        ...WEBPACK_CONFIG,
        mode: mode,
        target: 'web',
        // eslint-disable-next-line @typescript-eslint/no-unsafe-member-access
        devtool: mode === 'development' ? 'inline-source-map' : undefined,
        entry: {
            chemiscope: './src/index.ts',
            'chemiscope-app': './app/app.ts',
        },
        output: {
            filename: '[name].min.js',
            library: 'Chemiscope',
            // Use UMD modules for the main entry points
            libraryTarget: 'umd',
            path: path.resolve('./dist'),
        },
        // we know we have a very large `chemiscope.min.js`, we don't need
        // webpack to tell us everytime.
        performance: {
            hints: false,
            maxEntrypointSize: 512000,
            maxAssetSize: 512000,
        },
    };
};

export default config;
