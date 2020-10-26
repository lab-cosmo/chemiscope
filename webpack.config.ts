import path from 'path';
import * as webpack from 'webpack';

import { execSync } from 'child_process';

const GIT_VERSION = execSync('git describe --tags --dirty').toString().trim();

// eslint-disable-next-line @typescript-eslint/no-explicit-any
type ConfigFn = (env: unknown, argv: any) => webpack.Configuration;
const config: ConfigFn = (env, argv) => {
    if (!('mode' in argv)) {
        throw Error('please specify the build mode');
    }
    return {
        // eslint-disable-next-line @typescript-eslint/no-unsafe-member-access
        devtool: argv.mode === 'development' ? 'inline-source-map' : undefined,
        devServer: {
            injectClient: false,
        },
        entry: {
            chemiscope: './src/index.ts',
            'jsmol-widget': './src/structure/widget.ts',
        },
        module: {
            rules: [
                { test: /\.ts?$/, use: ['ts-loader', './utils/webpack-assert-message.ts'] },
                { test: /\.css?$/, use: ['style-loader', 'css-loader'] },
                { test: /\.html?$/, loader: 'html-loader', options: { minimize: true } },
                { test: /\.svg?$/, loader: 'html-loader', options: { minimize: true } },
                { test: /\.js?$/, use: ['ify-loader'] },
            ],
        },
        output: {
            filename: '[name].min.js',
            library: 'Chemiscope',
            libraryTarget: 'umd',
            path: path.resolve(__dirname, 'dist'),
        },
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
    };
};

export default config;
