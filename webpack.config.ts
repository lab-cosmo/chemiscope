import path from 'path';
import webpack from 'webpack';

import {execSync} from 'child_process';

const GIT_VERSION = execSync('git describe --tags --dirty').toString().trim();

const config: webpack.Configuration = {
    entry: {
        'chemiscope': './src/index.ts',
        'jsmol-widget': './src/structure/widget.ts',
    },
    mode: 'development',
    module: {
        rules: [
            { test: /\.ts?$/, loader: 'ts-loader' },
            { test: /\.css?$/, use: ['style-loader', 'css-loader'] },
            { test: /\.html?$/, loader: 'html-loader', options: { minimize: true } },
            { test: /\.svg?$/, loader: 'html-loader', options: { minimize: true } },
            { test: /\.js?$/, use: ['ify-loader'] },
        ],
    },
    output: {
        filename: '[name].min.js',
        library: 'Chemiscope',
        libraryTarget: 'var',
        path: path.resolve(__dirname, 'dist'),
    },
    plugins: [
        new webpack.DefinePlugin({
            CHEMISCOPE_GIT_VERSION: `"${GIT_VERSION}"`,
        }),
    ],
    resolve: {
        extensions: ['.js', '.ts'],
        modules: ['./node_modules'],
    },
};

export default config;
