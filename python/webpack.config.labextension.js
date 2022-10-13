/* eslint-disable */

const webpack = require('webpack');
const execSync = require('child_process').execSync;

const GIT_VERSION = execSync('git describe --tags --dirty').toString().trim();

module.exports = {
    plugins: [
        new webpack.DefinePlugin({
            CHEMISCOPE_GIT_VERSION: `"${GIT_VERSION}"`,
        }),
    ],
    resolve: {
        extensions: ['.js', '.ts'],
    },
    module: {
        rules: [
            { test: /\.ts$/, use: ['ts-loader', './utils/webpack-assert-message.js'] },
            { test: /\.html\.in$/, loader: 'raw-loader' },
            { test: /\.svg$/, loader: 'raw-loader' },
            // this is required by plotly, since we are building our own bundle
            { test: /\.js$/, use: ['ify-loader'] },
        ],
    },
};
