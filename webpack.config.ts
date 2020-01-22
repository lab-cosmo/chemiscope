import path from "path";
import webpack from "webpack";

const config: webpack.Configuration = {
    mode: "development",
    entry: "./src/index.ts",
    output: {
        path: path.resolve(__dirname, "dist"),
        filename: "chemiscope.min.js",
        libraryTarget: "var",
        library: "Chemiscope",
    },
    resolve: {
        modules: ['./node_modules'],
        extensions: ['.js', '.ts'],
    },
    module: {
        rules: [
            { test: /\.ts?$/, loader: 'ts-loader' },
            { test: /\.css?$/, use: ['style-loader', 'css-loader'] },
            { test: /\.html?$/, loader: 'html-loader', options: { minimize: true } },
            // from https://github.com/plotly/plotly-webpack
            { test: /\.js?$/, use: ['ify-loader', 'transform-loader?plotly.js/tasks/compress_attributes.js'] },
        ],
    },
};

export default config;
