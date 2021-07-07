// Karma configuration
// Generated on Wed Jul 07 2021 14:58:48 GMT+0200 (Central European Summer Time)

module.exports = function(config) {
  config.set({

    // base path that will be used to resolve all patterns (eg. files, exclude)
    basePath: '',


    // frameworks to use
    // available frameworks: https://www.npmjs.com/search?q=keywords:karma-adapter
    frameworks: ['mocha', 'chai', 'karma-typescript'],


    // list of files / patterns to load in the browser
    files: [
      { pattern: "app/**/*.ts" },
      { pattern: "src/**/*.ts" },
      { pattern: "tests/**/*.test.ts" }
    ],


    // list of files / patterns to exclude
    exclude: [
    ],


    // preprocess matching files before serving them to the browser
    // available preprocessors: https://www.npmjs.com/search?q=keywords:karma-preprocessor
    preprocessors: {
      "src/**/*.ts": ["karma-typescript"]
    },

    // test results reporter to use
    // possible values: 'dots', 'progress'
    // available reporters: https://www.npmjs.com/search?q=keywords:karma-reporter
    reporters: ['progress', 'karma-typescript'],


    // web server port
    port: 9876,


    // enable / disable colors in the output (reporters and logs)
    colors: true,


    // level of logging
    // possible values: config.LOG_DISABLE || config.LOG_ERROR || config.LOG_WARN || config.LOG_INFO || config.LOG_DEBUG
    logLevel: config.LOG_INFO,


    // enable / disable watching file and executing tests whenever any file changes
    autoWatch: true,


    // start these browsers
    // available browser launchers: https://www.npmjs.com/search?q=keywords:karma-launcher
    browsers: ['Chrome'],


    // Continuous Integration mode
    // if true, Karma captures browsers, runs the tests and exits
    singleRun: false,

    // Concurrency level
    // how many browser instances should be started simultaneously
    concurrency: Infinity,

    // Karma Typescript compiler options
    karmaTypescriptConfig: {
      bundlerOptions: {
        resolve: {
          directories: ["src", "node_modules", "src/utils"]
        }
      },
      compilerOptions: {
        module: 'commonjs',
        moduleResolution: 'node',
        paths: {
          'src': ['./src/utils'], 'utils/*': ['./src/utils/*']
        },
        baseUrl: "./"
      },
      tsconfig: './tsconfig.json',
    },
  })
}
