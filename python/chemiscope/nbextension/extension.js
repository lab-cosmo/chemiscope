/* eslint-disable */

// Entry point for the notebook bundle containing custom model definitions.

define(function () {
    'use strict';

    window['requirejs'].config({
        map: {
            '*': {
                'chemiscope-widget': 'nbextensions/chemiscope-widget/chemiscope-widget.min',
            },
        },
    });

    // Export the required load_ipython_extension function
    return {
        load_ipython_extension: function () {},
    };
});
