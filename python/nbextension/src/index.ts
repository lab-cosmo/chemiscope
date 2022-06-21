/* eslint-disable */

// Entry point for the notebook bundle containing custom model definitions.
//
// Setup notebook base URL
//
// Some static assets may be required by the custom widget javascript. The base
// url for the notebook is not known at build time and is therefore computed
// dynamically.
(window as any).__webpack_public_path__ =
    (document.querySelector('body')!.getAttribute('data-base-url') as string) +
    'nbextensions/chemiscope-widget';

export * from './widget';

// Prevent Chemiscope's Bootstrap v5 from replacing Jupyter's Bootstrap v3.
(window as any).$.fn.modal.noConflict();
