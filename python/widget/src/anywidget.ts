import type { RenderProps } from '@anywidget/types';
import * as jquery from 'jquery';

import { ChemiscopeView, MapView, StructureView } from './widget';

// Some hosts/plugins expect `window.$` to exist; ESM imports do not set it.
(window as unknown as { $: unknown }).$ = jquery;

/**
 * anywidget entry point. A single ES module is loaded for all three widget modes
 * (default, structure, map); the `mode` trait synced from Python selects which
 * view to instantiate.
 */
function render({ model, el }: RenderProps): () => void {
    const mode = model.get('mode') as string;

    let view: ChemiscopeView | StructureView | MapView;
    if (mode === 'structure') {
        view = new StructureView(model, el);
    } else if (mode === 'map') {
        view = new MapView(model, el);
    } else {
        view = new ChemiscopeView(model, el);
    }

    view.render();

    // anywidget calls the returned function to clean up when the view is removed
    return () => view.dispose();
}

export default { render };
