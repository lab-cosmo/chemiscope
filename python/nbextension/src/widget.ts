import { DOMWidgetView } from '@jupyter-widgets/base';
import { addWarningHandler, getByID } from '../../../src/utils';

// Import the CSS
import './widget.css';
import './chemiscope-bootstrap.css';

import { DefaultVisualizer } from '../../../src/index';
import { Dataset } from '../../../src/dataset';

/**
 * The [[ChemiscopeView]] class renders the Chemiscope App as a widget
 * in the Jupyter Notebook output window when instantiated from
 * the Chemiscope Python package.
 */
export class ChemiscopeView extends DOMWidgetView {
    private visualizer?: DefaultVisualizer;

    public render(): void {
        // this function works by first rendering the widget inside `this.el`,
        // and then inserting this.el inside the HTML document.
        const element = this.el;

        // HACK: resize Plotly when inserted. It currently render itself at full
        // width and then needs to be resized. For more on this, see
        // https://github.com/cosmo-epfl/chemiscope/pull/181#discussion_r693005307
        element.addEventListener('DOMNodeInserted', () => {
            window.dispatchEvent(new Event('resize'));
        });

        // handle warnings
        addWarningHandler((message) => {
            const display = getByID('warning-display', element);
            display.style.display = 'block';
            display.getElementsByTagName('p')[0].innerText = message;
        });

        element.innerHTML = `
        <div class="chemiscope-bootstrap">
            <div class="alert alert-warning" role="alert" id="warning-display" style="display: none">
                <button type="button" class="close" onclick="document.getElementById('warning-display').style.display = 'none';">
                    <span aria-hidden="true">&times;</span>
                </button>
                <p></p>
            </div>

            <div id="chemiscope-widget-container">
              <div id="chemiscope-meta-and-map">
                <div id="chemiscope-meta"></div>
                <div id="chemiscope-map" ></div>
              </div>
              <div id="chemiscope-structure-and-info">
                <div id="chemiscope-structure"></div>
                <div id="chemiscope-info"></div>
              </div>
            </div>
          </div>
        </div>`;

        const config = {
            map: element.querySelector('#chemiscope-map') as HTMLElement,
            info: element.querySelector('#chemiscope-info') as HTMLElement,
            meta: element.querySelector('#chemiscope-meta') as HTMLElement,
            structure: element.querySelector('#chemiscope-structure') as HTMLElement,
            maxStructureViewers: 4,
        };

        const data = JSON.parse(this.model.get('data')) as Dataset;
        void DefaultVisualizer.load(config, data).then((visualizer) => {
            this.visualizer = visualizer;
        });
    }

    public remove(): unknown {
        if (this.visualizer !== undefined) {
            this.visualizer.remove();
        }

        return super.remove();
    }
}
