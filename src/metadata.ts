/**
 * @packageDocumentation
 * @module main
 */

import markdown from 'markdown-it';

import { Metadata } from './dataset';
import { getElement } from './utils';
import * as styles from './styles';
import Modal from './modal';

import INFO_SVG from './static/info.svg';

function generateModal(metadata: Metadata): string {
    const md = markdown({
        html: false,
        linkify: true,
        typographer: true,
    });

    // Deal with missing metadata
    let description = 'No description for this dataset.';
    if (metadata.description !== undefined) {
        description = metadata.description;
    }

    let authors = 'No authors for this dataset.';
    if (metadata.authors !== undefined) {
        authors = '<ul class="chsp-authors-list">';
        for (const author of metadata.authors) {
            // remove enclosing <p> tag added by markdown rendering
            const render = md.render(author).slice(3, -5);
            authors += `<li>${render}</li>`;
        }
        authors += '</ul>';
    }

    let ref = 'No references for this dataset.';
    if (metadata.references !== undefined) {
        ref = '<ul>';
        for (const reference of metadata.references) {
            ref += `<li>${md.render(reference)}</li>`;
        }
        ref += '</ul>';
    }

    return `<div class="modal chemiscope-modal fade" tabindex='-1'>
        <div class="modal-dialog modal-lg">
            <div class="modal-content">
                <div class="modal-header chsp-modal-header">
                    <h4 class="modal-title">${metadata.name}</h4>
                    <button type="button" class="btn-close" data-bs-dismiss="modal" aria-label="Close"></button>
                </div>
                <div class="modal-body">
                    <div>${md.render(description)}</div>
                    <h5 style="margin-top: 1em;">Authors</h5>
                    ${authors}
                    <h5 style="margin-top: 1em;">References</h5>
                    ${ref}
                </div>
            </div>
        </div>
    </div>`;
}

/**
 * The [[MetadataPanel]] class displays information regarding the dataset: name,
 * authors, references, description, etc.
 *
 * By default, only the name is shown, and clicking on it reveals a modal with
 * the other information.
 */
export class MetadataPanel {
    /// Shadow root for isolation
    private _shadow: ShadowRoot;
    /// The HTML element serving as root element for the name
    private _name: HTMLElement;
    /// The HTML element containing the modal for the current dataset
    private _modal: Modal;

    /**
     * Create a new [[MetadataPanel]] inside the HTML element with the given
     * `id`.
     *
     * @param element HTML element or HTML id of the DOM element where the name of the dataset should be inserted
     * @param metadata dataset metadata
     */
    constructor(element: string | HTMLElement, metadata: Metadata) {
        // sanitize HTML, all the other field will go through markdown
        metadata.name = metadata.name.replace(/</g, '&lt;');
        metadata.name = metadata.name.replace(/>/g, '&gt;');

        const containerElement = getElement(element);
        const hostElement = document.createElement('div');

        hostElement.style.setProperty('height', '100%');
        containerElement.appendChild(hostElement);

        this._shadow = hostElement.attachShadow({ mode: 'open' });

        this._shadow.adoptedStyleSheets = [styles.bootstrap, styles.chemiscope];

        this._name = document.createElement('div');

        // this._root = document.createElement('div');
        // this._root.style.setProperty('height', '100%');
        this._shadow.appendChild(this._name);

        this._name.innerHTML = `<div class="chsp-meta">
            <span> ${metadata.name} </span>
            <div class='chsp-info-icon'>${INFO_SVG}</div>
        </div>`;

        const template = document.createElement('template');
        template.innerHTML = generateModal(metadata);
        const modalElement = template.content.firstChild as HTMLElement;

        this._modal = new Modal(modalElement);
        this._modal.shadow.adoptedStyleSheets = [styles.bootstrap, styles.chemiscope];

        this._name.onclick = () => this._modal.open();

        // Stop propagation of keydown events. This is required for the Jupyter
        // integration, otherwise Jupyter tries to interpret key press in the
        // modal as its own input
        modalElement.addEventListener('keydown', (event) => {
            event.stopPropagation();
        });
    }

    /**
     * Remove all HTML added by this [[EnvironmentInfo]] in the current document
     */
    public remove(): void {
        this._modal.remove();
        this._shadow.host.remove();
    }
}
