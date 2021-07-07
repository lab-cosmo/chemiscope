/**
 * @packageDocumentation
 * @module main
 */

import assert from 'assert';
import markdown from 'markdown-it';

import { Metadata } from './dataset';
import { generateGUID, getElement } from './utils';

import INFO_SVG from './static/info.svg';

function generateName(guid: string, name: string): string {
    return `<div data-toggle="modal" data-target="#${guid}">
    <span> ${name} </span> 
    <div class='chsp-info-icon'>${INFO_SVG}</div>
    </div>`;
}

function generateModal(guid: string, metadata: Metadata): string {
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

    return `<div id=${guid} class="modal fade" tabindex='-1'>
        <div class="modal-dialog modal-lg">
            <div class="modal-content">
                <div class="modal-header chsp-modal-header">
                    <h4 class="modal-title">${metadata.name}</h4>
                    <button type="button" class="close" data-dismiss="modal" aria-hidden="true">&times;</button>
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
    /// Global id of this panel
    private _guid: string;
    /// The HTML element serving as root element for the name
    private _name: HTMLElement;
    /// The HTML element containing the modal for the current dataset
    private _modal: HTMLElement;

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

        this._guid = `chsp-${generateGUID()}`;

        this._name = getElement(element);

        this._name.innerHTML = generateName(this._guid, metadata.name);
        this._name.classList.add('chsp-meta');

        const template = document.createElement('template');
        template.innerHTML = generateModal(this._guid, metadata);
        this._modal = template.content.firstChild as HTMLElement;
        document.body.appendChild(this._modal);
    }

    /**
     * Remove all HTML added by this [[EnvironmentInfo]] in the current document
     */
    public remove(): void {
        this._name.innerHTML = '';
        if (this._modal.classList.contains('show')) {
            const close = this._modal.querySelector('.close');
            assert(close !== null);
            (close as HTMLElement).click();
        }
        this._modal.remove();
    }
}
