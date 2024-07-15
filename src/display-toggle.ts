import { getElement } from './utils';
import * as styles from './styles';

export class DisplayToggle {
    /// Shadow root for isolation
    private _shadow: ShadowRoot;
    /// Toggle buttons container element
    private _toggleContainer: HTMLElement;
    /// Callback fired when the user changes the toggle value
    public onchange: (checked: boolean) => void;

    /**
     * Create a new {@link DisplayToggle} instance
     *
     * @param element HTML element or HTML id of the DOM element where the toggle will be attached
     * @param toggled flag indicating if the toggle should be checked
     * @param disabled flag indicating if the toggle should be disabled
     */
    constructor(
        element: string | HTMLElement,
        toggled: boolean = false,
        disabled: boolean = false
    ) {
        // Create a container element
        const containerElement = getElement(element);
        const hostElement = document.createElement('div');
        containerElement.appendChild(hostElement);

        // Attach to shadow
        this._shadow = hostElement.attachShadow({ mode: 'open' });
        this._shadow.adoptedStyleSheets = [styles.bootstrap, styles.chemiscope];

        // Create a toggle element
        this._toggleContainer = this._createToggleElement(toggled, disabled);
        this._shadow.appendChild(this._toggleContainer);

        // Set up events
        this.onchange = () => {};
    }

    /**
     * Create the HTML structure for the toggle element
     *
     * @param toggled flag indicating if the toggle should be checked
     * @param disabled flag indicating if the toggle should be disabled
     * @returns the container element of the toggle buttons
     */
    private _createToggleElement(toggled: boolean, disabled: boolean): HTMLElement {
        const toggleContainer = document.createElement('div');
        toggleContainer.innerHTML = `
            <div class="chsp-mode-toggle">
                <div class="btn-group-sm" role="group" aria-label="Mode toggle" style="display: ${disabled ? 'none' : 'block'}">
                    <button type="button" class="btn btn-outline-secondary ${toggled ? 'active' : ''}" id="structure-btn">Structures</button>
                    <button type="button" class="btn btn-outline-secondary ${!toggled ? 'active' : ''}" id="atom-btn">Atoms</button>
                </div>

                <!-- Loading Spinner -->
                <div id="loading-spinner" class="chsp-mode-spinner spinner-border text-secondary" role="status" style="display: none;">
                    <span class="visually-hidden">Loading...</span>
                </div>
            </div>
        `;

        // Handle structure button
        const structureBtn = toggleContainer.querySelector('#structure-btn') as HTMLButtonElement;
        structureBtn.onclick = () => this._handleToggle(true);

        // Handle atom button
        const atomBtn = toggleContainer.querySelector('#atom-btn') as HTMLButtonElement;
        atomBtn.onclick = () => this._handleToggle(false);
        return toggleContainer.firstElementChild as HTMLElement;
    }

    /**
     * Handle toggle button click
     * @param toggled flag indicating if the toggle should be checked
     */
    private _handleToggle(toggled: boolean): void {
        // Structure button
        const structureBtn = this._toggleContainer.querySelector(
            '#structure-btn'
        ) as HTMLButtonElement;
        structureBtn.classList.toggle('active', toggled);

        // Atom button
        const atomBtn = this._toggleContainer.querySelector('#atom-btn') as HTMLButtonElement;
        atomBtn.classList.toggle('active', !toggled);

        // Callback
        this.onchange(toggled);
    }

    /**
     * Remove HTML added by DisplayToggle in the current document
     */
    remove(): void {
        this._shadow.host.remove();
    }

    /**
     * Toggle loading spinner
     * @param visible flag to toggle
     */
    loader(visible: boolean): void {
        const loadingSpinner = this._toggleContainer.querySelector(
            '#loading-spinner'
        ) as HTMLDivElement;

        // Toggle show/hide elements
        loadingSpinner.style.display = visible ? 'inline-block' : 'none';
    }
}
