import { getElement } from './utils';
import * as styles from './styles';
import { DisplayTarget } from './indexer';

/**
 * The {@link DisplayTargetToggle} class represents a UI component for switching between two
 * display targets, structures and atom-centered environments, by clicking on the respective buttons.
 *
 * It automatically adjusts the height of the map element to ensure there is space for the toggle.
 * The component uses Shadow DOM to encapsulate styles.
 * A callback function (`onchange`) is triggered whenever the toggle value (structure/atom) is changed
 */
export class DisplayTargetToggle {
    /// Shadow root for isolation
    private _shadow: ShadowRoot;
    /// Toggle buttons container element
    private _toggleContainer: HTMLElement;
    /// Reference to the container element
    private _containerElement: HTMLElement;
    /// Callback fired when the user changes the toggle value
    public onchange: (target: DisplayTarget) => void;

    /**
     * Create a new {@link DisplayTargetToggle} instance
     *
     * @param element HTML element or HTML id of the DOM element where the toggle will be attached
     * @param target display target, either per environements or structures
     */
    constructor(element: string | HTMLElement, target: DisplayTarget) {
        // Create a container element
        this._containerElement = getElement(element);
        const hostElement = document.createElement('div');
        this._containerElement.appendChild(hostElement);

        // Attach to shadow
        this._shadow = hostElement.attachShadow({ mode: 'open' });
        this._shadow.adoptedStyleSheets = [styles.bootstrap, styles.chemiscope];

        // Create a toggle element
        this._toggleContainer = this._createToggleElement(target === 'atom');
        this._shadow.appendChild(this._toggleContainer);

        // Decrease size of map to get space for display target toggle
        this._containerElement.style.setProperty(
            'height',
            `calc(100% - ${this._toggleContainer.offsetHeight}px)`
        );

        // Set up events
        this.onchange = () => {};
    }

    /**
     * Create the HTML structure for the toggle element
     *
     * @param isPerAtom flag indicating if the atom target should be checked
     * @returns the container element of the toggle buttons
     */
    private _createToggleElement(isPerAtom: boolean): HTMLElement {
        const toggleContainer = document.createElement('div');
        toggleContainer.innerHTML = `
            <div class="chsp-target-toggle" title="Toggles between structure and atom-centered data">
                <!-- Spinner -->
                <div id="chsp-target-spinner" class="chsp-target-spinner spinner-border text-secondary" role="status" style="display: none;">
                    <span class="visually-hidden">Loading...</span>
                </div>

                <!-- Buttons -->
                <div class="btn-group-sm" role="group" aria-label="Target toggle">
                    <button type="button" class="btn btn-outline-secondary ${!isPerAtom ? 'active' : ''}" id="structure-btn">structure</button>
                    <button type="button" class="btn btn-outline-secondary ${isPerAtom ? 'active' : ''}" id="atom-btn">atom</button>
                </div>
            </div>
        `;

        // Handle structure button
        const structureBtn = toggleContainer.querySelector('#structure-btn') as HTMLButtonElement;
        structureBtn.onclick = () => this._select('structure');

        // Handle atom button
        const atomBtn = toggleContainer.querySelector('#atom-btn') as HTMLButtonElement;
        atomBtn.onclick = () => this._select('atom');
        return toggleContainer.firstElementChild as HTMLElement;
    }

    /**
     * Handle toggle button click
     * @param target flag indicating if the toggle should be checked
     */
    private _select(target: DisplayTarget): void {
        const isPerAtom = target === 'atom';

        // Activate/desactivate structure button
        const structureBtn = this._toggleContainer.querySelector(
            '#structure-btn'
        ) as HTMLButtonElement;
        structureBtn.classList.toggle('active', !isPerAtom);

        // Activate/desactivate atom button
        const atomBtn = this._toggleContainer.querySelector('#atom-btn') as HTMLButtonElement;
        atomBtn.classList.toggle('active', isPerAtom);

        // Callback
        this.onchange(target);
    }

    /**
     * Toggle loading spinner
     * @param visible flag to toggle
     */
    public loader(visible: boolean): void {
        // Show/hide spinnder
        const spinnerElement = this._toggleContainer.querySelector(
            '#chsp-target-spinner'
        ) as HTMLDivElement;
        spinnerElement.style.display = visible ? 'inline-block' : 'none';

        // Toggle button visibility or disable state
        const buttons = this._toggleContainer.querySelectorAll('.btn-outline-secondary');
        buttons.forEach((button) => {
            if (visible) {
                button.setAttribute('disabled', 'true');
            } else {
                button.removeAttribute('disabled');
            }
        });
    }

    /**
     * Remove HTML added by DisplayTargetToggle in the current document
     */
    public remove(): void {
        this._shadow.host.remove();

        // Reset containerElement height back to 100%
        this._containerElement.style.setProperty('height', '100%');
    }
}
