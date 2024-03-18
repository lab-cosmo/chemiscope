import assert from 'assert';

/**
 * Component for creating collapsible elements.
 *
 * Manipulated classes and styles are summarized here:
 *  - Hidden:         .collapse
 *  - Shown:          .collapse.show
 *  - Hiding/showing: .collapsing, width/height
 * The .collapse class should be present when calling the Collapse constructor.
 * See [the Bootstrap docs](https://getbootstrap.com/docs/5.2/components/collapse/)
 * for details.
 */
export default class Collapse {
    private _parent: HTMLElement | null = null;
    private _horizontal: boolean;
    private _property: 'width' | 'height';
    private _visible: boolean = false;

    /**
     * Create a Collapse instance.
     * @param _element The element to be made collapsible. It should have the
     *  .collapse class.
     */
    constructor(
        private _element: HTMLElement,
        root: HTMLElement
    ) {
        this._horizontal = this._element.classList.contains('horizontal');
        this._property = this._horizontal ? 'width' : 'height';

        const bsParent = this._element.dataset.bsParent;

        if (bsParent) {
            this._parent = root.querySelector(bsParent);
            assert(this._parent);
        }

        this._element.addEventListener('transitionend', (event) => {
            if (event.target === event.currentTarget && event.propertyName === this._property) {
                this._element.classList.replace('collapsing', 'collapse');

                if (this._visible) {
                    this._element.classList.add('show');
                    this._element.style.removeProperty(this._property);
                }
            }
        });
    }

    hide() {
        if (this._visible) {
            const size = this._element.getBoundingClientRect()[this._property];
            this._element.style.setProperty(this._property, `${size}px`);

            this._element.classList.replace('collapse', 'collapsing');
            this._element.classList.remove('show');

            // eslint-disable-next-line no-unused-expressions
            this._element.offsetHeight;
            this._element.style.removeProperty(this._property);

            this._visible = false;

            if (this._parent && parentsMap.get(this._parent) === this) {
                parentsMap.set(this._parent, null);
            }
        }
    }

    show() {
        if (!this._visible) {
            if (this._parent) {
                // If the parent has a visible collapse, hide it.
                parentsMap.get(this._parent)?.hide();

                // Set the visible collapse of the parent to the current collapse.
                parentsMap.set(this._parent, this);
            }

            const scrollProperty = this._horizontal ? 'scrollWidth' : 'scrollHeight';

            this._element.classList.replace('collapse', 'collapsing');
            this._element.style.setProperty(this._property, '0');

            const size = this._element[scrollProperty];
            this._element.style.setProperty(this._property, `${size}px`);

            this._visible = true;
        }
    }

    toggle(value: boolean = !this._visible) {
        if (value) {
            this.show();
        } else {
            this.hide();
        }
    }

    /**
     * Set an element as a trigger for the collapsible element.
     * @param el The element that should be set as a trigger.
     */
    addTrigger(el: HTMLElement) {
        el.addEventListener('click', (event) => {
            event.preventDefault();
            this.toggle();
        });
    }

    /**
     * Scan the DOM for elements with the [data-bs-toggle] attribute to make
     * them interactive with regards to this component.
     *
     * @param root The root element in which to scan.
     */
    static initialize(root: HTMLElement = document.body) {
        for (const el of root.querySelectorAll('[data-bs-toggle]')) {
            if (el instanceof HTMLElement && 'bsTarget' in el.dataset) {
                const target: HTMLElement | null = root.querySelector(
                    (el.dataset as { bsTarget: string }).bsTarget
                );

                if (target) {
                    getCollapse(target, root).addTrigger(el);
                }
            }
        }
    }
}

/**
 * Return the existing Collapse instance associated with an element, or a new
 * one if none exists. The WeakMap allows elements and instances to be
 * discarded once they are removed from the dom and have no more references.
 */
const collapsesMap = new WeakMap<HTMLElement, Collapse>();
const getCollapse = getFromMap(collapsesMap, (el, root: HTMLElement) => new Collapse(el, root));

function getFromMap<K extends object, V, A extends unknown[]>(
    map: WeakMap<K, V>,
    create: (key: K, ...args: A) => V
) {
    return (key: K, ...args: A): V => {
        if (!map.has(key)) {
            map.set(key, create(key, ...args));
        }

        const value = map.get(key);
        assert(value);

        return value;
    };
}

/**
 * Stores the open collapse for each parent, if any.
 */
const parentsMap = new WeakMap<HTMLElement, Collapse | null>();
