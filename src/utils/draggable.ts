/**
 * @packageDocumentation
 * @module utils
 */

import assert from 'assert';

interface DraggableElement extends HTMLElement {
    dragOffset?: {
        x: number;
        y: number;
    };
}

/** @hidden
 * Make `element` draggable, using `handle_query` as the drag handle.
 * @param  element      DOM element that should be draggable
 * @param  handle_query query (used in element.querySelector(handle_query))
 *                      for the drag handle
 */
export function makeDraggable(element: DraggableElement, handle_query: string): void {
    // adapted from https://stackoverflow.com/a/41737171
    const handle = element.querySelector(handle_query);
    if (handle === null) {
        throw Error(`no element matching '${handle_query}'`);
    }

    // on mouse move, change the position of the element
    const mousemove = (e: MouseEvent) => {
        const offset = element.dragOffset;
        assert(offset !== undefined);
        const left = e.clientX - offset.x;
        const top = e.clientY - offset.y;

        element.style.left = `${left}px`;
        element.style.top = `${top}px`;
    };

    // on mouse up, remove all registered events
    const mouseup = () => {
        document.removeEventListener('mousemove', mousemove);
        document.removeEventListener('mouseup', mouseup);
    };

    // on mouse down, register mouseup & mousemove events
    handle.addEventListener('mousedown', (e: Event) => {
        const event = e as MouseEvent;
        element.dragOffset = {
            x: event.clientX - element.offsetLeft,
            y: event.clientY - element.offsetTop,
        };
        element.style.left = `${element.offsetLeft}px`;
        element.style.top = `${element.offsetTop}px`;
        element.style.width = `${element.offsetWidth}px`;
        element.style.margin = '0';
        element.style.position = 'absolute';

        document.addEventListener('mousemove', mousemove, false);
        document.addEventListener('mouseup', mouseup, false);
    }, false);
}
