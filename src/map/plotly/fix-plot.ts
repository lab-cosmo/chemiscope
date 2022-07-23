import { PlotlyScatterElement } from './plotly-scatter';

/**
 * Fix a Plotly bug that causes incorrect mouse selection coordinates when
 * creating 2d plots inside shadow roots. The function adds mouse event
 * listeners to the plot and document and re-emits them with correct
 * coordinates.

 * Because some listeners are added to the document, they must be explicitly
 * removed (using the returned object) when the plot is destroyed. The fix is
 * automatically disabled when changing the plot to 3d, and enabled when
 * initializing or going back to 2d.
 *
 * See [Plotly.js issue #6108](https://github.com/plotly/plotly.js/issues/6108)
 * for details.
 *
 * @param plot The element to be fixed.
 * @returns A { disable, enable } object to toggle the fix.
 */
export default function fixPlot(plot: PlotlyScatterElement) {
    let movement: {
        data: {
            buttons: number;
            clientX: number;
            clientY: number;
            ctrlKey: boolean;
        };
        target: EventTarget;
    } | null = null;

    const mousedownListener = (event: MouseEvent | TouchEvent) => {
        if (event.isTrusted && event.target) {
            event.preventDefault();
            event.stopPropagation();

            // Store the event data. The event will be re-remited later on
            // as we cannot yet determine whether the coordinates need to be
            // corrected (if the user selected a point) or not (if the user
            // selected an area for zooming).
            movement = {
                data:
                    event instanceof MouseEvent
                        ? {
                              buttons: event.buttons,
                              clientX: event.clientX,
                              clientY: event.clientY,
                              ctrlKey: event.ctrlKey,
                          }
                        : {
                              buttons: 1,
                              clientX: event.touches[0].clientX,
                              clientY: event.touches[0].clientY,
                              ctrlKey: event.ctrlKey,
                          },
                target: event.target,
            };
        }
    };

    const mouseupListener = (event: MouseEvent | TouchEvent) => {
        if (movement && event.target) {
            event.preventDefault();
            event.stopImmediatePropagation();

            // 'composed' is necessary for the event to behave correctly
            // in a shadow root.
            const downEvent = new Event('mousedown', { composed: true });
            Object.assign(downEvent, {
                ...movement.data,
                // The correction is applied here and was determined empirically.
                clientX: movement.data.clientX - 50,
                clientY: movement.data.clientY - 50,
            });

            const target = movement.target;
            movement = null;
            target.dispatchEvent(downEvent);

            const upEvent = new Event('mouseup', { bubbles: true, composed: true });
            event.target.dispatchEvent(upEvent);
        }
    };

    const mousemoveListener = (event: MouseEvent | TouchEvent) => {
        if (movement) {
            event.preventDefault();
            event.stopImmediatePropagation();

            const ev = event instanceof MouseEvent ? event : event.changedTouches[0];

            // Additional logic to allow area selection, in which case the bug
            // doesn't happen and no fix is needed.
            if (
                Math.abs(ev.clientX - movement.data.clientX) > 5 ||
                Math.abs(ev.clientY - movement.data.clientY) > 5
            ) {
                const downEvent = new Event('mousedown', { composed: true });
                Object.assign(downEvent, movement.data);

                const target = movement.target;
                movement = null;
                target.dispatchEvent(downEvent);
            }
        }
    };

    let enabled = false;

    const enable = () => {
        enabled = true;

        plot.addEventListener('mousedown', mousedownListener, { capture: true });
        document.addEventListener('mouseup', mouseupListener);
        document.addEventListener('mousemove', mousemoveListener);

        plot.addEventListener('touchstart', mousedownListener, { capture: true });
        document.addEventListener('touchmove', mousemoveListener, { passive: false });
        document.addEventListener('touchend', mouseupListener, { passive: false });
    };

    const disable = () => {
        enabled = false;

        plot.removeEventListener('mousedown', mousedownListener, { capture: true });
        document.removeEventListener('mouseup', mouseupListener);
        document.removeEventListener('mousemove', mousemoveListener);

        plot.removeEventListener('touchstart', mousedownListener, { capture: true });
        document.removeEventListener('touchmove', mousemoveListener);
        document.removeEventListener('touchend', mouseupListener);
    };

    // Toggle the fix when the plot type changes.
    plot.on('plotly_restyle', ([update]) => {
        if ('type' in update) {
            if (enabled && update.type === 'scatter3d') {
                disable();
            }
            if (!enabled && update.type === 'scattergl') {
                enable();
            }
        }
    });

    // Enable the fix if the plot is 2d.
    if (plot._fullData[0].type === 'scattergl') {
        enable();
    }

    return { disable, enable };
}
