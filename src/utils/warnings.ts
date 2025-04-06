/**
 * @packageDocumentation
 * @module utils
 */

/** A {@link WarningHandler} is called when a new warning is emitted */
export type WarningHandler = (message: string) => void;
export type WarningSource = object | null;

/** A minimalistic class to handle warnings */
export class Warnings {
    private handlersList: WarningHandler[] = [
        // eslint-disable-next-line no-console
        (message) => console.warn(message),
    ];

    /** Adds a warning handler to the list */
    public add(handler: WarningHandler) : void {
        this.handlersList.push(handler);
    }

    /** Sends a warning message though all handlers */
    public send(message: string) : void{
        for (const cb of this.handlersList) {
            console.log(this.handlersList.length, ': using handler ', cb)
            cb(message);
        }
    }
}


/** List of registered warnings handlers */
const WARNINGS_HANDLERS: [WarningHandler, WarningSource][] = [
    // eslint-disable-next-line no-console
    [(message) => console.warn(message), null],
];

/** @hidden
 * Send a warning to the user with the given `message`
 */
export function sendWarning(message: string, source: WarningSource=null): void {
    for (const cb of WARNINGS_HANDLERS) {
        if (source === null || cb[1] === null || source === cb[1]) {
            cb[0](message);
        }
    }
}

/** Register the `handler` function to be called when a warning is emitted */
export function addWarningHandler(handler: WarningHandler, source: WarningSource=null): void {
    WARNINGS_HANDLERS.push([handler, source]);
}

/** Remove warning `handler`s matching a given source */
export function removeWarningHandlers(source: WarningSource): void {
    console.log("removing warning handlers for ", source);
    for (let i = WARNINGS_HANDLERS.length - 1; i >= 0; i--) {
        const [_, handlerSource] = WARNINGS_HANDLERS[i];
        if (handlerSource === source) {
            WARNINGS_HANDLERS.splice(i, 1);
        }
    }
}

