/**
 * @packageDocumentation
 * @module utils
 */

/** A {@link WarningHandler} is called when a new warning is emitted */
export type WarningHandler = (message: string) => void;

export type ColorModeMessages = {
    [key: string]: {
        allValuesNaN: string;
        someValuesNaN: string;
    };
};

/** List of registered warnings handlers */
const WARNINGS_HANDLERS: WarningHandler[] = [
    // eslint-disable-next-line no-console
    (message) => console.warn(message),
];

/** @hidden
 * Send a warning to the user with the given `message`
 */
export function sendWarning(message: string): void {
    for (const cb of WARNINGS_HANDLERS) {
        cb(message);
    }
}

/** Register the `handler` function to be called when a warning is emitted */
export function addWarningHandler(handler: WarningHandler): void {
    WARNINGS_HANDLERS.push(handler);
}
