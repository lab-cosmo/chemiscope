/**
 * @packageDocumentation
 * @module utils
 */

/** A [[WarningHandler]] is called when a new warning is emitted */
export type WarningHandler = (message: string) => void;

/** List of registered warnings handlers */
const WARNINGS_HANDLERS: WarningHandler[] = [
    (message) => console.warn(message),
];

/** @hidden
 * Send a warning to the user with the given `message`
 */
export function sendWarning(message: string) {
    for (const cb of WARNINGS_HANDLERS) {
        cb(message);
    }
}

/** Register the `handler` function to be called when a warning is emitted */
export function addWarningHandler(handler: WarningHandler) {
    WARNINGS_HANDLERS.push(handler)
}
