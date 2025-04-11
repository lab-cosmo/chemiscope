/**
 * @packageDocumentation
 * @module utils
 */

/** A {@link WarningHandler} is called when a new warning is emitted */
export type WarningHandler = (message: string, timeout?: number) => void;

/** A minimalistic class to handle warnings */
export class Warnings {
    // set to zero for no timeout, to a negative value to
    // disable warnings. handlers can use the value given here
    // to hide the warnings after a given time
    public defaultTimeout: number = 0;

    private handlersList: WarningHandler[] = [
        // eslint-disable-next-line no-console
        (message) => console.warn(message),
    ];

    /** Adds a warning handler to the list */
    public addHandler(handler: WarningHandler): void {
        this.handlersList.push(handler);
    }

    /** Sends a warning message though all handlers */
    public sendMessage(message: string, userTimeout: number = 0): void {
        const timeout = userTimeout || this.defaultTimeout;
        if (timeout < 0) {
            return;
        }

        for (const cb of this.handlersList) {
            cb(message, timeout);
        }
    }
}
