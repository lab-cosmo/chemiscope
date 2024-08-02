/**
 * @packageDocumentation
 * @module utils
 */

/**
 * Represents the log levels available in the logger
 */
type Level = 'info' | 'error' | 'warn';

/**
 * A simple logging utility that supports logging messages at different levels (info, warn, error)
 * and allows adding custom handlers for each log level.
 */
class Logger {
    /**
     * Stores the handlers for each log level. Each log level maps to a set of handler functions
     * that will be invoked when a message is logged at that level.
     */
    private _handlers = new Map<Level, Set<(message: string) => void>>();

    /**
     * Logs an info message
     * @param message - the info message to log
     */
    info(message: string) {
        this._log('info', message);
    }

    /**
     * Logs a warning message
     * @param message - the warning message to log
     */
    warn(message: string) {
        this._log('warn', message);
    }

    /**
     * Logs an error
     * @param error - the error to log
     */
    error(error: unknown) {
        const message = error instanceof Error ? `${error.message}\n${error.stack}` : String(error);
        this._log('error', message);
    }

    /**
     * Logs a message at a specified level and invokes any registered handlers for that level
     * @param level - the level at which to log the message
     * @param message - the message to log
     */
    private _log(level: Level, message: string) {
        // Direct console logging
        // eslint-disable-next-line no-console
        console[level](message);

        // Additional processing of the display if handler was provided
        const handlers = this._handlers.get(level);
        if (handlers) {
            handlers.forEach((callback) => callback(message));
        }
    }

    /**
     * Adds a handler for a specific log level
     * @param level - the log level to add the handler for
     * @param callback - the handler function to call when a log of the specified level occurs
     */
    addHandler(level: Level, callback: (message: string | Error) => void) {
        if (!this._handlers.has(level)) {
            this._handlers.set(level, new Set());
        }
        this._handlers.get(level)!.add(callback);
    }

    /**
     * Removes a handler for a specific log level.
     * @param level - the log level to remove the handler for
     * @param callback - the handler function to remove
     */
    removeHandler(level: Level, callback: (message: string) => void) {
        this._handlers.get(level)?.delete(callback);
    }
}

/**
 * Instance of logging to be used across project
 */
export const logger = new Logger();
