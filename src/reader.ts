/**
 * @packageDocumentation
 * @module main
 */

import assert from 'assert';

import { Dataset, Structure, UserStructure } from './dataset';

// Reading large datasets into a single string and calling `JSON.parse` can
// exhaust memory or exceed the maximal length of strings (~512 MiB). Instead,
// `readDataset` reads the file incrementally (decompressing it on the fly if
// needed), and parses each structure and each property separately.

// drop the consumed part of the parse buffer past this size
const COMPACT_THRESHOLD = 1 << 16;
// progress updates are throttled to keep the displayed numbers readable
const PROGRESS_THROTTLE_MS = 100;
// how often reading structures lets the browser render; each pause takes about
// a frame, so pausing too often slows down loading
const YIELD_INTERVAL_MS = 250;

/** @hidden Progress of {@link readDataset} */
export interface ReadProgress {
    /** number of bytes read from the file (compressed bytes for gzip files) */
    bytesRead: number;
    /** size of the file in bytes */
    bytesTotal: number;
    /** number of structures read so far */
    structures: number;
    /** key of the section of the dataset being read, e.g. `structures` */
    section: string;
}

/** @hidden Options for {@link readDataset} */
export interface ReadOptions {
    /**
     * What to keep in the dataset for each structure, given its raw JSON text,
     * its number of atoms and its index, e.g. to store structures outside of
     * memory. If this returns a promise, reading waits for it. Without this
     * callback, all structures are parsed together and kept in the dataset.
     */
    structure?: (
        raw: string,
        size: number,
        index: number
    ) => Structure | UserStructure | Promise<Structure | UserStructure>;
    /**
     * Called regularly with the reading progress. If this returns a promise,
     * reading waits for it, e.g. to let the page render the progress.
     */
    onProgress?: (progress: ReadProgress) => void | Promise<void>;
}

/**
 * @hidden
 * Read a chemiscope dataset from a JSON file, which can be gzip compressed.
 * Values can be `NaN`, `Infinity` or `-Infinity`, as written by Python.
 */
export async function readDataset(file: Blob, options: ReadOptions = {}): Promise<Dataset> {
    // count bytes as they are read, for the progress
    let bytesRead = 0;
    const counter = new TransformStream<Uint8Array, Uint8Array>({
        transform(chunk, controller) {
            bytesRead += chunk.byteLength;
            controller.enqueue(chunk);
        },
    });

    let stream = file.stream().pipeThrough(counter);
    if (await isGzipped(file)) {
        if (typeof DecompressionStream === 'undefined') {
            throw new Error('reading gzip files requires DecompressionStream, unavailable here');
        }
        const gunzip = new DecompressionStream('gzip') as unknown as ReadableWritablePair<
            Uint8Array,
            Uint8Array
        >;
        stream = stream.pipeThrough(gunzip);
    }
    const parser = new StreamingJSONParser(stream.getReader());

    const structures: (Structure | UserStructure)[] = [];
    let section = '';
    let lastProgress = 0;
    const reportProgress = async (force: boolean) => {
        const now = Date.now();
        if (
            options.onProgress !== undefined &&
            (force || now - lastProgress > PROGRESS_THROTTLE_MS)
        ) {
            lastProgress = now;
            const progress = { bytesRead, bytesTotal: file.size, structures: structures.length };
            await options.onProgress({ ...progress, section });
        }
    };

    let lastYield = Date.now();
    const yieldRegularly = (): Promise<void> | undefined => {
        if (Date.now() - lastYield < YIELD_INTERVAL_MS) {
            return undefined;
        }
        return (async () => {
            await reportProgress(false);
            await new Promise<void>((resolve) => setTimeout(resolve, 0));
            lastYield = Date.now();
        })();
    };

    const addStructure = (raw: string): Promise<void> | undefined => {
        assert(options.structure !== undefined);
        const index = structures.length;
        const structure = options.structure(raw, structureSize(raw, index), index);

        if (structure instanceof Promise) {
            return structure.then((s) => {
                structures.push(s);
                return yieldRegularly();
            });
        }
        structures.push(structure);
        return yieldRegularly();
    };

    const dataset: Record<string, unknown> = {};
    await parser.readObject(async (key) => {
        section = key;
        await reportProgress(true);

        if (key === 'structures' && options.structure !== undefined) {
            await parser.readArray(addStructure);
            dataset.structures = structures;
        } else if (key === 'properties') {
            // each property is parsed separately, to stay below the limits of strings
            const properties: Record<string, unknown> = {};
            await parser.readObject(async (name) => {
                properties[name] = parseJsonWithNaN(await parser.readValue());
            });
            dataset.properties = properties;
        } else {
            dataset[key] = parseJsonWithNaN(await parser.readValue());
        }
    });

    return dataset as unknown as Dataset;
}

// non-finite values written by Python, which are not valid JSON
const NON_FINITE = new Map([
    ['***NaN***', NaN],
    ['***Infinity***', Infinity],
    ['***-Infinity***', -Infinity],
]);

/**
 * @hidden
 * Parse JSON `text`, which can contain `NaN`, `Infinity` and `-Infinity` as
 * written by Python.
 */
export function parseJsonWithNaN(text: string): unknown {
    // strings can contain "NaN" (e.g. the chemical formula "NaNO3"), so only use
    // the much slower parsing below if the text is not valid JSON
    try {
        return JSON.parse(text) as unknown;
    } catch (error) {
        if (!text.includes('NaN') && !text.includes('Infinity')) {
            throw error;
        }
    }

    // replace non-finite values outside of strings by placeholder strings
    let out = '';
    let chunkStart = 0;
    for (let i = 0; i < text.length; i++) {
        const char = text[i];
        if (char === '"') {
            const end = makeValueScanner()(text, i);
            i = (end < 0 ? text.length : end) - 1;
        } else if (char === 'N' || char === 'I' || char === '-') {
            for (const value of ['NaN', 'Infinity', '-Infinity']) {
                if (text.startsWith(value, i)) {
                    out += text.slice(chunkStart, i) + `"***${value}***"`;
                    i += value.length - 1;
                    chunkStart = i + 1;
                    break;
                }
            }
        }
    }
    out += text.slice(chunkStart);

    return JSON.parse(out, (_, value: unknown) =>
        typeof value === 'string' && NON_FINITE.has(value) ? NON_FINITE.get(value) : value
    ) as unknown;
}

/** Check if `file` starts with the gzip magic number */
async function isGzipped(file: Blob): Promise<boolean> {
    const header = new Uint8Array(await file.slice(0, 2).arrayBuffer());
    return header[0] === 0x1f && header[1] === 0x8b;
}

// `"size"` is usually the first key of a structure
const SIZE_FIRST = /^\s*\{\s*"size"\s*:\s*(\d+)\s*[,}]/;

/** Get the number of atoms from the raw JSON text of a structure */
function structureSize(raw: string, index: number): number {
    const match = SIZE_FIRST.exec(raw);
    const size = match !== null ? Number(match[1]) : (parseJsonWithNaN(raw) as Structure).size;
    if (!Number.isInteger(size)) {
        throw new Error(`structure ${index} is missing an integer "size"`);
    }
    return size;
}

const QUOTE = 0x22;
const BACKSLASH = 0x5c;
const OPEN_BRACE = 0x7b;
const OPEN_BRACKET = 0x5b;

function isWhitespace(char: number): boolean {
    return char === 0x20 || char === 0x0a || char === 0x0d || char === 0x09;
}

// characters that change the state of the scanner, outside and inside of strings
const STRUCTURAL_CHARS = /["[\]{}]/g;
const STRING_CHARS = /["\\]/g;

/**
 * Create a function finding the end of the JSON object, array or string
 * starting at index `from` of `buf`, returning the index just past its end, or
 * -1 if the value is incomplete. The function keeps track of where it stopped,
 * and can be called again with the next part of the value.
 */
function makeValueScanner(): (buf: string, from: number) => number {
    let depth = 0;
    let inString = false;
    let escape = false;

    // jump directly to the next character that matters, which is much faster
    // than looking at every character of large values
    return function scan(buf: string, from: number): number {
        let i = from;
        if (escape) {
            if (i >= buf.length) {
                return -1;
            }
            // skip the escaped character
            escape = false;
            i++;
        }

        while (true) {
            const search = inString ? STRING_CHARS : STRUCTURAL_CHARS;
            search.lastIndex = i;
            if (!search.test(buf)) {
                return -1;
            }
            i = search.lastIndex - 1;
            const char = buf.charCodeAt(i);

            if (inString) {
                if (char === BACKSLASH) {
                    if (i + 1 >= buf.length) {
                        escape = true;
                        return -1;
                    }
                    i++;
                } else {
                    inString = false;
                    if (depth === 0) {
                        return i + 1;
                    }
                }
            } else if (char === QUOTE) {
                inString = true;
            } else if (char === OPEN_BRACE || char === OPEN_BRACKET) {
                depth++;
            } else {
                depth--;
                if (depth === 0) {
                    return i + 1;
                }
            }
            i++;
        }
    };
}

/** Incremental JSON tokenizer over a stream of bytes */
class StreamingJSONParser {
    private buf = '';
    private pos = 0;
    private eof = false;
    private decoder = new TextDecoder('utf-8');

    constructor(private reader: ReadableStreamDefaultReader<Uint8Array>) {}

    /**
     * Read the object at the current position, calling `onEntry` with each
     * key. `onEntry` must read the corresponding value.
     */
    public async readObject(onEntry: (key: string) => Promise<void>): Promise<void> {
        await this.expect('{', 'an object');
        while (true) {
            const char = await this.peek();
            if (char === '}') {
                this.pos++;
                return;
            } else if (char === ',') {
                this.pos++;
                continue;
            } else if (char !== '"') {
                throw new Error(`expected an object key, found ${JSON.stringify(char)}`);
            }

            const key = JSON.parse(await this.readValue()) as string;
            await this.expect(':', '":" after an object key');
            await onEntry(key);
        }
    }

    /**
     * Read the array at the current position, calling `onElement` with the raw
     * JSON text of each element. If `onElement` returns a promise, it is
     * awaited before reading the next element.
     */
    public async readArray(onElement: (raw: string) => Promise<void> | undefined): Promise<void> {
        await this.expect('[', 'an array');
        while (true) {
            // elements already in the buffer are read synchronously, since
            // awaiting for each one is slow for arrays with many small elements
            const char = this.skipWhitespace() ?? (await this.peek());
            if (char === ']') {
                this.pos++;
                return;
            } else if (char === ',') {
                this.pos++;
                continue;
            }

            const raw = this.readBufferedValue() ?? (await this.readValue());
            const pending = onElement(raw);
            if (pending !== undefined) {
                await pending;
            }
        }
    }

    /** Read the raw JSON text of the object, array or string at the current position */
    public async readValue(): Promise<string> {
        const first = await this.peek();
        if (first !== '{' && first !== '[' && first !== '"') {
            throw new Error(`expected an object, array or string, found ${JSON.stringify(first)}`);
        }

        // large values span many chunks of the stream, keep the parts read so
        // far separately instead of growing the buffer
        const scan = makeValueScanner();
        const parts: string[] = [];
        let start = this.pos;
        while (true) {
            const end = scan(this.buf, this.pos);
            if (end >= 0) {
                parts.push(this.buf.slice(start, end));
                this.pos = end;
                return parts.join('');
            }

            parts.push(this.buf.slice(start));
            this.buf = '';
            this.pos = 0;
            start = 0;
            if (!(await this.pull())) {
                throw new Error('unexpected end of JSON data');
            }
        }
    }

    /** Return the next non-whitespace character, without consuming it */
    private async peek(): Promise<string> {
        if (this.pos > COMPACT_THRESHOLD) {
            // drop the consumed part of the buffer
            this.buf = this.buf.slice(this.pos);
            this.pos = 0;
        }

        while (true) {
            const char = this.skipWhitespace();
            if (char !== undefined) {
                return char;
            }
            if (!(await this.pull())) {
                throw new Error('unexpected end of JSON data');
            }
        }
    }

    /** Consume the next non-whitespace character, which should be `char` */
    private async expect(char: string, description: string): Promise<void> {
        const next = await this.peek();
        if (next !== char) {
            throw new Error(`expected ${description}, found ${JSON.stringify(next)}`);
        }
        this.pos++;
    }

    /**
     * Skip whitespace in the buffer, and return the next character, or
     * undefined if the buffer is exhausted
     */
    private skipWhitespace(): string | undefined {
        while (this.pos < this.buf.length && isWhitespace(this.buf.charCodeAt(this.pos))) {
            this.pos++;
        }
        return this.pos < this.buf.length ? this.buf[this.pos] : undefined;
    }

    /**
     * Read the value at the current position if it is complete in the buffer,
     * and return undefined without consuming anything otherwise
     */
    private readBufferedValue(): string | undefined {
        const first = this.buf[this.pos];
        if (first !== '{' && first !== '[' && first !== '"') {
            // let `readValue` report the error
            return undefined;
        }
        const end = makeValueScanner()(this.buf, this.pos);
        if (end < 0) {
            return undefined;
        }
        const value = this.buf.slice(this.pos, end);
        this.pos = end;
        return value;
    }

    /** Read the next chunk of the stream, returning false at its end */
    private async pull(): Promise<boolean> {
        if (this.eof) {
            return false;
        }
        const { value, done } = await this.reader.read();
        if (done) {
            // flush the end of a multi-byte character
            this.buf += this.decoder.decode();
            this.eof = true;
            return false;
        }
        this.buf += this.decoder.decode(value, { stream: true });
        return true;
    }
}
