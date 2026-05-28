import type { Dataset, Environment, Property, Structure, UserStructure } from '../src/dataset';

// structures per IndexedDB record
const BATCH_SIZE = 256;

// batch writes kept in flight
const MAX_INFLIGHT_WRITES = 32;

// only drop the consumed buffer prefix past this size
const COMPACT_THRESHOLD = 1 << 16;

// unique id for the store name
function uid(): string {
    if (typeof crypto !== 'undefined' && typeof crypto.randomUUID === 'function') {
        return crypto.randomUUID();
    }
    return Math.random().toString(36).slice(2) + Date.now().toString(36);
}

export interface StreamingDataset {
    // dataset with lightweight structure placeholders (`{ size, data: index }`)
    dataset: Dataset;
    // on-demand structure loader
    loadStructure: (index: number) => Promise<Structure>;
    // release the backing IndexedDB store
    release: () => Promise<void>;
}

export interface StreamingProgress {
    // compressed bytes read so far
    bytesRead: number;
    // total compressed size, if known
    bytesTotal: number;
    // structures stored so far
    structures: number;
}

// parse JSON, tolerating Python's bare `NaN`
export function parseJsonWithNaN(text: string): unknown {
    return JSON.parse(text.replace(/\bNaN\b/g, '"***NaN***"'), (_key, value: unknown) =>
        value === '***NaN***' ? NaN : value
    ) as unknown;
}

// like Number(token) but throws on non-numbers instead of silently returning
// 0 (for "") or NaN
function parseNumberToken(token: string): number {
    if (token === 'NaN') {
        return NaN;
    }
    const value = Number(token);
    if (token === '' || Number.isNaN(value)) {
        throw new Error(`invalid number in dataset: ${JSON.stringify(token)}`);
    }
    return value;
}

// raw structure JSON in IndexedDB, batched, with a one-batch read cache
class StructureStore {
    // db names are `${NAME_PREFIX}${timestamp}-${uid}`
    private static readonly NAME_PREFIX = 'chemiscope-stream-';
    private static readonly ORPHAN_MAX_AGE_MS = 60 * 60 * 1000; // 1 hour

    private db: IDBDatabase;
    private readonly storeName = 'structures';

    // which batch is currently held in `cache` (-1 = none)
    private cacheBatch = -1;
    private cache: string[] = [];

    private constructor(db: IDBDatabase) {
        this.db = db;
    }

    public static async open(): Promise<StructureStore> {
        if (typeof indexedDB === 'undefined') {
            throw new Error('IndexedDB is not available in this environment');
        }

        // reap leftovers from crashed sessions before adding ours
        await StructureStore.sweepOrphans();

        // unique name so back-to-back loads never collide
        const dbName = `${StructureStore.NAME_PREFIX}${Date.now()}-${uid()}`;

        const db = await new Promise<IDBDatabase>((resolve, reject) => {
            const req = indexedDB.open(dbName, 1);
            req.onupgradeneeded = () => req.result.createObjectStore('structures');
            req.onsuccess = () => resolve(req.result);
            req.onerror = () => reject(req.error ?? new Error('failed to open IndexedDB'));
        });

        return new StructureStore(db);
    }

    // drop stores left behind by sessions that never called release()
    private static async sweepOrphans(): Promise<void> {
        const factory = indexedDB as IDBFactory & {
            databases?: () => Promise<{ name?: string }[]>;
        };
        if (typeof factory.databases !== 'function') {
            return; // e.g. Firefox
        }

        try {
            const dbs = await factory.databases();
            const now = Date.now();

            // keep only our own databases that are old enough to be orphans
            const stale = dbs
                .map((d) => d.name ?? '')
                .filter((name) => {
                    if (!name.startsWith(StructureStore.NAME_PREFIX)) {
                        return false;
                    }
                    const ts = Number(name.slice(StructureStore.NAME_PREFIX.length).split('-')[0]);
                    return Number.isFinite(ts) && now - ts > StructureStore.ORPHAN_MAX_AGE_MS;
                });

            await Promise.all(
                stale.map(
                    (name) =>
                        new Promise<void>((resolve) => {
                            const req = indexedDB.deleteDatabase(name);
                            req.onsuccess = req.onerror = req.onblocked = () => resolve();
                        })
                )
            );
        } catch {
            // cleanup is best-effort, ignore failures
        }
    }

    // persist one batch under its index
    public putBatch(batchIndex: number, raws: string[]): Promise<void> {
        return new Promise((resolve, reject) => {
            const tx = this.db.transaction(this.storeName, 'readwrite');
            tx.objectStore(this.storeName).put(raws, batchIndex);
            tx.oncomplete = () => resolve();
            tx.onerror = () => reject(tx.error ?? new Error('IndexedDB write failed'));
        });
    }

    public async get(index: number): Promise<Structure> {
        const batchIndex = Math.floor(index / BATCH_SIZE);

        // only hit IndexedDB when we move to a different batch
        if (batchIndex !== this.cacheBatch) {
            this.cache = await new Promise<string[]>((resolve, reject) => {
                const tx = this.db.transaction(this.storeName, 'readonly');
                const req = tx.objectStore(this.storeName).get(batchIndex);
                req.onsuccess = () => resolve((req.result as string[]) ?? []);
                req.onerror = () => reject(req.error ?? new Error('IndexedDB read failed'));
            });
            this.cacheBatch = batchIndex;
        }

        const raw = this.cache[index % BATCH_SIZE];
        if (raw === undefined) {
            throw new Error(`structure ${index} not found in store (out of range or released)`);
        }
        return parseJsonWithNaN(raw) as Structure;
    }

    public async release(): Promise<void> {
        const name = this.db.name;
        this.db.close();

        // the database is single-use, so drop it entirely
        await new Promise<void>((resolve) => {
            const req = indexedDB.deleteDatabase(name);
            req.onsuccess = req.onerror = req.onblocked = () => resolve();
        });
    }
}

// finds where one `{...}`/`[...]` value ends, or -1 if not all read yet
function makeDepthScanner(): (buf: string, from: number) => number {
    let depth = 0;
    let inString = false;
    let escape = false;
    let started = false;

    return function scan(buf: string, from: number): number {
        for (let i = from; i < buf.length; i++) {
            const char = buf[i];

            // inside a string, only an unescaped quote ends it
            if (inString) {
                if (escape) {
                    escape = false;
                } else if (char === '\\') {
                    escape = true;
                } else if (char === '"') {
                    inString = false;
                }
                continue;
            }

            if (char === '"') {
                inString = true;
            } else if (char === '{' || char === '[') {
                depth++;
                started = true;
            } else if (char === '}' || char === ']') {
                depth--;
                if (started && depth === 0) {
                    return i + 1; // closed the top-level value
                }
            }
        }
        return -1;
    };
}

const WHITESPACE = new Set([' ', '\t', '\n', '\r']);

/**
 * stream a gzipped or plain chemiscope dataset from a `File`, routing structures
 * to IndexedDB and parsing the remaining (small) sections in memory
 */
export async function loadDatasetStreaming(
    file: File,
    onProgress?: (progress: StreamingProgress) => void
): Promise<StreamingDataset> {
    const isGzipped = file.name.endsWith('.gz');
    if (isGzipped && typeof DecompressionStream === 'undefined') {
        throw new Error(
            'streaming load of gzipped files requires DecompressionStream (unavailable in this browser)'
        );
    }

    // count compressed bytes (pre-decompression) so progress is against file.size
    let bytesRead = 0;
    const counting = new TransformStream<Uint8Array, Uint8Array>({
        transform(chunk, controller) {
            bytesRead += chunk.byteLength;
            controller.enqueue(chunk);
        },
    });

    let byteStream: ReadableStream<Uint8Array> = file.stream().pipeThrough(counting);
    if (isGzipped) {
        const gunzip = new DecompressionStream('gzip') as unknown as ReadableWritablePair<
            Uint8Array,
            Uint8Array
        >;
        byteStream = byteStream.pipeThrough(gunzip);
    }

    const reader = byteStream.getReader();
    const decoder = new TextDecoder('utf-8');

    const store = await StructureStore.open();

    // sliding buffer over the decoded stream
    let buf = '';
    let pos = 0;
    let eof = false;

    // pull the next chunk onto `buf`
    async function pull(): Promise<boolean> {
        if (eof) {
            return false;
        }
        const { value, done } = await reader.read();
        if (done) {
            buf += decoder.decode(); // flush any trailing partial multibyte char
            eof = true;
            return false;
        }
        buf += decoder.decode(value, { stream: true });
        return true;
    }

    // drop the already-read prefix so the buffer doesn't grow without bound
    function compact(): void {
        if (pos > 0) {
            buf = buf.slice(pos);
            pos = 0;
        }
    }

    async function skipWhitespace(): Promise<string | undefined> {
        while (true) {
            while (pos < buf.length && WHITESPACE.has(buf[pos])) {
                pos++;
            }
            if (pos < buf.length) {
                return buf[pos];
            }
            if (!(await pull())) {
                return undefined;
            }
        }
    }

    // read a JSON string at `pos`
    async function readKey(): Promise<string> {
        let end = -1;
        while (true) {
            let escape = false;

            // look for the closing quote respecting escapes
            for (let i = pos + 1; i < buf.length; i++) {
                const char = buf[i];
                if (escape) {
                    escape = false;
                } else if (char === '\\') {
                    escape = true;
                } else if (char === '"') {
                    end = i;
                    break;
                }
            }

            if (end >= 0) {
                break;
            }
            if (!(await pull())) {
                throw new Error('unterminated JSON string');
            }
        }

        const raw = buf.slice(pos, end + 1);
        pos = end + 1;
        return JSON.parse(raw) as string;
    }

    // read one whole object/array value as raw text. `pos` must sit on the opening bracket
    async function captureValue(): Promise<string> {
        const scan = makeDepthScanner();
        const parts: string[] = [];
        let start = pos;

        while (true) {
            const end = scan(buf, pos);
            if (end >= 0) {
                parts.push(buf.slice(start, end));
                pos = end;
                return parts.join('');
            }

            // value runs past the buffer: stash what we have then pull more
            parts.push(buf.slice(start));
            pos = buf.length;
            compact();
            start = 0;
            if (!(await pull())) {
                throw new Error('unterminated JSON value');
            }
        }
    }

    // read a bare token up to the next delimiter
    async function readBareToken(): Promise<string> {
        let end = -1;
        while (true) {
            for (let i = pos; i < buf.length; i++) {
                const char = buf[i];
                if (char === ',' || char === ']' || char === '}' || WHITESPACE.has(char)) {
                    end = i;
                    break;
                }
            }
            if (end >= 0) {
                break;
            }
            if (!(await pull())) {
                end = buf.length; // token runs to EOF
                break;
            }
        }

        const token = buf.slice(pos, end);
        pos = end;
        return token;
    }

    async function readAnyValue(): Promise<unknown> {
        const char = await skipWhitespace();

        if (char === '"') {
            return readKey();
        }
        if (char === '{' || char === '[') {
            return parseJsonWithNaN(await captureValue());
        }

        // bare literal or number
        const token = await readBareToken();
        if (token === 'true') {
            return true;
        }
        if (token === 'false') {
            return false;
        }
        if (token === 'null') {
            return null;
        }
        return parseNumberToken(token);
    }

    // stream a property's `values` array one element at a time
    async function parseValuesArray(): Promise<(number | string | number[])[]> {
        if ((await skipWhitespace()) !== '[') {
            throw new Error('"values" must be an array');
        }
        pos++;

        const out: (number | string | number[])[] = [];
        while (true) {
            if (pos > COMPACT_THRESHOLD) {
                compact();
            }

            const char = await skipWhitespace();
            if (char === undefined) {
                throw new Error('unterminated "values" array');
            }
            if (char === ']') {
                pos++;
                break;
            }
            if (char === ',') {
                pos++;
                continue;
            }

            if (char === '"') {
                out.push(await readKey());
            } else if (char === '[') {
                // multidimensional property
                out.push(parseJsonWithNaN(await captureValue()) as number[]);
            } else {
                out.push(parseNumberToken(await readBareToken()));
            }

            if (out.length % 50000 === 0) {
                reportProgress();
            }
        }
        return out;
    }

    // parse one property: stream its `values`, read the rest of the fields normally
    async function parseProperty(): Promise<Property> {
        if ((await skipWhitespace()) !== '{') {
            throw new Error('property must be an object');
        }
        pos++;

        const prop: Record<string, unknown> = {};
        while (true) {
            const char = await skipWhitespace();
            if (char === '}') {
                pos++;
                break;
            }
            if (char === ',') {
                pos++;
                continue;
            }
            if (char !== '"') {
                throw new Error('expected property field name');
            }

            const field = await readKey();
            if ((await skipWhitespace()) !== ':') {
                throw new Error('expected ":" after field');
            }
            pos++;

            // `values` is the big one we stream, everything else is small
            if (field === 'values') {
                prop.values = await parseValuesArray();
            } else {
                prop[field] = await readAnyValue();
            }
        }
        return prop as unknown as Property;
    }

    async function parseEnvironments(): Promise<Environment[]> {
        if ((await skipWhitespace()) !== '[') {
            throw new Error('"environments" must be an array');
        }
        pos++;

        const environments: Environment[] = [];
        while (true) {
            if (pos > COMPACT_THRESHOLD) {
                compact();
            }

            const char = await skipWhitespace();
            if (char === undefined) {
                throw new Error('unterminated "environments" array');
            }
            if (char === ']') {
                pos++;
                break;
            }
            if (char === ',') {
                pos++;
                continue;
            }

            // each entry is a small `{structure, center, cutoff}` object
            environments.push(parseJsonWithNaN(await captureValue()) as Environment);

            if (environments.length % 50000 === 0) {
                reportProgress();
            }
        }
        return environments;
    }

    // parse the whole `properties` object, one property at a time
    async function parseProperties(): Promise<Record<string, Property>> {
        if ((await skipWhitespace()) !== '{') {
            throw new Error('"properties" must be an object');
        }
        pos++;

        const properties: Record<string, Property> = {};
        while (true) {
            const char = await skipWhitespace();
            if (char === '}') {
                pos++;
                break;
            }
            if (char === ',') {
                pos++;
                continue;
            }
            if (char !== '"') {
                throw new Error('expected property name');
            }

            const name = await readKey();
            if ((await skipWhitespace()) !== ':') {
                throw new Error('expected ":" after property');
            }
            pos++;

            properties[name] = await parseProperty();
            compact();
        }
        return properties;
    }

    const placeholders: UserStructure[] = [];
    function reportProgress(): void {
        if (onProgress) {
            onProgress({ bytesRead, bytesTotal: file.size, structures: placeholders.length });
        }
    }

    // structure batches on their way to IndexedDB
    let batch: string[] = [];
    let batchIndex = 0;
    let pendingWrites: Promise<void>[] = [];

    // latches the first failed write so it aborts the load instead of becoming an
    // unhandled rejection
    let writeError: Error | undefined;
    function throwIfWriteFailed(): void {
        const error = writeError;
        if (error !== undefined) {
            throw error;
        }
    }

    // queue the current batch; only blocks once MAX_INFLIGHT_WRITES are outstanding
    async function flushBatch(): Promise<void> {
        throwIfWriteFailed();
        if (batch.length === 0) {
            return;
        }

        pendingWrites.push(
            store.putBatch(batchIndex, batch).catch((error: unknown) => {
                if (writeError === undefined) {
                    writeError = error instanceof Error ? error : new Error(String(error));
                }
            })
        );
        batchIndex++;
        batch = [];

        // apply backpressure once too many writes are in flight
        if (pendingWrites.length >= MAX_INFLIGHT_WRITES) {
            const group = pendingWrites;
            pendingWrites = [];
            await Promise.all(group);
            throwIfWriteFailed();
        }
    }

    // flush the last partial batch and wait for every outstanding write
    async function drainWrites(): Promise<void> {
        await flushBatch();
        const group = pendingWrites;
        pendingWrites = [];
        await Promise.all(group);
        throwIfWriteFailed();
    }

    // walk the `structures` array, sending each element to the store
    async function streamStructures(): Promise<void> {
        if ((await skipWhitespace()) !== '[') {
            throw new Error('"structures" must be an array');
        }
        pos++;

        while (true) {
            const char = await skipWhitespace();
            if (char === undefined) {
                throw new Error('unexpected EOF in "structures"');
            }
            if (char === ']') {
                pos++;
                break;
            }
            if (char === ',') {
                pos++;
                continue;
            }

            const raw = await captureValue();

            const sizeMatch = /"size"\s*:\s*(\d+)/.exec(raw);
            if (sizeMatch === null) {
                throw new Error(`structure ${placeholders.length} is missing an integer "size"`);
            }
            const size = parseInt(sizeMatch[1], 10);

            placeholders.push({ size, data: placeholders.length });
            batch.push(raw);

            if (batch.length >= BATCH_SIZE) {
                await flushBatch();
            }
            if (pos > COMPACT_THRESHOLD) {
                compact();
            }
            if (placeholders.length % 5000 === 0) {
                reportProgress();
            }
        }

        await drainWrites();
    }

    // release the on-disk store on any failure or it's orphaned with partial data
    try {
        const sections: Record<string, string> = {};
        let properties: Record<string, Property> = {};
        let environments: Environment[] | undefined;

        // top-level object: dispatch each key to the right parser
        if ((await skipWhitespace()) !== '{') {
            throw new Error('dataset must be a JSON object');
        }
        pos++;

        while (true) {
            const char = await skipWhitespace();
            if (char === '}') {
                pos++;
                break;
            }
            if (char === ',') {
                pos++;
                continue;
            }
            if (char !== '"') {
                throw new Error(`expected object key, found ${JSON.stringify(char)}`);
            }

            const key = await readKey();
            if ((await skipWhitespace()) !== ':') {
                throw new Error('expected ":" after key');
            }
            pos++;
            await skipWhitespace();

            // big sections stream; small ones are captured raw and parsed below
            if (key === 'structures') {
                await streamStructures();
            } else if (key === 'properties') {
                properties = await parseProperties();
                compact();
            } else if (key === 'environments') {
                environments = await parseEnvironments();
                compact();
            } else {
                sections[key] = await captureValue();
                compact();
            }
        }

        // assemble the dataset: small sections + streamed properties + placeholders
        const dataset = {} as Dataset;
        for (const [key, raw] of Object.entries(sections)) {
            (dataset as unknown as Record<string, unknown>)[key] = parseJsonWithNaN(raw);
        }
        dataset.properties = properties;
        dataset.structures = placeholders;
        if (environments !== undefined) {
            dataset.environments = environments;
        }

        reportProgress();

        return {
            dataset,
            loadStructure: (index: number) => store.get(index),
            release: () => store.release(),
        };
    } catch (error) {
        await store.release().catch(() => {});
        throw error;
    }
}
