import type { Dataset, Environment, Property, Structure, UserStructure } from '../src/dataset';

// how many structures go into a single store record
const BATCH_SIZE = 2048;
// how many store writes can run in the background before we wait
const MAX_PENDING_WRITES = 32;
// drop the consumed part of the parse buffer past this size
const COMPACT_THRESHOLD = 1 << 16;
// progress updates are throttled to keep the displayed numbers readable
const PROGRESS_THROTTLE_MS = 100;

const WHITESPACE = new Set([' ', '\t', '\n', '\r']);

export interface StreamingDataset {
    dataset: Dataset;
    loadStructure: (index: number) => Promise<Structure>;
    release: () => Promise<void>;
}

export interface StreamingProgress {
    bytesRead: number;
    bytesTotal: number;
    structures: number;
    phase:
        | 'structures'
        | 'flushing'
        | 'properties'
        | 'environments'
        | 'finalizing'
        | 'preparing';
}

export function parseJsonWithNaN(text: string): unknown {
    if (!text.includes('NaN')) {
        return JSON.parse(text) as unknown;
    }

    let out = '';
    let chunkStart = 0;
    let inString = false;
    let escape = false;
    for (let i = 0; i < text.length; i++) {
        const c = text[i];
        if (inString) {
            if (escape) {
                escape = false;
            } else if (c === '\\') {
                escape = true;
            } else if (c === '"') {
                inString = false;
            }
            continue;
        }
        if (c === '"') {
            inString = true;
            continue;
        }
        if (c === 'N' && text[i + 1] === 'a' && text[i + 2] === 'N') {
            out += text.slice(chunkStart, i) + '"***NaN***"';
            i += 2;
            chunkStart = i + 1;
        }
    }
    out += text.slice(chunkStart);
    return JSON.parse(out, (_key, value: unknown) =>
        value === '***NaN***' ? NaN : value
    ) as unknown;
}

function makeDepthScanner(): (buf: string, from: number) => number {
    let depth = 0;
    let inString = false;
    let escape = false;
    let started = false;

    return function scan(buf: string, from: number): number {
        for (let i = from; i < buf.length; i++) {
            const char = buf[i];

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
                    return i + 1;
                }
            }
        }
        return -1;
    };
}

function extractSize(raw: string, index: number): number {
    let inString = false;
    let escape = false;

    for (let i = 0; i < raw.length; i++) {
        const c = raw[i];
        if (inString) {
            if (escape) {
                escape = false;
            } else if (c === '\\') {
                escape = true;
            } else if (c === '"') {
                inString = false;
            }
            continue;
        }
        if (c !== '"') {
            continue;
        }
        if (!raw.startsWith('"size"', i)) {
            inString = true;
            continue;
        }

        let j = i + 6;
        while (j < raw.length && WHITESPACE.has(raw[j])) {
            j++;
        }
        if (raw[j] !== ':') {
            // a string value that happens to be "size" — skip it
            i += 5;
            continue;
        }
        j++;
        while (j < raw.length && WHITESPACE.has(raw[j])) {
            j++;
        }
        const start = j;
        while (j < raw.length && raw[j] >= '0' && raw[j] <= '9') {
            j++;
        }
        if (j > start) {
            return parseInt(raw.slice(start, j), 10);
        }
        break;
    }
    throw new Error(`structure ${index} is missing an integer "size"`);
}

class StructureStore {
    private static readonly NAME_PREFIX = 'chemiscope-stream-';
    private static readonly ORPHAN_MAX_AGE_MS = 60 * 60 * 1000;

    private db: IDBDatabase;
    private readonly storeName = 'structures';
    private cacheIndex = -1;
    private cache: string[] = [];

    private constructor(db: IDBDatabase) {
        this.db = db;
    }

    public static async open(): Promise<StructureStore> {
        if (typeof indexedDB === 'undefined') {
            throw new Error('browser storage is not available in this environment');
        }

        await StructureStore.sweepOrphans();

        const dbName = `${StructureStore.NAME_PREFIX}${Date.now()}-${Math.random().toString(36).slice(2)}`;
        const db = await new Promise<IDBDatabase>((resolve, reject) => {
            const req = indexedDB.open(dbName, 1);
            req.onupgradeneeded = () => req.result.createObjectStore('structures');
            req.onsuccess = () => resolve(req.result);
            req.onerror = () => reject(req.error ?? new Error('failed to open browser storage'));
        });

        return new StructureStore(db);
    }

    private static async sweepOrphans(): Promise<void> {
        const factory = indexedDB as IDBFactory & {
            databases?: () => Promise<{ name?: string }[]>;
        };
        if (typeof factory.databases !== 'function') {
            // not available on Firefox
            return;
        }

        try {
            const dbs = await factory.databases();
            const now = Date.now();
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
            // cleanup
        }
    }

    public putBatch(batchIndex: number, raws: string[]): Promise<void> {
        return new Promise((resolve, reject) => {
            const tx = this.db.transaction(this.storeName, 'readwrite');
            tx.objectStore(this.storeName).put(raws, batchIndex);
            tx.oncomplete = () => resolve();
            tx.onerror = () => reject(tx.error ?? new Error('store write failed'));
        });
    }

    public async get(index: number): Promise<Structure> {
        const batchIndex = Math.floor(index / BATCH_SIZE);

        if (batchIndex !== this.cacheIndex) {
            this.cache = await new Promise<string[]>((resolve, reject) => {
                const tx = this.db.transaction(this.storeName, 'readonly');
                const req = tx.objectStore(this.storeName).get(batchIndex);
                req.onsuccess = () => resolve((req.result as string[]) ?? []);
                req.onerror = () => reject(req.error ?? new Error('store read failed'));
            });
            this.cacheIndex = batchIndex;
        }

        const raw = this.cache[index % BATCH_SIZE];
        if (raw === undefined) {
            throw new Error(`structure ${index} not found in store (out of range or released)`);
        }
        return parseJsonWithNaN(raw) as Structure;
    }

    public async release(): Promise<void> {
        const name = this.db.name;
        this.cacheIndex = -1;
        this.cache = [];
        this.db.close();
        await new Promise<void>((resolve) => {
            const req = indexedDB.deleteDatabase(name);
            req.onsuccess = req.onerror = req.onblocked = () => resolve();
        });
    }
}

/** Streaming JSON tokenizer over a byte ReadableStream */
class StreamingJSONParser {
    private buf = '';
    private pos = 0;
    private eof = false;
    private decoder = new TextDecoder('utf-8');

    constructor(private reader: ReadableStreamDefaultReader<Uint8Array>) {}

    /** drop the consumed prefix so the buffer doesn't grow without bound */
    private compact(): void {
        if (this.pos > 0) {
            this.buf = this.buf.slice(this.pos);
            this.pos = 0;
        }
    }

    public consume(): void {
        this.pos++;
    }

    /** return the next non-whitespace character without advancing past it */
    public async peek(): Promise<string | undefined> {
        if (this.pos > COMPACT_THRESHOLD) {
            this.compact();
        }
        while (true) {
            while (this.pos < this.buf.length && WHITESPACE.has(this.buf[this.pos])) {
                this.pos++;
            }
            if (this.pos < this.buf.length) {
                return this.buf[this.pos];
            }
            if (!(await this.pull())) {
                return undefined;
            }
        }
    }

    /** read a JSON string starting at the current " and consume past it */
    public async readKey(): Promise<string> {
        let end = -1;
        let escape = false;
        let from = this.pos + 1;
        while (true) {
            for (let i = from; i < this.buf.length; i++) {
                const char = this.buf[i];
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
            from = this.buf.length;
            if (!(await this.pull())) {
                throw new Error('unterminated JSON string');
            }
        }

        const raw = this.buf.slice(this.pos, end + 1);
        this.pos = end + 1;
        return JSON.parse(raw) as string;
    }

    public async captureValue(): Promise<string> {
        const scan = makeDepthScanner();
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
            this.pos = this.buf.length;
            this.compact();
            start = 0;
            if (!(await this.pull())) {
                throw new Error('unterminated JSON value');
            }
        }
    }

    private async pull(): Promise<boolean> {
        if (this.eof) {
            return false;
        }
        const { value, done } = await this.reader.read();
        if (done) {
            this.buf += this.decoder.decode(); // flush trailing partial multibyte char
            this.eof = true;
            return false;
        }
        this.buf += this.decoder.decode(value, { stream: true });
        return true;
    }
}

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

    // count compressed bytes so progress is against file.size
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

    const parser = new StreamingJSONParser(byteStream.getReader());
    const store = await StructureStore.open();

    // progress and phase

    const placeholders: UserStructure[] = [];
    let phase: StreamingProgress['phase'] = 'structures';
    let lastProgressMs = 0;

    function reportProgress(force: boolean = false): void {
        if (!onProgress) {
            return;
        }
        const now = typeof performance !== 'undefined' ? performance.now() : Date.now();
        if (!force && now - lastProgressMs < PROGRESS_THROTTLE_MS) {
            return;
        }
        lastProgressMs = now;
        onProgress({
            bytesRead,
            bytesTotal: file.size,
            structures: placeholders.length,
            phase,
        });
    }

    async function setPhase(next: StreamingProgress['phase']): Promise<void> {
        phase = next;
        reportProgress(true);

        await new Promise<void>((resolve) => setTimeout(resolve, 32));
    }

    // batching and writing structures to the store

    let batch: string[] = [];
    let batchIndex = 0;
    let pendingWrites: Promise<void>[] = [];
    let writeError: Error | undefined;

    function throwIfWriteFailed(): void {
        if (writeError !== undefined) {
            throw writeError;
        }
    }

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

        // pause when too many writes are already pending
        if (pendingWrites.length >= MAX_PENDING_WRITES) {
            const group = pendingWrites;
            pendingWrites = [];
            await Promise.all(group);
            throwIfWriteFailed();
        }
    }

    async function parseProperties(): Promise<Record<string, Property>> {
        if ((await parser.peek()) !== '{') {
            throw new Error('"properties" must be an object');
        }
        parser.consume();

        const properties: Record<string, Property> = {};
        while (true) {
            const char = await parser.peek();
            if (char === '}') {
                parser.consume();
                break;
            }
            if (char === ',') {
                parser.consume();
                continue;
            }
            if (char !== '"') {
                throw new Error('expected property name');
            }

            const name = await parser.readKey();
            if ((await parser.peek()) !== ':') {
                throw new Error('expected ":" after property');
            }
            parser.consume();
            if ((await parser.peek()) !== '{') {
                throw new Error(`property "${name}" must be an object`);
            }

            properties[name] = parseJsonWithNaN(await parser.captureValue()) as Property;
            reportProgress();
        }
        return properties;
    }

    async function parseEnvironments(): Promise<Environment[]> {
        if ((await parser.peek()) !== '[') {
            throw new Error('"environments" must be an array');
        }
        const raw = await parser.captureValue();
        return parseJsonWithNaN(raw) as Environment[];
    }

    async function streamStructures(): Promise<void> {
        if ((await parser.peek()) !== '[') {
            throw new Error('"structures" must be an array');
        }
        parser.consume();

        while (true) {
            const char = await parser.peek();
            if (char === undefined) {
                throw new Error('unexpected EOF in "structures"');
            }
            if (char === ']') {
                parser.consume();
                break;
            }
            if (char === ',') {
                parser.consume();
                continue;
            }

            const raw = await parser.captureValue();
            const size = extractSize(raw, placeholders.length);
            placeholders.push({ size, data: placeholders.length });
            batch.push(raw);

            if (batch.length >= BATCH_SIZE) {
                await flushBatch();
            }
            if (placeholders.length % 5000 === 0) {
                reportProgress();
                await new Promise<void>((resolve) => setTimeout(resolve, 0));
            }
        }

        // wait for every outstanding write before properties/environments
        await setPhase('flushing');
        await flushBatch();
        await Promise.all(pendingWrites);
        pendingWrites = [];
        await setPhase('finalizing');
        throwIfWriteFailed();
    }

    // top-level dispatch

    try {
        const sections: Record<string, string> = {};
        let properties: Record<string, Property> = {};
        let environments: Environment[] | undefined;

        if ((await parser.peek()) !== '{') {
            throw new Error('dataset must be a JSON object');
        }
        parser.consume();

        while (true) {
            const char = await parser.peek();
            if (char === '}') {
                parser.consume();
                break;
            }
            if (char === ',') {
                parser.consume();
                continue;
            }
            if (char !== '"') {
                throw new Error(`expected object key, found ${JSON.stringify(char)}`);
            }

            const key = await parser.readKey();
            if ((await parser.peek()) !== ':') {
                throw new Error('expected ":" after key');
            }
            parser.consume();

            // huge sections are handled with dedicated streaming parsers;
            // small ones are captured as raw text and parsed at the end
            if (key === 'structures') {
                await setPhase('structures');
                await streamStructures();
            } else if (key === 'properties') {
                await setPhase('properties');
                properties = await parseProperties();
            } else if (key === 'environments') {
                await setPhase('environments');
                environments = await parseEnvironments();
            } else {
                await parser.peek();
                sections[key] = await parser.captureValue();
            }
        }

        await setPhase('finalizing');

        const dataset = {} as Dataset;
        for (const [key, raw] of Object.entries(sections)) {
            (dataset as unknown as Record<string, unknown>)[key] = parseJsonWithNaN(raw);
        }
        dataset.properties = properties;
        dataset.structures = placeholders;
        if (environments !== undefined) {
            dataset.environments = environments;
        }

        await setPhase('preparing');

        return {
            dataset,
            loadStructure: (index: number) => store.get(index),
            release: () => store.release(),
        };
    } catch (error) {
        await store.release().catch(() => {});
        throw writeError ?? error;
    }
}
