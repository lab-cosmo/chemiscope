import type { Dataset, Environment, Property, Structure, UserStructure } from '../../src/dataset';
import { BATCH_SIZE, StructureStore } from './store';
import { StreamingJSONParser, parseJsonWithNaN } from './parser';

export { parseJsonWithNaN };

const MAX_PENDING_WRITES = 32;
const COMPACT_THRESHOLD = 1 << 16;
const PROGRESS_THROTTLE_MS = 100;

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

const WHITESPACE = new Set([' ', '\t', '\n', '\r']);

function extractSize(raw: string, index: number): number {
    const missing = () => new Error(`structure ${index} is missing an integer "size"`);
    let depth = 0;
    let inString = false;
    let escape = false;
    let expectKey = false;

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
        if (c === '"') {
            if (depth === 1 && expectKey && raw.startsWith('"size"', i)) {
                let j = i + 6;
                while (j < raw.length && WHITESPACE.has(raw[j])) {
                    j++;
                }
                if (raw[j] !== ':') {
                    throw missing();
                }
                j++;
                while (j < raw.length && WHITESPACE.has(raw[j])) {
                    j++;
                }
                const start = j;
                while (j < raw.length && raw[j] >= '0' && raw[j] <= '9') {
                    j++;
                }
                if (j === start || raw[j] === '.' || raw[j] === 'e' || raw[j] === 'E') {
                    throw missing();
                }
                return parseInt(raw.slice(start, j), 10);
            }
            inString = true;
            expectKey = false;
        } else if (c === '{') {
            depth++;
            expectKey = depth === 1;
        } else if (c === '[') {
            depth++;
            expectKey = false;
        } else if (c === '}' || c === ']') {
            depth--;
            expectKey = false;
        } else if (c === ',' && depth === 1) {
            expectKey = true;
        } else if (c === ':' && depth === 1) {
            expectKey = false;
        }
    }
    throw missing();
}

/**
 * Stream a chemiscope dataset from a File: the huge structures array is
 * offloaded to browser storage and replaced in the returned dataset with
 * lightweight placeholders, while everything else is parsed in memory.
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

        // let the UI repaint the new message before we resume parsing
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

    async function drainWrites(): Promise<void> {
        await flushBatch();
        const group = pendingWrites;
        pendingWrites = [];
        await Promise.all(group);
        await setPhase('finalizing');
        throwIfWriteFailed();
    }

    // chemiscope-specific section parsers

    async function parseValuesArray(): Promise<(number | string | number[])[]> {
        if ((await parser.peek()) !== '[') {
            throw new Error('"values" must be an array');
        }
        const raw = await parser.captureValue();
        return parseJsonWithNaN(raw) as (number | string | number[])[];
    }

    async function parseProperty(): Promise<Property> {
        if ((await parser.peek()) !== '{') {
            throw new Error('property must be an object');
        }
        parser.consume();

        const prop: Record<string, unknown> = {};
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
                throw new Error('expected property field name');
            }

            const field = await parser.readKey();
            if ((await parser.peek()) !== ':') {
                throw new Error('expected ":" after field');
            }
            parser.consume();

            if (field === 'values') {
                prop.values = await parseValuesArray();
            } else {
                prop[field] = await parser.readAnyValue();
            }
        }
        return prop as unknown as Property;
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

            properties[name] = await parseProperty();
            parser.compact();
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
            if (parser.position > COMPACT_THRESHOLD) {
                parser.compact();
            }
            if (placeholders.length % 5000 === 0) {
                reportProgress();
            }
        }

        await setPhase('flushing');
        await drainWrites();
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
            await parser.peek();

            // huge sections are handled with dedicated streaming parsers;
            // small ones are captured as raw text and parsed at the end
            if (key === 'structures') {
                await setPhase('structures');
                await streamStructures();
            } else if (key === 'properties') {
                await setPhase('properties');
                properties = await parseProperties();
                parser.compact();
            } else if (key === 'environments') {
                await setPhase('environments');
                environments = await parseEnvironments();
                parser.compact();
            } else {
                sections[key] = await parser.captureValue();
                parser.compact();
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
        // release the store on failure, otherwise it's orphaned with partial data
        await store.release().catch(() => {});
        throw writeError ?? error;
    }
}
