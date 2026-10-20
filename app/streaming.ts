/**
 * Streaming dataset loader for the web app.
 *
 * Large datasets are read incrementally with `readDataset` from the library
 * (see `src/reader.ts`). The bulk of a dataset is its `structures` array: the
 * raw JSON text of each structure is written (in batches) to an IndexedDB
 * store, and only a lightweight `UserStructure` placeholder ({ size, data:
 * index }) is kept in memory. This keeps the resident set small even for files
 * larger than RAM.
 *
 * The returned `loadStructure(index)` callback lazily reads a structure back
 * from IndexedDB and parses it on demand (chemiscope only materializes a
 * structure when it is actually displayed). `release()` drops the IndexedDB
 * store; only one streamed dataset is kept alive per page (see `activeRelease`),
 * and the store is also released on `pagehide` and when a new dataset loads.
 *
 * Usage from the app:
 *   if (shouldUseStreaming(file)) {
 *       const { dataset, loadStructure } = await loadDatasetStreaming(file, (p) =>
 *           setLoadingProgress(progressMessage(p))
 *       );
 *       await visualizer.load({ loadStructure }, dataset);
 *   }
 */

import type { Dataset, Structure } from '../src/dataset';
import { ReadProgress, parseJsonWithNaN, readDataset } from '../src/reader';

// how many structures go into a single store record
const BATCH_SIZE = 2048;
// how many store writes can run in the background before we wait
const MAX_PENDING_WRITES = 32;

export interface StreamingDataset {
    dataset: Dataset;
    loadStructure: (index: number) => Promise<Structure>;
    release: () => Promise<void>;
}

export interface StreamingProgress extends ReadProgress {
    /** set after reading the file, while finishing to store the structures */
    phase?: 'flushing' | 'preparing';
}

// files at/above this size (MB) stream; smaller ones load fully in memory
const STREAMING_THRESHOLD_MB = 75;

/** Whether a file should be loaded through the streaming path. */
export function shouldUseStreaming(file: File): boolean {
    const isGzipped = file.name.endsWith('.gz');
    const sizeMB = file.size / (1024 * 1024);
    const canStream =
        typeof indexedDB !== 'undefined' &&
        (!isGzipped || typeof DecompressionStream !== 'undefined');
    return sizeMB >= STREAMING_THRESHOLD_MB && canStream;
}

/** Human-readable loading message for a streaming progress update. */
export function progressMessage(progress: StreamingProgress): string {
    if (progress.phase === 'flushing') {
        return 'finalizing structure data…';
    } else if (progress.phase === 'preparing') {
        return 'preparing visualization…';
    }

    const percent =
        progress.bytesTotal > 0
            ? ` (${((100 * progress.bytesRead) / progress.bytesTotal).toFixed(0)}%)`
            : '';

    switch (progress.section) {
        case 'structures':
            return `loaded ${progress.structures.toLocaleString()} structures${percent}`;
        case 'properties':
            return 'parsing properties…';
        case 'environments':
            return 'parsing environments…';
        default:
            return 'reading dataset…';
    }
}

// Only one streamed dataset is kept alive per page; its store is released when a
// new dataset is loaded or when the tab goes away.
let activeRelease: (() => Promise<void>) | undefined;
let pagehideRegistered = false;

/** Release the currently-active streamed dataset's store, if any. */
export async function releaseActiveStream(): Promise<void> {
    if (activeRelease !== undefined) {
        const release = activeRelease;
        activeRelease = undefined;
        await release().catch(() => {});
    }
}

// drop the streaming store when the tab is going away
function registerPagehideOnce(): void {
    if (pagehideRegistered || typeof window === 'undefined') {
        return;
    }
    pagehideRegistered = true;
    window.addEventListener('pagehide', (event) => {
        if (event.persisted || activeRelease === undefined) {
            return;
        }
        void activeRelease().catch(() => {});
    });
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

/** Wait for the next frame, to let the page render */
function nextFrame(): Promise<void> {
    return new Promise<void>((resolve) => {
        if (typeof requestAnimationFrame !== 'undefined') {
            requestAnimationFrame(() => resolve());
        } else {
            setTimeout(resolve, 0);
        }
    });
}

export async function loadDatasetStreaming(
    file: File,
    onProgress?: (progress: StreamingProgress) => void
): Promise<StreamingDataset> {
    registerPagehideOnce();
    await releaseActiveStream();

    const store = await StructureStore.open();

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

    // write the remaining structures to the store, and wait for all writes
    async function finishWrites(): Promise<void> {
        await flushBatch();
        await Promise.all(pendingWrites);
        pendingWrites = [];
        throwIfWriteFailed();
    }

    try {
        let section = '';
        const dataset = await readDataset(file, {
            // keep placeholders in the dataset, and the structures in the store
            structure: (raw, size, index) => {
                batch.push(raw);
                const placeholder = { size, data: index };
                if (batch.length < BATCH_SIZE) {
                    return placeholder;
                }
                return flushBatch().then(() => placeholder);
            },
            onProgress: async (progress) => {
                if (progress.section !== section) {
                    if (section === 'structures') {
                        // finish storing structures before reading the rest
                        onProgress?.({ ...progress, phase: 'flushing' });
                        await nextFrame();
                        await finishWrites();
                    }
                    section = progress.section;
                    onProgress?.(progress);
                    // let the page render the new section before reading it
                    await nextFrame();
                } else {
                    onProgress?.(progress);
                }
            },
        });
        // in case structures were the last section of the file
        await finishWrites();

        const structures = dataset.structures.length;
        const done = { bytesRead: file.size, bytesTotal: file.size, structures, section };
        onProgress?.({ ...done, phase: 'preparing' });
        await nextFrame();

        const release = async (): Promise<void> => {
            // forget this store if it is still the active one
            if (activeRelease === release) {
                activeRelease = undefined;
            }
            await store.release();
        };
        activeRelease = release;

        return {
            dataset,
            loadStructure: (index: number) => store.get(index),
            release,
        };
    } catch (error) {
        await store.release().catch(() => {});
        throw writeError ?? error;
    }
}
