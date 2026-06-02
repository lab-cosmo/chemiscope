import type { Structure } from '../../src/dataset';
import { parseJsonWithNaN } from './parser';

export const BATCH_SIZE = 2048;

function uid(): string {
    if (typeof crypto !== 'undefined' && typeof crypto.randomUUID === 'function') {
        return crypto.randomUUID();
    }
    return Math.random().toString(36).slice(2) + Date.now().toString(36);
}

interface StoredBatch {
    offsets: number[];
    data: string;
}

/**
 * Stores structure JSON in browser storage so we don't keep millions of them
 * in memory. Writes are batched, recent reads are cached.
 */
export class StructureStore {
    private static readonly NAME_PREFIX = 'chemiscope-stream-';
    private static readonly ORPHAN_MAX_AGE_MS = 60 * 60 * 1000;
    private static readonly READ_CACHE_BATCHES = 4;

    private db: IDBDatabase;
    private readonly storeName = 'structures';
    private readCache = new Map<number, StoredBatch>();

    private constructor(db: IDBDatabase) {
        this.db = db;
    }

    public static async open(): Promise<StructureStore> {
        if (typeof indexedDB === 'undefined') {
            throw new Error('browser storage is not available in this environment');
        }

        await StructureStore.sweepOrphans();

        const dbName = `${StructureStore.NAME_PREFIX}${Date.now()}-${uid()}`;
        const db = await new Promise<IDBDatabase>((resolve, reject) => {
            const req = indexedDB.open(dbName, 1);
            req.onupgradeneeded = () => req.result.createObjectStore('structures');
            req.onsuccess = () => resolve(req.result);
            req.onerror = () => reject(req.error ?? new Error('failed to open browser storage'));
        });

        return new StructureStore(db);
    }

    // drop our own databases left behind by sessions that never called release()
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
            // cleanup is best-effort, ignore failures
        }
    }

    public putBatch(batchIndex: number, raws: string[]): Promise<void> {
        const offsets = new Array<number>(raws.length);
        let total = 0;
        for (let i = 0; i < raws.length; i++) {
            offsets[i] = total;
            total += raws[i].length;
        }
        const batch: StoredBatch = { offsets, data: raws.join('') };

        return new Promise((resolve, reject) => {
            const tx = this.db.transaction(this.storeName, 'readwrite');
            tx.objectStore(this.storeName).put(batch, batchIndex);
            tx.oncomplete = () => resolve();
            tx.onerror = () => reject(tx.error ?? new Error('store write failed'));
        });
    }

    public async get(index: number): Promise<Structure> {
        const batchIndex = Math.floor(index / BATCH_SIZE);
        let batch = this.readCache.get(batchIndex);

        if (batch === undefined) {
            const fetched = await new Promise<StoredBatch | undefined>((resolve, reject) => {
                const tx = this.db.transaction(this.storeName, 'readonly');
                const req = tx.objectStore(this.storeName).get(batchIndex);
                req.onsuccess = () => resolve(req.result as StoredBatch | undefined);
                req.onerror = () => reject(req.error ?? new Error('store read failed'));
            });
            if (fetched === undefined) {
                throw new Error(
                    `structure ${index} not found in store (out of range or released)`
                );
            }
            batch = fetched;
            if (this.readCache.size >= StructureStore.READ_CACHE_BATCHES) {
                // Map iterates in insertion order so the first key is the oldest
                const oldest = this.readCache.keys().next().value;
                if (oldest !== undefined) {
                    this.readCache.delete(oldest);
                }
            }
        } else {
            this.readCache.delete(batchIndex);
        }
        this.readCache.set(batchIndex, batch);

        const indexInBatch = index % BATCH_SIZE;
        if (indexInBatch >= batch.offsets.length) {
            throw new Error(`structure ${index} not found in store (out of range or released)`);
        }
        const start = batch.offsets[indexInBatch];
        const end =
            indexInBatch + 1 < batch.offsets.length
                ? batch.offsets[indexInBatch + 1]
                : batch.data.length;
        return parseJsonWithNaN(batch.data.slice(start, end)) as Structure;
    }

    public async release(): Promise<void> {
        const name = this.db.name;
        this.readCache.clear();
        this.db.close();
        await new Promise<void>((resolve) => {
            const req = indexedDB.deleteDatabase(name);
            req.onsuccess = req.onerror = req.onblocked = () => resolve();
        });
    }
}
