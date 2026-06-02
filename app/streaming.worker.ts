/// <reference lib="WebWorker" />

import { loadDatasetStreaming, StreamingDataset, StreamingProgress } from './streaming';
import type { Structure } from '../src/dataset';

declare const self: DedicatedWorkerGlobalScope;

type InMessage =
    | { type: 'load'; file: File }
    | { type: 'loadStructure'; index: number; requestId: number }
    | { type: 'release' };

type OutMessage =
    | { type: 'progress'; progress: StreamingProgress }
    | { type: 'done'; dataset: unknown }
    | { type: 'structure'; requestId: number; structure: Structure }
    | { type: 'structureError'; requestId: number; message: string }
    | { type: 'released' }
    | { type: 'error'; message: string };

function post(msg: OutMessage): void {
    self.postMessage(msg);
}

let streamingDataset: StreamingDataset | undefined;

self.onmessage = async (event: MessageEvent<InMessage>) => {
    const msg = event.data;
    try {
        if (msg.type === 'load') {
            streamingDataset = await loadDatasetStreaming(msg.file, (progress) => {
                post({ type: 'progress', progress });
            });
            post({ type: 'done', dataset: streamingDataset.dataset });
        } else if (msg.type === 'loadStructure') {
            if (streamingDataset === undefined) {
                post({
                    type: 'structureError',
                    requestId: msg.requestId,
                    message: 'loadStructure called before load completed',
                });
                return;
            }
            try {
                const structure = await streamingDataset.loadStructure(msg.index);
                post({ type: 'structure', requestId: msg.requestId, structure });
            } catch (error) {
                post({
                    type: 'structureError',
                    requestId: msg.requestId,
                    message: error instanceof Error ? error.message : String(error),
                });
            }
        } else if (msg.type === 'release') {
            if (streamingDataset !== undefined) {
                await streamingDataset.release();
                streamingDataset = undefined;
            }
            post({ type: 'released' });
        }
    } catch (error) {
        post({
            type: 'error',
            message: error instanceof Error ? error.message : String(error),
        });
    }
};
