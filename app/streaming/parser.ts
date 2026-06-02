const WHITESPACE = new Set([' ', '\t', '\n', '\r']);

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
        if (
            c === 'N' &&
            text[i + 1] === 'a' &&
            text[i + 2] === 'N' &&
            !isWordChar(text[i - 1]) &&
            !isWordChar(text[i + 3])
        ) {
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

function isWordChar(c: string | undefined): boolean {
    return (
        c !== undefined &&
        ((c >= 'a' && c <= 'z') || (c >= 'A' && c <= 'Z') || (c >= '0' && c <= '9') || c === '_')
    );
}

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

/** Streaming JSON tokenizer over a byte ReadableStream */
export class StreamingJSONParser {
    private buf = '';
    private pos = 0;
    private eof = false;
    private decoder = new TextDecoder('utf-8');

    constructor(private reader: ReadableStreamDefaultReader<Uint8Array>) {}

    public get position(): number {
        return this.pos;
    }

    /** drop the consumed prefix so the buffer doesn't grow without bound */
    public compact(): void {
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

    /** read one balanced `{...}` or `[...]` as raw text */
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

    /** read a bare token (number / true / false / null) up to the next delimiter */
    public async readBareToken(): Promise<string> {
        let end = -1;
        while (true) {
            for (let i = this.pos; i < this.buf.length; i++) {
                const char = this.buf[i];
                if (char === ',' || char === ']' || char === '}' || WHITESPACE.has(char)) {
                    end = i;
                    break;
                }
            }
            if (end >= 0) {
                break;
            }
            if (!(await this.pull())) {
                end = this.buf.length;
                break;
            }
        }

        const token = this.buf.slice(this.pos, end);
        this.pos = end;
        return token;
    }

    /** dispatch on the next non-whitespace char and return whatever value comes next */
    public async readAnyValue(): Promise<unknown> {
        const char = await this.peek();
        if (char === '"') {
            return this.readKey();
        }
        if (char === '{' || char === '[') {
            return parseJsonWithNaN(await this.captureValue());
        }
        const token = await this.readBareToken();
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
