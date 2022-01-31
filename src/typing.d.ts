declare module '*.html' {
    const content: string;
    export default content;
}

declare module '*.svg' {
    const content: string;
    export default content;
}

declare module '*.css';

declare module '*.css?sheet' {
    const content: CSSStyleSheet;
    export default content;
}
