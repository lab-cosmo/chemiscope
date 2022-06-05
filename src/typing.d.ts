declare module '*.html' {
    const content: string;
    export default content;
}

declare module '*.svg' {
    const content: string;
    export default content;
}

declare module '*.css';

declare module '!css-loader?exportType=css-style-sheet!*.css' {
    const content: CSSStyleSheet;
    export default content;
}
