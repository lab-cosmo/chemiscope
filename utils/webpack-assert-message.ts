/**
 * A webpack loader to use the assertion expression as error message if the
 * assertion fails.
 *
 * ```
 * assert(foo && bar.baz == "test");
 * ```
 * is transformed to
 * ```
 * assert(foo && bar.baz == "test", "foo && bar.baz == \"test\"");
 * ```
 *
 * The replacement mechanism is not very robust since it is based on regex
 * instead of a full AST, but it should be good enough to work with simple,
 * single line assertions.
 */
export default function (source: string): string {
    return source.replace(/assert\((.*)\)/g, (_, expr) => {
        const escaped = expr.replace(/"/g, '\\"');
        return `assert(${expr}, "${escaped}")`;
    });
}
