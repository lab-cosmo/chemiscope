// Minimal ambient declarations to satisfy TypeScript in environments
// where `streamlit-component-lib` does not provide its own type defs.
declare module "streamlit-component-lib" {
  export const Streamlit: any;
  export type RenderData = any;
  export default Streamlit;
}
