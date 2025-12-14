declare module "streamlit-component-lib" {
  export interface RenderData {
    args: Record<string, any>;
    disabled: boolean;
    theme?: {
      base: "light" | "dark";
      primaryColor: string;
      backgroundColor: string;
      secondaryBackgroundColor: string;
      textColor: string;
      font: string;
    };
  }

  export interface StreamlitAPI {
    readonly RENDER_EVENT: string;
    events: {
      readonly RENDER_EVENT: string;
      addEventListener(event: string, callback: (event: CustomEvent<RenderData>) => void): void;
    };

    setComponentReady(): void;
    setComponentValue(value: any): void;
    setFrameHeight(height?: number): void;
  }

  export const Streamlit: StreamlitAPI;
  export default Streamlit;
}
