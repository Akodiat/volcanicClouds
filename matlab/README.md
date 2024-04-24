Follow this procedure to generate js and wasm files:
https://se.mathworks.com/matlabcentral/fileexchange/69973-generatejavascriptusingmatlabcoder

 1. Launch the MATLAB Coder app and follow the instructions to add the foo.m function as the entry-point function.

 2. Use the Autodefine Input Types to specify the inputs types as double(1x3).

 3. In the Generate Code panel, set the parameters as shown:
    Build type — Dynamic Library
    Langauge — C++
    Interface style — Functions
    Hardware Board — None - Select device below
    Device vendor — Google
    Device type — V8 Engine
    Toolchain — Emscripten v2.0.26 | gmake (64-bit Windows)

 4. In More Settigns > All Settings, find and set the *Post-code-generation command* parameter to: `wasm.coder.postcodegen.addCLinkage(buildInfo); wasm.coder.postcodegen.registerExportedFunctions(buildInfo)`. This adds C linkage to C++ functions to prevent name mangling and registers the entrypoint functions to be exported and accessible to the JavaScript environment.

  5. Click the GENERATE button to start the code generation process. When the build process completes, the output files, foo.js and foo.wasm, get generated into the codegen/dll/foo directory.

  6. Move the output files to the root of the matlab folder:

```shell
cp codegen/dll/*/*.js codegen/dll/*/*.wasm .
```