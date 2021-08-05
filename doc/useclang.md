# How to use clang and VSCode

From a developer point of view, we highly recommend `clang` over other compilers.
The toolset provided is very convenient and the errors are easier to track.

### configure your C/C++ vscode extension

- install the Intelisens extension
- create a new configuration profile (either open `c_cpp_properties.json` or use the UI: `cmd+P` + type "`>C/C++: Edit configurations (UI)`")
- my configuration is something like:
```json
{
    "name": "murphy-container-clang",
    "includePath": [
        "src",
        "src/clients",
        "src/core",
        "src/grid",
        "src/operator",
        "src/poisson",
        "src/time",
        "src/tools",
        "src/wavelet",
        "test/src",
        "/soft/flups/include",
        "/soft/p4est-github/include",
        "/usr/local/hdf5/include",
        "/usr/local/openmpi/include",
        "/usr/local/fftw/include",
        "/soft/googletest/include"
    ],
    "defines": [],
    "compilerPath": "/usr/bin/clang",
    "cStandard": "c17",
    "cppStandard": "c++17",
    "intelliSenseMode": "linux-clang-x64",
    "compileCommands": "${workspaceFolder}/compile_commands.json",
    "browse": {
        "path": [
            "${workspaceFolder}/src/**",
            "${workspaceFolder}/test/src/**"
        ],
        "limitSymbolsToIncludedHeaders": true,
        // "databaseFilename": ""
    }
}
```
- choose the configuration: `cmd+P` + type "`>C/C++: Select a configuration`")

Now the IntelliSense of VSCode is using clang to guide you and generate error/warnings


### configure the extension Clang-Tidy
Clang-Tidy is helpfull to track and detect all the c++ errors (like objects declaration etc).

VSCode has an extension that allows you to visualize the output directly in the editor.
- install the extension if not already done (should be done automatically as part of the `.devcontainer.json` file)
- you need to generate a _compilation database_ for clang (a file `compile_commands.json`). As we now have clang in the docker the easiest is by far to [use clang directly](https://sarcasm.github.io/notes/dev/compilation-database.html#clang)
- we implemented a rule in the makefile to do so: `make compdb` (or `make compdb_full` to include the tests)

### Clang Optimization viewer

To run (and view the optimization) analysis, add the flags `-fsave-optimization-record` to the compilation.

Then, simply use the tool from llvm to create html webpages (run the command directly in the Docker):
```
python3 /usr/lib/llvm-11/share/opt-viewer/opt-viewer.py --output-dir opt_reports build/*.opt.yaml
```

### Clang resources:

- [clang pragma's](https://clang.llvm.org/docs/LanguageExtensions.html#id30)
- [auto-vectorization](https://llvm.org/docs/Vectorizers.html)