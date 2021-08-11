# Debug murphy

Here is a non-exhaustive list of the debugging tips

### Use `-fsanitize=address`

Use the address sanitizer from `clang` by adding the `-fsanitize=address` to the compilation and linker flags.

To get a decent stack trace you should compile with `-O1` and `-fno-omit-frame-pointer`. To further improve the trace, use `-fno-optimize-sibling-calls`. 

For more information, refer to the [clang documentation](https://clang.llvm.org/docs/AddressSanitizer.html)

### Use `-fsanitize=undefined`
Use the undefined behavior sanitizer from `clang` by adding the `-fsanitize=undefined` to the compilation and linker flags.

To get a decent stack trace you should compile with `-O1` and `-fno-omit-frame-pointer`. Also, run with `UBSAN_OPTIONS=print_stacktrace=1`.

For more information, refer to the [clang documentation](https://clang.llvm.org/docs/AddressSanitizer.html)

### Use `valgrind`
The classical vagrind is usefull but sometimes not enough (see the sanitize). To run it:

```bash
mpirun -n X valgrind ./murphy...
```

You can also add some options like `--show-leak-kinds=all` `--leak-check=full` `--track-origins=yes` `--error-limit=no` to help detect the leaks.