# How to contribute?
----------------------------
### Commit messages
Use the [gitmoji guide](https://gitmoji.carloscuesta.me) to describe your commit purpose.

----------------------------
### Typing variables
To ease the remplacement of the doubles into floats and handle the different types of ints, we define 3 types:
- `sid_t`: small ID types, for numbers aimed between `-127` and `127`.
- `lid_t`: local ID types, for every **local** number, aimed between `-2 147 483 648` and `2 147 483 648`. This type does not fit for memory types, use `size_t` instead
- `real_t` and it's pointer `real_p`: stands for floating points numbers (`double` or `float`).

No `int` declarations are used in the code, except for MPI rank-related numbers, which are `int` by the MPI standard.

----------------------------
### Style Guide
To ensure a consistent style accross murphy, we rely on the [Google C++ Style Guide](https://google.github.io/styleguide/cppguide.html#C++_Version).
We extensivelly use the `cpplint` tool to detect style errors, as well,
```
cpplint --filter=-whitespace/line_length,-runtime/printf myfile
```

However, some exceptions are made to the mentionned guide style:
- `typedef` defined types are named: `*_t` except `real_p` which is with `_p` to emphasize the pointer type (e.g. `real_t` for the real datatype, `sid_t` for the small integers, `lid_t` for the local ids, ...)
- the macro definitions acting as functions are name `m_myfunction`. The `m_` states for `Murphy`
- the files are names in lowercase with a `.cpp` and `.hpp` extension, alike the class it contains (e.g. class `Grid` in `grid.cpp`)
- callback functions used to interface with p4est start with `cback_`

