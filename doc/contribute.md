# How to contribute?
----------------------------
### Git management
We embrace the [Git Flow](https://nvie.com/posts/a-successful-git-branching-model/) approach to maintain a hierarchy and a clean repository.
It means that we have the following branches:
- `master` is the default working branch
- `develop` is the development branch, the source branch for every new development
- `dev-*` are the ongoing development branch
- `fix-*` fix a bug or solve an issue

To create a new development branch, simply do
```bash
# go on the develop branch
git checkout develop
# update your version
git pull
# create a new branch
git checkout -b dev-mynewfeature develop
```
To incorporate your development into the develop branch, use pull requests.
The automatic testing **must** be included and succeed in the code.

Additionnaly, for the commit we use the [gitmoji guide](https://gitmoji.carloscuesta.me) to describe your commit purpose.
It helps to automatically identify the reason of the commit. Here is a small non-exhaustive list

Action | Corresponding emoji
--------|-----------------------------
Solve a regular bug | :bug: `:bug:`
Sovle a critical bug | :ambulance: `:ambulance:`
Add documentation / comments / doxygen | :memo: `:memo:`
Compilation / makefile | :wrench: `:wrench:`
Docker | :whale: `:whale:`
Refactor | :recycle: `:recycle:`

<!-- ----------------------------
### Typing variables
To ease the remplacement of the doubles into floats and handle the different types of ints, we define 3 types:
- `sid_t`: small ID types, for numbers aimed between `-127` and `127`.
- `lid_t`: local ID types, for every **local** number, aimed between `-2 147 483 648` and `2 147 483 648`. This type does not fit for memory types, use `size_t` instead
- `real_t` and it's pointer `real_p`: stands for floating points numbers (`double` or `float`).

No `int` declarations are used in the code, except for MPI rank-related numbers, which are `int` by the MPI standard. -->

Create issues to keep track of the development and also to discuss any question you might have.
You can use the keywords **clos(e/es/ed)**, **resolv(e/es/ed)** and **fix(/es/ed)**


----------------------------
### Style Guide
To ensure a consistent style accross murphy, we rely on the [Google C++ Style Guide](https://google.github.io/styleguide/cppguide.html#C++_Version).
Additionally, we use two tools to help the formating:
- `cpplint` tool to detect style errors,
```bash
cpplint --filter=-whitespace/line_length,-runtime/printf myfile
```
- `clang-format` to help the formating in VSCode. The file `.clang-format` should be automatically detected and loaded in VSCode, giving you access Immediately to the formating behavior. You regenerate that file using 
```bash
# brew install clang-format if needed
clang-format -style='{BasedOnStyle: Google, ColumnLimit: 0, IndentWidth: 4, AlignConsecutiveAssignments: true, AlignConsecutiveDeclarations: true, AlignEscapedNewlines: true, AlignOperands: true }' -dump-config > .clang-format
```

However, some exceptions are made to the mentionned guide style:
- `typedef` defined types are named: `*_t` except `real_p` which is with `_p` to emphasize the pointer type (e.g. `real_t` for the real datatype, `sid_t` for the small integers, `lid_t` for the local ids, ...)
- the macro definitions acting as functions are name `m_myfunction`. The `m_` states for `Murphy`
- the files are names in lowercase with a `.cpp` and `.hpp` extension, alike the class it contains (e.g. class `Grid` in `grid.cpp`)
- callback functions used to interface with p4est start with `cback_`


----------------------------
### When to use what - table

To ease the development of the code, we use a few custom types and macros: here is a summary

|`What`? | When
--------|-----------------------------
| **MACROS**
`m_verb` | displays information in the body of a function or when some diag of the function are not improtant
`m_log` | displays important information at the end of a function, to give info to the user
`m_min` | returns the min of two numbers
`m_max` | returns the max of two numbers
`m_sign` | returns the sign
`m_pos` | position of a 3D point in the domain
`m_pos_relative` | position of a 3D point in the block
`m_begin` | starts every function to have a timer in full verbose mode
`m_end` | ends every function to have a timer in full verbose mode
| **MEMORY**
`m_calloc` | allocate aligned memory
`m_free` | free aligned memory
`m_isaligned` | check the alignment
`m_assume_aligned` | tells the compiler the array is aligned (and check the assertion in debug mode)
`m_blockmemsize` | returns the number of elements in one block
`m_zeroidx` | shifts the memory to the element `(0,0,0)`
| **TYPES**
`rank_t` | a MPI-rank
`lda_t` | a leading dimension numbers (typically 0,1,2)
`iblock_t` | a local block index
`iface_t` | a face id (0 to 26)
`level_t` | a level
`qdrt_t`| a `p8est_quadrant_t` variable
`real_t` | a floating point number
`data_ptr` | a memory adress of the element (0,0,0)
`mem_ptr` | the raw memory address
`const_mem_ptr` | the constant raw memory address

