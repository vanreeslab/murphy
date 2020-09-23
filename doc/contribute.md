# How to contribute?
----------------------------
### Git management
We embrace the [Git Flow](https://nvie.com/posts/a-successful-git-branching-model/) approach to maintain a hierarchy and a clean repository.
It means that we have the following branches:
- `master` is the default working branch
- `develop` is the development branch, the source branch for every new development
- `dev-*` are the ongoing development branch

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
Solve a regular bug | :bug:
Sovle a critical bug | :ambulance:
Add documentation / comments / doxygen | :memo:
Compilation / makefile | :wrench:

<!-- ----------------------------
### Typing variables
To ease the remplacement of the doubles into floats and handle the different types of ints, we define 3 types:
- `sid_t`: small ID types, for numbers aimed between `-127` and `127`.
- `lid_t`: local ID types, for every **local** number, aimed between `-2 147 483 648` and `2 147 483 648`. This type does not fit for memory types, use `size_t` instead
- `real_t` and it's pointer `real_p`: stands for floating points numbers (`double` or `float`).

No `int` declarations are used in the code, except for MPI rank-related numbers, which are `int` by the MPI standard. -->



----------------------------
### When to use what - table

To ease the development of the code, we use a few custom types and macros: here is a summary

`what`? | When
--------|-----------------------------
`m_verb` | to display information in the body of a function or when some diag of the fucntion are not improtant
`m_log` | to display important information at the end of a function, to give info to the user
`rank_t` | to store MPI-rank ids
`lda_` | to store leading dimension numbers (typically 0,1,2)


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
