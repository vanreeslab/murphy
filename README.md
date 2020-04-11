# murphy

Anything that can ~~go wrong~~ **scale strong** will ~~go wrong~~ **scale strong**


#### Possible compilation flag
- ```-DLOG_ALLRANKS``` will enable log on every processor. By default, only the master logs
- ```-DVERBOSE``` enable extended logs
- ```-DNDEBUG``` disable the assertion checks and the debug comments

#### Style Guide
To ensure a consistent style accross murphy, we rely on the [Google C++ Style Guide](https://google.github.io/styleguide/cppguide.html#C++_Version).
We also rely on the `cpplint` tool to detect style errors:
```
cpplint --filter=-whitespace/line_length,-runtime/printf myfile
```

However, some exceptions are made to the mentionned guide style:
- `typedef` defined types are named: `*_t` except `real_p` which is with `_p` to emphasize the pointer type (e.g. `real_t` for the real datatype, `sid_t` for the small integers, `lid_t` for the local ids, ...)
- the macro definitions acting as functions are name `m_myfunction`. The `m_` states for `Murphy`
- the files are names in lowercase with a `.cpp` and `.hpp` extension, alike the class it contains (e.g. class `Grid` in `grid.cpp`)
- Callback functions used to interface with p4est start with `cback_`