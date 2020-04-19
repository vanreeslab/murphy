# Implementation choices and Doxygen documentation

### Documentation and implementation strategies
To get the doxygen documentation, go to the `doc` folder and enter
```
doxygen
```
You will need `dot` to get a visual graphs. It is a part of Graphivz.
You can install it using homebrew: `brew install graphviz`.


### Misc

We define two positions of the memory pointer:
- the absolute position, which is the pointer to the memory allocated
- the 0-pointer, which is the pointer to the effective memory, i.e. withouth any ghost points. # Implementation choices