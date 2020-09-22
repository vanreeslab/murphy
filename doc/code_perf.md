# Code perfomance


Here is a few tips I learned that may impact the code performance.
I used the [compiler explorer app](https://godbolt.org/) that can also be run locally!

#### The inline and macro replacement
Since c++11, `constexpr` provide a way to define compile-time **evaluation** (not replacement but real evaluation!).
However, I couldn't find a way to enforce function inlining, even with `constexpr`.
Therefore the macro as kept and I tried to minimize the overlapping with the variables etc:
- renaming and typing every variable used
- unique name prefixed by the macro name `my_macro_` and postfix by `_`.