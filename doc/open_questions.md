# Open questions

---------------------------
### [closed :wastebasket:] Operators 
The definition of so many operators is difficult to follow.
You have operators, then implementation of operators, it that really usefull? 
On one side, it helps to decrease the lines of code, ease the writting of new operators. On the other side, it's a bit messy to have all those operators defined independently.

##### :+1: hosen solution
I reduced the useless inheritance to the operators by use the `DoOp_F_` function as much as possible.

---------------------------
### [closed :wastebasket:] Operators that require ghost points
The ghost object is owned by the grid, so the grid has to handle it.
On the other side, those operators need the ghost information and due to communication and computation overlap, we need to split all the operations: send the ghost, compute the inner operation, receive, compute the outer operation.

##### possible solution: the Grid solution

This overlapping implementation can be done in the grid, as a similar approach is used in the GhostPull implementation. However, this does not fit with the operator overall approach. as we want to deleguate every implementation for the opertor side.


##### chosen solution: the Operator solution

We implemented another appraoch, we define yet another operator.
We can't implement it in the `operator.hpp` file, because the operators only rely on `ForestGrid` and not on the full `Grid`. So, they don't have the knowledge of what is a ghost, etc. This choice has been made to avoid include loops and because operators where intended to be virtual classes.

I then implemented yet another operator, which inheritates from `OperatorF2F` and which implements this overlap.
His name is `Stencil`. This would remove the ghost overlapping from the grid, which is okay. I don't want to implement a ghost operator as the user doesn't have to know about this overlapping approach. At worst, the user can call the anny stencil operator with a `nullptr` as target field.

So, we have twice the same implementation (once for the `GhostPull` function, another one for the `Stencil` class), this not very general and I am not very happy about it.

##### :+1: hosen solution
I chose to build a `Stencil` class which defines the overlapping and to inheritate it with templated laplacian function. This can be used to implement a rotational or any other stencil operator