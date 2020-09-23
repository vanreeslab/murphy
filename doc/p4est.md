## Notes on P4EST library

--------------------------------------------
### Miscellaneous notes
- the quadrants are indexed using the type `p4est_locidx_t`.
- the constant `P8EST_MAXLEVEL` defines the number of max level, ranging from `0` to `P8EST_QMAXLEVEL`.

--------------------------------------------
### Ghosts and Mirrors
#### Definitions
- **ghost** the ghosting is seen as an external layer of 1 quadrant around the local quadrant space.
For any procs, this list of quadrant within the external layer is called `ghost`.
- **mirrors** from the other side, one block which is in the ghost list of someone else is called a `mirror`.

#### Mirrors
The access to the mirrors is done through several layers of arrays as one block might be a ghost of several other ranks.
Within the ghost structure, the array
1. `mirrors` stores the list of the quadrant that are a mirror. The qadrants are not full but they `piggy3` structure is filled:
    - `piggy3.local_num`: the local number accross all the trees
    - `piggy3.which_tree`: the tree number.
2. `mirror_proc_mirrors` stores an index of each mirrors, ordered by destination rank. E.g. considering 5 mirrors: `0 2 3 4 2 0 1 3`.
3. `mirror_proc_offset` is a cummulative array that stores the starting and ending position in the array `mirros_proc_mirrors` of each destination proc. Still considering the example above, if three ranks are there, and we are rank 0, we might have: `0 0 4 8`, indicating that rank `1` needs mirrors `0`, `2`, `3`, `4` and that rank `2` will get mirror `2`, `0`, `1` and `3`

As an example, here is the code used to access one mirror that need to go to the rank `ir`:

```C++
// get the first and the last index for the current rank
lid_t first  = ghost->mirror_proc_offsets[ir];
lid_t last   = ghost->mirror_proc_offsets[ir + 1];
// for each mirror that needs to be send to rank ir
for(lid_t bid = send_first; bid< send_last; bid++){
    // get the corresponding mirror object
    p8est_quadrant_t* mirror = p8est_quadrant_array_index(&ghost->mirrors, ghost->mirror_proc_mirrors[bid]);
    // rest of the code
    ...
}
```

#### Ghosts
Ghosts are simplier as they do not have duplicates among the other ranks.
Again, we can access them using different arrays:
1. `proc_offsets` stores a cummulative list of the ghost index for a given rank
2. `tree_offsets` stores a cummulative list of the ghost indexes for a given rank (not very clear to me)

Then, a ghost that comes from rank `ir`, can simply be accessed using

```C++
lid_t first = ghost->proc_offsets[ir];
lid_t last  = ghost->proc_offsets[ir + 1];
for (lid_t bid = recv_first; bid < recv_last; bid++) {
    // get the ghost quad
    p8est_quadrant_t* ghost = p8est_quadrant_array_index(&ghost->ghosts, bid);
     // rest of the code
    ...
}
```
From my understanding and the accumulated experience, the `piggy3` member is filled and the corresponding `level` is present as well.
At least, they use that assumption as well in the function `p4est_ghost_exchange_custom_levels_begin` (`p4est_ghost.c`, l. 2499).