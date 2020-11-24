#ifndef SRC_NAVIER_STOKES_HPP_
#define SRC_NAVIER_STOKES_HPP_

#include "murphy.hpp"

static char        ns_doc[] = "Navier-Stokes";
extern struct argp extern_ns_argp; //!< promise the declaration of a struct argp somewhere

class NavierStokes {
    explicit NavierStokes();
};

#endif  // SRC_NAVIER_STOKES_HPP_