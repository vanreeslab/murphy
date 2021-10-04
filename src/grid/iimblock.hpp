#ifndef SRC_IIM_IIMBLOCK_HPP_
#define SRC_IIM_IIMBLOCK_HPP_

#include "grid/gridblock.hpp"

class IIMBlock: public GridBlock {
  protected:
    short_t number_of_intersections_ = -1;
    bidx_t*  mind_ = nullptr; // closest affected point with a negative level set value 
    short_t* axis_ = nullptr; // 0, 1, or 2 to indicate an x, y, or z gridline
    short_t* sign_ = nullptr; // -1 or 1 to indicate orientation along gridline

    #ifndef NDEBUG
        bool is_allocated = false;
    #endif

  public:

    IIMBlock(const real_t length, const real_t xyz[3], const sid_t level) : 
        GridBlock(length, xyz, level) {}
    
    void AllocateIIM(short_t npts) {
        #ifndef NDEBUG
            m_assert(!is_allocated, "Attempt to allocate IIM without freeing");
            is_allocated = true;
        #endif
        mind_ = new bidx_t[npts];
        axis_ = new short_t[npts];
        sign_ = new short_t[npts];
        number_of_intersections_ = npts;
    }

    void FreeIIM() {
        #ifndef NDEBUG
            m_assert(is_allocated, "Attempt to free IIM without allocating");
            is_allocated = false;
        #endif
        delete[] mind_;
        delete[] axis_;
        delete[] sign_;
        number_of_intersections_ = -1;
    }

    short_t number_of_intersections() { return number_of_intersections_; }
    bidx_t*  mind() { return mind_; }   
    short_t* axis() { return axis_; }
    short_t* sign() { return sign_; }
};

#endif