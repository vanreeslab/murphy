#ifndef SRC_CORE_DATA_HPP
#define SRC_CORE_DATA_HPP

#include "core/types.hpp"
#include "core/macros.hpp"

template <typename T>
class m_Ptr {
    T* ptr_ = std::nullptr;
#ifndef NDEBUG
    size_t size_         = 0;
    bool   is_allocated_ = false;
#endif

    explicit m_Ptr() noexcept {};
    ~m_Ptr() noexcept {
#ifndef NDEBUG
        m_assert(!is_allocated_, "you must free the memory first!");
#endif
    }

    [[nodiscard]] T* ptr() noexcept const {
        return ptr_;
    }
#ifndef NDEBUG
    [[nodiscard]] size_t size() noexcept const {
        return size_;
    }
#endif

    void Allocate(const size_t size) {
        ptr_ = reinterpret_cast<T*>(m_calloc(size * sizeof(T)));
#ifndef NDEBUG
        size_         = size;
        is_allocated_ = true;
#endif
    }

    void Free() {
    }
};

struct m_Layout {
    size_t stride = 0;
    size_t gs[3]  = {0, 0, 0};
};

template <typename T>
class m_Data {
    T*     data_;
    size_t stride_;

    explicit m_Data(const m_Layout& layout, const m_Ptr& ptr) noexcept {
        m_assert((M_ALIGNMENT % sizeof(T)) == 0, "the alignement must be a mutiple of %d", sizeof(T));
        const size_t chunk = M_ALIGNMENT / sizeof(T);
        const size_t mod   = (layout.stride % chunk);

        // get the actual memory stride!
        stride_ = layout.stride - mod + (mod > 0) * chunk;

        // get the sifted pointer
        const size_t offset = layout.gs[0] + stride_ * (layout.gs[1] + layout.stride * layout.gs[2]);
        m_assert((offset % chunk) == 0, "the offset must be a modulo of the chunk!");
        data_ = layout.ptr + offset;
    }
};

#endif