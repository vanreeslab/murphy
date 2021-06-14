#include "wavelet/interpolating_wavelet.hpp"

//-----------------------------------------------------------------------------
// Wavelet 2.0
template <>
constexpr short_t len_ha_<2, 0> = 1;
template <>
constexpr short_t len_ga_<2, 0> = 3;
template <>
constexpr short_t len_js_<2, 0> = 1;
template <>
constexpr short_t len_ks_<2, 0> = 3;
template <>
constexpr real_t ha_<2, 0>[1] = {1.0};
template <>
constexpr real_t ga_<2, 0>[3] = {-0.5, 1.0, -0.5};
template <>
constexpr real_t js_<2, 0>[1] = {1.0};
template <>
constexpr real_t ks_<2, 0>[3] = {0.5, 1.0, 0.5};
template<>
constexpr real_t eps_c<2,0> = 7.0;

//-----------------------------------------------------------------------------
// Wavelet 2.2
template <>
constexpr short_t len_ha_<2, 2> = 5;
template <>
constexpr short_t len_ga_<2, 2> = 3;
template <>
constexpr short_t len_js_<2, 2> = 3;
template <>
constexpr short_t len_ks_<2, 2> = 5;
template <>
constexpr real_t ha_<2, 2>[5] = {-0.125, 0.25, 0.75, 0.25, -0.125};
template <>
constexpr real_t ga_<2, 2>[3] = {-0.5, 1.0, -0.5};
template <>
constexpr real_t js_<2, 2>[3] = {-1.0 / 4.0, 1.0, -1.0 / 4.0};
template <>
constexpr real_t ks_<2, 2>[5] = {-1.0 / 8.0, 1.0 / 2.0, 3.0 / 4.0, 1.0 / 2.0, -1.0 / 8.0};
template<>
constexpr real_t eps_c<2,2> = 7.0;

//-----------------------------------------------------------------------------
// Wavelet 4.0
template <>
constexpr short_t len_ha_<4, 0> = 1;
template <>
constexpr short_t len_ga_<4, 0> = 7;
template <>
constexpr short_t len_js_<4, 0> = 1;
template <>
constexpr short_t len_ks_<4, 0> = 7;
template <>
constexpr real_t ha_<4, 0>[1] = {1.0};
template <>
constexpr real_t ga_<4, 0>[7] = {1.0 / 16.0, 0.0, -9.0 / 16.0, 1.0, -9.0 / 16.0, 0.0, 1.0 / 16.0};
template <>
constexpr real_t js_<4, 0>[1] = {1.0};
template <>
constexpr real_t ks_<4, 0>[7] = {-1.0 / 16.0, 0.0, 9.0 / 16.0, 1.0, 9.0 / 16.0, 0.0, -1.0 / 16.0};
template <>
constexpr real_t eps_c<4, 0> = 9.4375;

//-----------------------------------------------------------------------------
// Wavelet 4.2
template <>
constexpr short_t len_ha_<4, 2> = 9;
template <>
constexpr short_t len_ga_<4, 2> = 7;
template <>
constexpr short_t len_js_<4, 2> = 3;
template <>
constexpr short_t len_ks_<4, 2> = 9;
template <>
constexpr real_t ha_<4, 2>[9] = {1.0 / 64.0, 0.0, -1.0 / 8.0, 1.0 / 4.0, 23.0 / 32.0, 1.0 / 4.0, -1.0 / 8.0, 0.0, 1.0 / 64.0};
template <>
constexpr real_t ga_<4, 2>[7] = {1.0 / 16.0, 0.0, -9.0 / 16.0, 1.0, -9.0 / 16.0, 0.0, 1.0 / 16.0};
template <>
constexpr real_t js_<4, 2>[3] = {-1.0 / 4.0, 1.0, -1.0 / 4.0};
template <>
constexpr real_t ks_<4, 2>[9] = {1.0 / 64.0, -1.0 / 16.0, -1.0 / 8.0, 9.0 / 16.0, 23.0 / 32.0, 9.0 / 16.0, -1.0 / 8.0, -1.0 / 16.0, 1.0 / 64.0};
template <>
constexpr real_t eps_c<4, 2> = 9.4375;

// //-----------------------------------------------------------------------------
// // Wavelet 4.4
// template <>
// constexpr short_t len_ha_<4, 4> = 13;
// template <>
// constexpr short_t len_ga_<4, 4> = 4;
// template <>
// constexpr short_t len_js_<4, 4> = 7;
// template <>
// constexpr short_t len_ks_<4, 4> = 13;
// template <>
// constexpr real_t ha_<4, 4>[13] = {-1.0 / 512.0, 0.0, 9.0 / 256.0, -1.0 / 32.0, -63.0 / 512.0, 9.0 / 32.0, 87.0 / 128.0, 9.0 / 32.0, -63.0 / 512.0, -1.0 / 32.0, 9.0 / 256.0, 0.0, -1.0 / 512.0};
// template <>
// constexpr real_t ga_<4, 4>[4] = {-0.0625, 0.5625, 0.5625, -0.0625};
// template <>
// constexpr real_t js_<4, 4>[7] = {1.0 / 32.0, 0.0, -9.0 / 32.0, 1.0, -9.0 / 32.0, 0.0, 1.0 / 32.0};
// template <>
// constexpr real_t ks_<4, 4>[13] = {-1.0 / 512.0, 0.0, 9.0 / 256.0, -1.0 / 16.0, -63.0 / 512.0, 9 / 16.0, 87.0 / 128.0, 9.0 / 16.0, -63.0 / 512.0, -1.0 / 16.0, 9.0 / 256.0, 0.0, -1.0 / 512.0};

//-----------------------------------------------------------------------------
// Wavelet 6.0
template <>
constexpr short_t len_ha_<6, 0> = 1;
template <>
constexpr short_t len_ga_<6, 0> = 11;
template <>
constexpr short_t len_js_<6, 0> = 1;
template <>
constexpr short_t len_ks_<6, 0> = 11;
template <>
constexpr real_t ha_<6, 0>[1] = {1.0};
template <>
constexpr real_t ga_<6, 0>[11] = {-3.0 / 256.0, 0.0, 25.0 / 256.0, 0.0, -75.0 / 128.0, 1.0, -75.0 / 128.0, 0.0, 25.0 / 256.0, 0.0, -3.0 / 256.0};
template <>
constexpr real_t js_<6, 0>[1] = {1.0};
template <>
constexpr real_t ks_<6, 0>[11] = {3.0 / 256.0, 0.0, -25.0 / 256.0, 0.0, 75.0 / 128.0, 1.0, 75.0 / 128.0, 0.0, -25.0 / 256.0, 0.0, 3.0 / 256.0};
template <>
constexpr real_t eps_c<6, 0> = 10.973388671875;

//-----------------------------------------------------------------------------
// Wavelet 6.2
template <>
constexpr short_t len_ha_<6, 2> = 13;
template <>
constexpr short_t len_ga_<6, 2> = 11;
template <>
constexpr short_t len_js_<6, 2> = 3;
template <>
constexpr short_t len_ks_<6, 2> = 13;
template <>
constexpr real_t ha_<6, 2>[13] = {-3.0 / 1024.0, 0.0, 11.0 / 512.0, 0.0, -125.0 / 1024.0, 1.0 / 4.0, 181.0 / 256.0, 1.0 / 4.0, -125.0 / 1024.0, 0.0, 11.0 / 512.0, 0.0, -3.0 / 1024.0};
template <>
constexpr real_t ga_<6, 2>[11] = {-3.0 / 256.0, 0.0, 25.0 / 256.0, 0.0, -75.0 / 128.0, 1.0, -75.0 / 128.0, 0.0, 25.0 / 256.0, 0.0, -3.0 / 256.0};
template <>
constexpr real_t js_<6, 2>[3] = {-1.0 / 4.0, 1.0, -1.0 / 4.0};
template <>
constexpr real_t ks_<6, 2>[13] = {-3.0 / 1024.0, 3.0 / 256.0, 11.0 / 512.0, -25.0 / 256.0, -125.0 / 1024.0, 75.0 / 128.0, 181.0 / 256.0, 75.0 / 128.0, -125.0 / 1024.0, -25.0 / 256.0, 11.0 / 512.0, 3.0 / 256.0, -3.0 / 1024.0};
template <>
constexpr real_t eps_c<6, 2> = 10.973388671875;

// //-----------------------------------------------------------------------------
// // Wavelet 6.4
// template <>
// constexpr short_t len_ha_<6, 4> = 17;
// template <>
// constexpr short_t len_ga_<6, 4> = 6;
// template <>
// constexpr short_t len_js_<6, 4> = 7;
// template <>
// constexpr short_t len_ks_<6, 4> = 17;
// template <>
// constexpr real_t ha_<6, 4>[17] = {3.0 / 8192.0, 0.0, -13.0 / 2048.0, 0.0, 87.0 / 2048.0, -1.0 / 32.0, -243.0 / 2048.0, 9.0 / 32.0, 2721.0 / 4096.0, 9.0 / 32.0, -243.0 / 2048.0, -1.0 / 32.0, 87.0 / 2048.0, 0.0, -13.0 / 2048.0, 0.0, 3.0 / 8192.0};
// template <>
// constexpr real_t ga_<6, 4>[6] = {3.0 / 256.0, -25.0 / 256.0, 75.0 / 128.0, 75.0 / 128.0, -25.0 / 256.0, 3.0 / 256.0};
// template <>
// constexpr real_t js_<6, 4>[7] = {1.0 / 32.0, 0.0, -9.0 / 32.0, 1.0, -9.0 / 32.0, 0.0, 1.0 / 32.0};
// template <>
// constexpr real_t ks_<6, 4>[17] = {3.0 / 8192.0, 0.0, -13.0 / 2048.0, 3.0 / 256.0, 87.0 / 2048.0, -25.0 / 256.0, -243.0 / 2048.0, 75.0 / 128.0, 2721.0 / 4096.0, 75.0 / 128.0, -243.0 / 2048.0, -25.0 / 256.0, 87.0 / 2048.0, 3.0 / 256.0, -13.0 / 2048.0, 0.0, 3.0 / 8192.0};
