/*
Copyright(c) 2021, Tetsuya Takahashi and Christopher Batty

Permission is hereby granted, free of charge, to any person obtaining a copy
of this softwareand associated documentation files(the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and /or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions :

The above copyright noticeand this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/

#ifndef FRACTION_HPP
#define FRACTION_HPP
#include <array>
#include "Eigen/Dense"

namespace Fraction
{

//=================================================================================================
//		
//=================================================================================================

using iType = int;
using fType = double;

using S4 = std::array<fType, 4>;
using S6 = std::array<fType, 6>;
using S8 = std::array<fType, 8>;
using Vec3s = Eigen::Matrix<fType, 3, 1>;

constexpr fType iso_value = 0.0;
const Vec3s v0(0.0, 0.0, 0.0);
const Vec3s v1(1.0, 0.0, 0.0);
const Vec3s v2(1.0, 1.0, 0.0);
const Vec3s v3(0.0, 1.0, 0.0);
const Vec3s v4(0.0, 0.0, 1.0);
const Vec3s v5(1.0, 0.0, 1.0);
const Vec3s v6(1.0, 1.0, 1.0);
const Vec3s v7(0.0, 1.0, 1.0);

template<typename T>
constexpr T square(const T x) { return x * x; }

constexpr fType get_len_frac(fType a, fType b) {
	if (a < 0 && b < 0) return 1.0;
	if (a < 0 && 0 <= b) return a / (a - b);
	if (0 <= a && b < 0) return b / (b - a);
	return 0;
}

//=================================================================================================
//		ms-style
//=================================================================================================

constexpr iType get_ms_table_index(const S4& v)
{
	iType table_index = 0;
	if (v[0] < iso_value) table_index |= 1;
	if (v[1] < iso_value) table_index |= 2;
	if (v[2] < iso_value) table_index |= 4;
	if (v[3] < iso_value) table_index |= 8;
	return table_index;
}

constexpr fType get_ms_area(const S4& v)
{
	const auto get_split_case = [](const S4& v) {
		if (v[0] + v[1] + v[2] + v[3] < 0) {
			return 1.0 - 0.5 * (
				((1.0 - get_len_frac(v[0], v[1])) * (1.0 - get_len_frac(v[1], v[2]))) +
				((1.0 - get_len_frac(v[3], v[0])) * (1.0 - get_len_frac(v[3], v[2]))));
		}
		else {
			return 0.5 * (get_len_frac(v[0], v[1]) * get_len_frac(v[0], v[3]) + get_len_frac(v[2], v[1]) * get_len_frac(v[2], v[3]));
		}
	};
	switch (const iType table_index = get_ms_table_index(v); table_index) {
	case 0: return 0.0;
	case 1: return 0.5 * get_len_frac(v[0], v[1]) * get_len_frac(v[0], v[3]);
	case 2: return 0.5 * get_len_frac(v[1], v[0]) * get_len_frac(v[1], v[2]);
	case 3: return 0.5 * (get_len_frac(v[0], v[3]) + get_len_frac(v[1], v[2]));
	case 4: return 0.5 * get_len_frac(v[2], v[1]) * get_len_frac(v[2], v[3]);
	case 5: return get_split_case(v);
	case 6: return 0.5 * (get_len_frac(v[1], v[0]) + get_len_frac(v[2], v[3]));
	case 7: return 1.0 - 0.5 * (1.0 - get_len_frac(v[0], v[3])) * (1.0 - get_len_frac(v[2], v[3]));
	case 8: return 0.5 * get_len_frac(v[0], v[3]) * get_len_frac(v[2], v[3]);
	case 9: return 0.5 * (get_len_frac(v[0], v[1]) + get_len_frac(v[3], v[2]));
	case 10: return get_split_case({ v[1], v[2], v[3], v[0] });
	case 11: return 1.0 - 0.5 * (1.0 - get_len_frac(v[1], v[2])) * (1.0 - get_len_frac(v[2], v[3]));
	case 12: return 0.5 * (get_len_frac(v[0], v[3]) + get_len_frac(v[2], v[1]));
	case 13: return 1.0 - 0.5 * (1.0 - get_len_frac(v[0], v[1])) * (1.0 - get_len_frac(v[1], v[2]));
	case 14: return 1.0 - 0.5 * (1.0 - get_len_frac(v[0], v[1])) * (1.0 - get_len_frac(v[0], v[3]));
	case 15: return 1.0;
	default: return 0.0;//dummy
	}
}

constexpr fType get_ms_len(const S4& v)
{
	const auto get_split_case = [](const S4& v) {
		if (v[0] + v[1] + v[2] + v[3] < 0) {
			return
				std::sqrt(square(1.0 - get_len_frac(v[0], v[1])) + square(1.0 - get_len_frac(v[1], v[2]))) +
				std::sqrt(square(1.0 - get_len_frac(v[3], v[0])) + square(1.0 - get_len_frac(v[3], v[2])));
		}
		else {
			return
				std::sqrt(square(get_len_frac(v[0], v[1])) + square(get_len_frac(v[0], v[3]))) +
				std::sqrt(square(get_len_frac(v[2], v[1])) + square(get_len_frac(v[2], v[3])));
		}
	};
	switch (const iType table_index = get_ms_table_index(v); table_index) {
	case 0: return 0.0;
	case 1: return std::sqrt(square(get_len_frac(v[0], v[1])) + square(get_len_frac(v[0], v[3])));
	case 2: return std::sqrt(square(get_len_frac(v[1], v[0])) + square(get_len_frac(v[1], v[2])));
	case 3: return std::sqrt(1.0 + square(get_len_frac(v[0], v[3]) - get_len_frac(v[1], v[2])));
	case 4: return std::sqrt(square(get_len_frac(v[2], v[1])) + square(get_len_frac(v[2], v[3])));
	case 5: return get_split_case(v);
	case 6: return std::sqrt(1.0 + square(get_len_frac(v[1], v[0]) - get_len_frac(v[2], v[3])));
	case 7: return std::sqrt(square(1.0 - get_len_frac(v[0], v[3])) + square(1.0 - get_len_frac(v[2], v[3])));
	case 8: return std::sqrt(square(get_len_frac(v[0], v[3])) + square(get_len_frac(v[2], v[3])));
	case 9: return std::sqrt(1.0 + square(get_len_frac(v[0], v[1]) - get_len_frac(v[3], v[2])));
	case 10: return get_split_case({ v[1], v[2], v[3], v[0] });
	case 11: return std::sqrt(square(1.0 - get_len_frac(v[1], v[2])) + square(1.0 - get_len_frac(v[2], v[3])));
	case 12: return std::sqrt(1.0 + square(get_len_frac(v[0], v[3]) - get_len_frac(v[2], v[1])));
	case 13: return std::sqrt(square(1.0 - get_len_frac(v[0], v[1])) + square(1.0 - get_len_frac(v[1], v[2])));
	case 14: return std::sqrt(square(1.0 - get_len_frac(v[0], v[1])) + square(1.0 - get_len_frac(v[0], v[3])));
	case 15: return 0.0;
	default: return 0.0;//dummy
	}
}

//=================================================================================================
//		mc-style
//=================================================================================================

constexpr iType get_mc_table_index(const S8& v)
{
	iType table_index = 0;
	if (v[0] < iso_value) table_index |= 1;
	if (v[1] < iso_value) table_index |= 2;
	if (v[2] < iso_value) table_index |= 4;
	if (v[3] < iso_value) table_index |= 8;
	if (v[4] < iso_value) table_index |= 16;
	if (v[5] < iso_value) table_index |= 32;
	if (v[6] < iso_value) table_index |= 64;
	if (v[7] < iso_value) table_index |= 128;
	return table_index;
}

template<typename T, iType N>
constexpr std::array<T, 8> get_rotated_vals(const std::array<T, 8>& v)
{
	if constexpr (N == 0) return v;
	else if constexpr (N == 1) return { v[4], v[5], v[1], v[0], v[7], v[6], v[2], v[3] };
	else if constexpr (N == 2) return { v[7], v[6], v[5], v[4], v[3], v[2], v[1], v[0] };
	else if constexpr (N == 3) return { v[3], v[2], v[6], v[7], v[0], v[1], v[5], v[4] };
	else if constexpr (N == 4) return { v[4], v[0], v[3], v[7], v[5], v[1], v[2], v[6] };
	else if constexpr (N == 5) return { v[5], v[4], v[7], v[6], v[1], v[0], v[3], v[2] };
	else if constexpr (N == 6) return { v[1], v[5], v[6], v[2], v[0], v[4], v[7], v[3] };
	else if constexpr (N == 7) return { v[3], v[0], v[1], v[2], v[7], v[4], v[5], v[6] };
	else if constexpr (N == 8) return { v[2], v[3], v[0], v[1], v[6], v[7], v[4], v[5] };
	else if constexpr (N == 9) return { v[1], v[2], v[3], v[0], v[5], v[6], v[7], v[4] };
	else if constexpr (N == 10) return { v[0], v[4], v[5], v[1], v[3], v[7], v[6], v[2] };
	else if constexpr (N == 11) return { v[0], v[3], v[7], v[4], v[1], v[2], v[6], v[5] };
	else if constexpr (N == 12) return { v[2], v[1], v[5], v[6], v[3], v[0], v[4], v[7] };
	else if constexpr (N == 13) return { v[5], v[1], v[0], v[4], v[6], v[2], v[3], v[7] };
	else if constexpr (N == 14) return { v[5], v[6], v[2], v[1], v[4], v[7], v[3], v[0] };
	else if constexpr (N == 15) return { v[7], v[3], v[2], v[6], v[4], v[0], v[1], v[5] };
	else if constexpr (N == 16) return { v[2], v[6], v[7], v[3], v[1], v[5], v[4], v[0] };
	else if constexpr (N == 17) return { v[7], v[4], v[0], v[3], v[6], v[5], v[1], v[2] };
	else if constexpr (N == 18) return { v[1], v[0], v[4], v[5], v[2], v[3], v[7], v[6] };
	else if constexpr (N == 19) return { v[3], v[7], v[4], v[0], v[2], v[6], v[5], v[1] };
	else if constexpr (N == 20) return { v[6], v[7], v[3], v[2], v[5], v[4], v[0], v[1] };
	else if constexpr (N == 21) return { v[6], v[2], v[1], v[5], v[7], v[3], v[0], v[4] };
	else if constexpr (N == 22) return { v[4], v[7], v[6], v[5], v[0], v[3], v[2], v[1] };
	else if constexpr (N == 23) return { v[6], v[5], v[4], v[7], v[2], v[1], v[0], v[3] };
	return v; //dummy
}

template<iType N, bool inside>
constexpr Vec3s get_e(const S8& v)
{
	if constexpr (inside) {
		if constexpr (N == 0) return Vec3s(get_len_frac(v[0], v[1]), 0.0, 0.0);
		else if constexpr (N == 1) return Vec3s(1.0, get_len_frac(v[1], v[2]), 0.0);
		else if constexpr (N == 2) return Vec3s(get_len_frac(v[2], v[3]), 1.0, 0.0);
		else if constexpr (N == 3) return Vec3s(0.0, get_len_frac(v[0], v[3]), 0.0);
		else if constexpr (N == 4) return Vec3s(get_len_frac(v[4], v[5]), 0.0, 1.0);
		else if constexpr (N == 5) return Vec3s(1.0, get_len_frac(v[5], v[6]), 1.0);
		else if constexpr (N == 6) return Vec3s(get_len_frac(v[6], v[7]), 1.0, 1.0);
		else if constexpr (N == 7) return Vec3s(0.0, get_len_frac(v[4], v[7]), 1.0);
		else if constexpr (N == 8) return Vec3s(0.0, 0.0, get_len_frac(v[0], v[4]));
		else if constexpr (N == 9) return Vec3s(1.0, 0.0, get_len_frac(v[1], v[5]));
		else if constexpr (N == 10) return Vec3s(1.0, 1.0, get_len_frac(v[2], v[6]));
		else if constexpr (N == 11) return Vec3s(0.0, 1.0, get_len_frac(v[3], v[7]));
	}
	else if constexpr (!inside) {
		if constexpr (N == 0) return Vec3s(1.0 - get_len_frac(v[0], v[1]), 0.0, 0.0);
		else if constexpr (N == 1) return Vec3s(1.0, 1.0 - get_len_frac(v[1], v[2]), 0.0);
		else if constexpr (N == 2) return Vec3s(1.0 - get_len_frac(v[2], v[3]), 1.0, 0.0);
		else if constexpr (N == 3) return Vec3s(0.0, 1.0 - get_len_frac(v[0], v[3]), 0.0);
		else if constexpr (N == 4) return Vec3s(1.0 - get_len_frac(v[4], v[5]), 0.0, 1.0);
		else if constexpr (N == 5) return Vec3s(1.0, 1.0 - get_len_frac(v[5], v[6]), 1.0);
		else if constexpr (N == 6) return Vec3s(1.0 - get_len_frac(v[6], v[7]), 1.0, 1.0);
		else if constexpr (N == 7) return Vec3s(0.0, 1.0 - get_len_frac(v[4], v[7]), 1.0);
		else if constexpr (N == 8) return Vec3s(0.0, 0.0, 1.0 - get_len_frac(v[0], v[4]));
		else if constexpr (N == 9) return Vec3s(1.0, 0.0, 1.0 - get_len_frac(v[1], v[5]));
		else if constexpr (N == 10) return Vec3s(1.0, 1.0, 1.0 - get_len_frac(v[2], v[6]));
		else if constexpr (N == 11) return Vec3s(0.0, 1.0, 1.0 - get_len_frac(v[3], v[7]));
	}
	return v0;//dummy
}

template<iType N, iType M = 0>
constexpr fType get_mc_vol_case(const S8& v)
{
	constexpr fType scale = 1.0 / 6.0;
	const auto sv = [](const Vec3s& v0, const Vec3s& v1, const Vec3s& v2) {//signed vol
		return ((v1 - v0).cross(v0 - v2)).dot(v0);
	};
	const auto pattern1 = [](const S4& v) {
		return get_len_frac(v[0], v[1]) * get_len_frac(v[0], v[2]) * get_len_frac(v[0], v[3]);
	};
	const auto pattern2 = [](const S6& v) {
		const fType e1 = get_len_frac(v[3], v[4]);
		return ((get_len_frac(v[0], v[1]) + e1) * get_len_frac(v[1], v[2]) + (e1 * get_len_frac(v[4], v[5])));
	};
	const auto pattern21 = [](const S4& v) {
		return (1.0 - get_len_frac(v[0], v[1])) * (1.0 - get_len_frac(v[0], v[2])) * (1.0 - get_len_frac(v[0], v[3]));
	};

	if constexpr (N == 0) {
		return 0.0;
	}
	else if constexpr (N == 1) {
		return scale * pattern1({ v[0], v[1], v[3], v[4] });
	}
	else if constexpr (N == 2) {
		return scale * pattern2({ v[3], v[0], v[4], v[2], v[1], v[5] });
	}
	else if constexpr (N == 3) {
		return scale * (pattern1({ v[0], v[1], v[3], v[4] }) + pattern1({ v[5], v[1], v[6], v[4] }));
	}
	else if constexpr (N == 4) {
		return scale * (pattern1({ v[0], v[1], v[3], v[4] }) + pattern1({ v[6], v[7], v[5], v[2] }));
	}
	else if constexpr (N == 5) {
		const fType e3 = get_len_frac(v[0], v[3]);
		const fType vol0 = get_len_frac(v[1], v[5]) *
			(1.0 - 0.5 * (1.0 - get_len_frac(v[0], v[1])) * (1.0 - e3)) / 3.0;
		const fType e11 = get_len_frac(v[3], v[7]);
		return vol0 + scale * ((get_len_frac(v[2], v[6]) + e11) + e3 * e11);
	}
	else if constexpr (N == 6) {
		return scale * (pattern2({ v[3], v[0], v[4], v[2], v[1], v[5] }) + pattern1({ v[6], v[7], v[5], v[2] }));
	}
	else if constexpr (N == 7) {
		return scale * (pattern1({ v[4], v[0], v[5], v[7] }) + pattern1({ v[1], v[0], v[2], v[5] }) +
			pattern1({ v[6], v[7], v[5], v[2] }));
	}
	else if constexpr (N == 8) {
		return scale * (2.0 * get_len_frac(v[0], v[4]) + get_len_frac(v[1], v[5]) +
			2.0 * get_len_frac(v[2], v[6]) + get_len_frac(v[3], v[7]));
	}
	else if constexpr (N == 9) {
		const Vec3s e0 = get_e<0, true>(v), e1 = get_e<1, false>(v), e6 = get_e<6, true>(v),
			e7 = get_e<7, false>(v), e8 = get_e<8, true>(v), e10 = get_e<10, true>(v);
		const fType _e10 = get_len_frac(v[2], v[6]);
		const fType _e6 = get_len_frac(v[6], v[7]);
		const fType area =
			(1.0 - 0.5 * (1.0 - _e6) * (1.0 - _e10)) / 3.0 +
			scale * get_len_frac(v[1], v[2]) * _e10 +
			scale * (_e6 * get_len_frac(v[4], v[7]));
		if constexpr (M == 0) {
			return area + scale * (sv(e10, e7, e6) + sv(e1, e7, e10) + sv(e1, e8, e7) + sv(e1, e0, e8));
		}
		else if constexpr (M == 1) {
			return area + scale * (sv(e10, e1, e0) + sv(e6, e10, e0) + sv(e6, e0, e8) + sv(e6, e8, e7));
		}
	}
	else if constexpr (N == 10) {
		return scale * (pattern2({ v[4], v[7], v[6], v[0], v[3], v[2] }) + pattern2({ v[2], v[1], v[0], v[6], v[5], v[4] }));
	}
	else if constexpr (N == 11) {
		const Vec3s e0 = get_e<0, true>(v), e1 = get_e<1, false>(v), e5 = get_e<5, false>(v),
			e6 = get_e<6, false>(v), e8 = get_e<8, true>(v), e11 = get_e<11, true>(v);
		const fType _e5 = get_len_frac(v[5], v[6]);
		const fType _e6 = get_len_frac(v[6], v[7]);
		const fType area =
			scale * (get_len_frac(v[1], v[2]) + _e5) +
			(1.0 - 0.5 * (1.0 - _e6) * (1.0 - get_len_frac(v[3], v[7]))) / 3.0 +
			scale * (_e5 * _e6);
		return area + scale * (sv(e0, e8, e11) + sv(e0, e11, e5) + sv(e0, e5, e1) + sv(e5, e11, e6));
	}
	else if constexpr (N == 12) {
		const fType e0 = get_len_frac(v[0], v[1]);
		const fType e9 = get_len_frac(v[1], v[5]);
		return get_len_frac(v[3], v[7]) * (1.0 - 0.5 * (1.0 - e0) * (1.0 - get_len_frac(v[0], v[3]))) / 3.0 +
			scale * ((e9 + get_len_frac(v[2], v[6])) + e0 * e9 + pattern1({ v[4], v[0], v[5], v[7] }));
	}
	else if constexpr (N == 13) {
		return scale * (pattern1({ v[4], v[0], v[5], v[7] }) + pattern1({ v[1], v[0], v[2], v[5] }) +
			pattern1({ v[6], v[7], v[5], v[2] }) + pattern1({ v[3], v[0], v[2], v[7] }));
	}
	else if constexpr (N == 14) {
		const Vec3s e0 = get_e<0, false>(v), e3 = get_e<3, false>(v), e6 = get_e<6, true>(v),
			e7 = get_e<7, false>(v), e9 = get_e<9, true>(v), e10 = get_e<10, true>(v);
		const fType _e10 = get_len_frac(v[2], v[6]);
		const fType _e6 = get_len_frac(v[6], v[7]);
		const fType area = scale * (get_len_frac(v[1], v[5]) + _e10) +
			scale * (_e6 * get_len_frac(v[4], v[7])) +
			(1.0 - 0.5 * (1.0 - _e6) * (1.0 - _e10)) / 3.0;
		return area + scale * (sv(e0, e3, e7) + sv(e0, e7, e10) + sv(e0, e10, e9) + sv(e6, e10, e7));
	}
	else if constexpr (N == 15) {
		const Vec3s e0 = get_e<0, true>(v), e1 = get_e<1, false>(v),
			e6 = get_e<6, true>(v), e7 = get_e<7, false>(v),
			e8 = get_e<8, true>(v), e10 = get_e<10, true>(v);
		const fType _e10 = get_len_frac(v[2], v[6]);
		const fType _e6 = get_len_frac(v[6], v[7]);
		return (1.0 - 0.5 * (1.0 - _e6) * (1.0 - _e10)) / 3.0 +
			scale * get_len_frac(v[1], v[2]) * _e10 +
			scale * (_e6 * get_len_frac(v[4], v[7])) +
			scale * (sv(e1, e6, e10) + sv(e1, e7, e6) + sv(e1, e0, e7) + sv(e8, e7, e0)) +
			scale * pattern1({ v[5], v[4], v[1], v[6] });
	}
	else if constexpr (N == 16) {
		const Vec3s e1 = get_e<1, false>(v), e3 = get_e<3, false>(v), e5 = get_e<5, true>(v),
			e6 = get_e<6, true>(v), e8 = get_e<8, false>(v), e9 = get_e<9, false>(v),
			e10 = get_e<10, true>(v);
		const fType _e10 = get_len_frac(v[2], v[6]);
		const fType _e5 = get_len_frac(v[5], v[6]);
		const fType _e6 = get_len_frac(v[6], v[7]);
		return scale * (sv(e1, e3, e6) + sv(e1, e6, e10) + sv(e3, e8, e6) + sv(e5, e6, e9) + sv(e8, e9, e6)) +
			scale * (get_len_frac(v[1], v[2]) * _e10 + _e5 * get_len_frac(v[1], v[5])) +
			(1.0 - 0.5 * (1.0 - _e6) * (1.0 - _e10)) / 3.0 +
			(1.0 - 0.5 * (1.0 - _e6) * (1.0 - _e5)) / 3.0;
	}
	else if constexpr (N == 17) {
		const fType e3 = get_len_frac(v[0], v[3]);
		const fType e11 = get_len_frac(v[3], v[7]);
		return 1.0 - ((1.0 - get_len_frac(v[1], v[5])) * (1.0 - 0.5 * get_len_frac(v[0], v[1]) * e3) / 3.0 +
			scale * ((1.0 - get_len_frac(v[2], v[6])) + (1.0 - e11) + (1.0 - e3) * (1.0 - e11)));
	}
	else if constexpr (N == 18) {
		return 1.0 - scale * (pattern21({ v[0], v[1], v[3], v[4] }) + pattern21({ v[6], v[5], v[2], v[7] }));
	}
	else if constexpr (N == 19) {
		const Vec3s e0 = get_e<0, false>(v), e3 = get_e<3, false>(v), e4 = get_e<4, true>(v),
			e5 = get_e<5, false>(v), e8 = get_e<8, false>(v), e9 = get_e<9, true>(v);
		const fType _e5 = get_len_frac(v[5], v[6]);
		return 1.0 - (scale * (-sv(e8, e4, e5) - sv(e8, e5, e3) - sv(e9, e0, e5) - sv(e0, e3, e5)) +
			(0.5 * (1.0 - get_len_frac(v[1], v[5])) * (1.0 - _e5)) / 3.0 +
			(0.5 * (1.0 - _e5) * (1.0 - get_len_frac(v[4], v[5]))) / 3.0);
	}
	else if constexpr (N == 20) {
		const fType e1 = get_len_frac(v[1], v[2]);
		return 1.0 - scale * (((1.0 - get_len_frac(v[0], v[3])) +
			(1.0 - e1)) * (1.0 - get_len_frac(v[0], v[4])) +
			(1.0 - e1) * (1.0 - get_len_frac(v[1], v[5])));
	}
	else if constexpr (N == 21) {
		return 1.0 - scale * pattern21({ v[0], v[1], v[3], v[4] });
	}
	else if constexpr (N == 22) {
		return 1.0;
	}
	return 0.0;//dummy
}

template<iType N, iType M = 0>
constexpr fType get_mc_area_case(const S8& v)
{
	constexpr fType scale = 0.5;
	const auto a = [](const Vec3s& v0, const Vec3s& v1, const Vec3s& v2) { return ((v1 - v0).cross(v2 - v0)).norm(); };

	if constexpr (N == 0) {
		return 0.0;
	}
	else if constexpr (N == 1) {
		const Vec3s e0 = get_e<0, true>(v), e3 = get_e<3, true>(v), e8 = get_e<8, true>(v);
		return scale * (a(e0, e8, e3));
	}
	else if constexpr (N == 2) {
		const Vec3s e1 = get_e<1, true>(v), e3 = get_e<3, true>(v), e8 = get_e<8, true>(v), e9 = get_e<9, true>(v);
		return scale * (a(e1, e8, e3) + a(e9, e8, e1));
	}
	else if constexpr (N == 3) {
		const Vec3s e0 = get_e<0, true>(v), e3 = get_e<3, true>(v), e4 = get_e<4, false>(v),
			e5 = get_e<5, true>(v), e8 = get_e<8, true>(v), e9 = get_e<9, false>(v);
		return scale * (a(e9, e5, e4) + a(e0, e8, e3));
	}
	else if constexpr (N == 4) {
		const Vec3s e0 = get_e<0, true>(v), e3 = get_e<3, true>(v), e5 = get_e<5, false>(v),
			e6 = get_e<6, false>(v), e8 = get_e<8, true>(v), e10 = get_e<10, false>(v);
		return scale * (a(e0, e8, e3) + a(e5, e10, e6));
	}
	else if constexpr (N == 5) {
		const Vec3s e0 = get_e<0, false>(v), e3 = get_e<3, false>(v), e9 = get_e<9, true>(v),
			e10 = get_e<10, true>(v), e11 = get_e<11, true>(v);
		return scale * (a(e3, e9, e0) + a(e3, e11, e9) + a(e11, e10, e9));
	}
	else if constexpr (N == 6) {
		const Vec3s e1 = get_e<1, true>(v), e3 = get_e<3, true>(v), e5 = get_e<5, false>(v),
			e6 = get_e<6, false>(v), e8 = get_e<8, true>(v), e9 = get_e<9, true>(v), e10 = get_e<10, false>(v);
		return scale * (a(e1, e8, e3) + a(e1, e9, e8) + a(e5, e10, e6));
	}
	else if constexpr (N == 7) {
		const Vec3s e0 = get_e<0, false>(v), e1 = get_e<1, true>(v), e4 = get_e<4, true>(v),
			e5 = get_e<5, false>(v), e6 = get_e<6, false>(v), e7 = get_e<7, true>(v),
			e8 = get_e<8, false>(v), e9 = get_e<9, true>(v), e10 = get_e<10, false>(v);
		return scale * (a(e1, e9, e0) + a(e5, e10, e6) + a(e8, e4, e7));
	}
	else if constexpr (N == 8) {
		const Vec3s e8 = get_e<8, true>(v), e9 = get_e<9, true>(v), e10 = get_e<10, true>(v), e11 = get_e<11, true>(v);
		return scale * (a(e9, e8, e10) + a(e10, e8, e11));
	}
	else if constexpr (N == 9) {
		const Vec3s e0 = get_e<0, true>(v), e1 = get_e<1, false>(v), e6 = get_e<6, true>(v),
			e7 = get_e<7, false>(v), e8 = get_e<8, true>(v), e10 = get_e<10, true>(v);
		if constexpr (M == 0) {
			return scale * (a(e10, e7, e6) + a(e1, e7, e10) + a(e1, e8, e7) + a(e1, e0, e8));
		}
		else if constexpr (M == 1) {
			return scale * (a(e10, e1, e0) + a(e6, e10, e0) + a(e6, e0, e8) + a(e6, e8, e7));
		}
	}
	else if constexpr (N == 10) {
		const Vec3s e0 = get_e<0, false>(v), e1 = get_e<1, true>(v), e2 = get_e<2, true>(v),
			e3 = get_e<3, false>(v), e4 = get_e<4, false>(v), e5 = get_e<5, true>(v),
			e6 = get_e<6, true>(v), e7 = get_e<7, false>(v);
		return scale * (a(e3, e6, e2) + a(e3, e7, e6) + a(e1, e5, e0) + a(e5, e4, e0));
	}
	else if constexpr (N == 11) {
		const Vec3s e0 = get_e<0, true>(v), e1 = get_e<1, false>(v), e5 = get_e<5, false>(v),
			e6 = get_e<6, false>(v), e8 = get_e<8, true>(v), e11 = get_e<11, true>(v);
		return scale * (a(e0, e8, e11) + a(e0, e11, e5) + a(e0, e5, e1) + a(e5, e11, e6));
	}
	else if constexpr (N == 12) {
		const Vec3s e0 = get_e<0, false>(v), e3 = get_e<3, false>(v), e4 = get_e<4, true>(v),
			e7 = get_e<7, true>(v), e8 = get_e<8, false>(v), e9 = get_e<9, true>(v),
			e10 = get_e<10, true>(v), e11 = get_e<11, true>(v);
		return scale * (a(e4, e7, e8) + a(e9, e0, e11) + a(e9, e11, e10) + a(e11, e0, e3));
	}
	else if constexpr (N == 13) {
		const Vec3s e0 = get_e<0, false>(v), e1 = get_e<1, true>(v), e2 = get_e<2, true>(v),
			e3 = get_e<3, false>(v), e4 = get_e<4, true>(v), e5 = get_e<5, false>(v),
			e6 = get_e<6, false>(v), e7 = get_e<7, true>(v), e8 = get_e<8, false>(v),
			e9 = get_e<9, true>(v), e10 = get_e<10, false>(v), e11 = get_e<11, true>(v);
		return scale * (a(e0, e1, e9) + a(e4, e7, e8) + a(e2, e3, e11) + a(e5, e10, e6));
	}
	else if constexpr (N == 14) {
		const Vec3s e0 = get_e<0, false>(v), e3 = get_e<3, false>(v), e6 = get_e<6, true>(v),
			e7 = get_e<7, false>(v), e9 = get_e<9, true>(v), e10 = get_e<10, true>(v);
		return scale * (a(e0, e3, e7) + a(e0, e7, e10) + a(e0, e10, e9) + a(e6, e10, e7));
	}
	else if constexpr (N == 15) {
		const Vec3s e0 = get_e<0, true>(v), e1 = get_e<1, false>(v), e4 = get_e<4, false>(v),
			e5 = get_e<5, true>(v), e6 = get_e<6, true>(v), e7 = get_e<7, false>(v),
			e8 = get_e<8, true>(v), e9 = get_e<9, false>(v), e10 = get_e<10, true>(v);
		return scale * (a(e1, e6, e10) + a(e1, e7, e6) + a(e1, e0, e7) + a(e8, e7, e0) + a(e9, e5, e4));
	}
	else if constexpr (N == 16) {
		const Vec3s e1 = get_e<1, false>(v), e3 = get_e<3, false>(v), e5 = get_e<5, true>(v),
			e6 = get_e<6, true>(v), e8 = get_e<8, false>(v), e9 = get_e<9, false>(v), e10 = get_e<10, true>(v);
		return scale * (a(e1, e3, e6) + a(e1, e6, e10) + a(e3, e8, e6) + a(e5, e6, e9) + a(e8, e9, e6));
	}
	else if constexpr (N == 17) {
		const Vec3s e0 = get_e<0, true>(v), e3 = get_e<3, true>(v), e9 = get_e<9, false>(v),
			e10 = get_e<10, false>(v), e11 = get_e<11, false>(v);
		return scale * (a(e3, e0, e9) + a(e3, e9, e11) + a(e11, e9, e10));
	}
	else if constexpr (N == 18) {
		const Vec3s e0 = get_e<0, false>(v), e3 = get_e<3, false>(v), e5 = get_e<5, true>(v),
			e6 = get_e<6, true>(v), e8 = get_e<8, false>(v), e10 = get_e<10, true>(v);
		return scale * (a(e0, e3, e8) + a(e5, e6, e10));
	}
	else if constexpr (N == 19) {
		const Vec3s e0 = get_e<0, false>(v), e3 = get_e<3, false>(v), e4 = get_e<4, true>(v),
			e5 = get_e<5, false>(v), e8 = get_e<8, false>(v), e9 = get_e<9, true>(v);
		return scale * (a(e8, e4, e5) + a(e8, e5, e3) + a(e9, e0, e5) + a(e0, e3, e5));
	}
	else if constexpr (N == 20) {
		const Vec3s e1 = get_e<1, false>(v), e3 = get_e<3, false>(v), e8 = get_e<8, false>(v), e9 = get_e<9, false>(v);
		return scale * (a(e1, e3, e8) + a(e9, e1, e8));
	}
	else if constexpr (N == 21) {
		const Vec3s e0 = get_e<0, false>(v), e3 = get_e<3, false>(v), e8 = get_e<8, false>(v);
		return scale * (a(e0, e3, e8));
	}
	else if constexpr (N == 22) {
		return 0.0;
	}
	return 0.0;//dummy
}

constexpr fType get_mc_vol(const S8& v)
{
	switch (const iType table_index = get_mc_table_index(v); table_index) {
	case 0: return get_mc_vol_case<0, 0>(get_rotated_vals<fType, 0>(v));
	case 1: return get_mc_vol_case<1, 0>(get_rotated_vals<fType, 0>(v));
	case 2: return get_mc_vol_case<1, 0>(get_rotated_vals<fType, 6>(v));
	case 3: return get_mc_vol_case<2, 0>(get_rotated_vals<fType, 0>(v));
	case 4: return get_mc_vol_case<1, 0>(get_rotated_vals<fType, 8>(v));
	case 5: return get_mc_vol_case<3, 0>(get_rotated_vals<fType, 12>(v));
	case 6: return get_mc_vol_case<2, 0>(get_rotated_vals<fType, 12>(v));
	case 7: return get_mc_vol_case<5, 0>(get_rotated_vals<fType, 7>(v));
	case 8: return get_mc_vol_case<1, 0>(get_rotated_vals<fType, 3>(v));
	case 9: return get_mc_vol_case<2, 0>(get_rotated_vals<fType, 7>(v));
	case 10: return get_mc_vol_case<3, 0>(get_rotated_vals<fType, 3>(v));
	case 11: return get_mc_vol_case<5, 0>(get_rotated_vals<fType, 8>(v));
	case 12: return get_mc_vol_case<2, 0>(get_rotated_vals<fType, 8>(v));
	case 13: return get_mc_vol_case<5, 0>(get_rotated_vals<fType, 9>(v));
	case 14: return get_mc_vol_case<5, 0>(get_rotated_vals<fType, 0>(v));
	case 15: return get_mc_vol_case<8, 0>(get_rotated_vals<fType, 0>(v));
	case 16: return get_mc_vol_case<1, 0>(get_rotated_vals<fType, 1>(v));
	case 17: return get_mc_vol_case<2, 0>(get_rotated_vals<fType, 10>(v));
	case 18: return get_mc_vol_case<3, 0>(get_rotated_vals<fType, 4>(v));
	case 19: return get_mc_vol_case<5, 0>(get_rotated_vals<fType, 13>(v));
	case 20: return get_mc_vol_case<4, 0>(get_rotated_vals<fType, 8>(v));
	case 21: return get_mc_vol_case<6, 0>(get_rotated_vals<fType, 4>(v));
	case 22: return get_mc_vol_case<6, 0>(get_rotated_vals<fType, 12>(v));
	case 23: return get_mc_vol_case<11, 0>(get_rotated_vals<fType, 8>(v));
	case 24: return get_mc_vol_case<3, 0>(get_rotated_vals<fType, 7>(v));
	case 25: return get_mc_vol_case<5, 0>(get_rotated_vals<fType, 17>(v));
	case 26: return get_mc_vol_case<7, 0>(get_rotated_vals<fType, 12>(v));
	case 27: return get_mc_vol_case<9, 1>(get_rotated_vals<fType, 19>(v));
	case 28: return get_mc_vol_case<6, 0>(get_rotated_vals<fType, 8>(v));
	case 29: return get_mc_vol_case<14, 1>(get_rotated_vals<fType, 9>(v));
	case 30: return get_mc_vol_case<12, 0>(get_rotated_vals<fType, 0>(v));
	case 31: return get_mc_vol_case<17, 0>(get_rotated_vals<fType, 22>(v));
	case 32: return get_mc_vol_case<1, 0>(get_rotated_vals<fType, 5>(v));
	case 33: return get_mc_vol_case<3, 0>(get_rotated_vals<fType, 0>(v));
	case 34: return get_mc_vol_case<2, 0>(get_rotated_vals<fType, 13>(v));
	case 35: return get_mc_vol_case<5, 0>(get_rotated_vals<fType, 1>(v));
	case 36: return get_mc_vol_case<3, 0>(get_rotated_vals<fType, 13>(v));
	case 37: return get_mc_vol_case<7, 0>(get_rotated_vals<fType, 4>(v));
	case 38: return get_mc_vol_case<5, 0>(get_rotated_vals<fType, 21>(v));
	case 39: return get_mc_vol_case<9, 1>(get_rotated_vals<fType, 14>(v));
	case 40: return get_mc_vol_case<4, 0>(get_rotated_vals<fType, 5>(v));
	case 41: return get_mc_vol_case<6, 0>(get_rotated_vals<fType, 7>(v));
	case 42: return get_mc_vol_case<6, 0>(get_rotated_vals<fType, 13>(v));
	case 43: return get_mc_vol_case<14, 0>(get_rotated_vals<fType, 8>(v));
	case 44: return get_mc_vol_case<6, 0>(get_rotated_vals<fType, 3>(v));
	case 45: return get_mc_vol_case<12, 1>(get_rotated_vals<fType, 9>(v));
	case 46: return get_mc_vol_case<11, 0>(get_rotated_vals<fType, 14>(v));
	case 47: return get_mc_vol_case<17, 0>(get_rotated_vals<fType, 5>(v));
	case 48: return get_mc_vol_case<2, 0>(get_rotated_vals<fType, 1>(v));
	case 49: return get_mc_vol_case<5, 0>(get_rotated_vals<fType, 18>(v));
	case 50: return get_mc_vol_case<5, 0>(get_rotated_vals<fType, 10>(v));
	case 51: return get_mc_vol_case<8, 0>(get_rotated_vals<fType, 13>(v));
	case 52: return get_mc_vol_case<6, 0>(get_rotated_vals<fType, 1>(v));
	case 53: return get_mc_vol_case<12, 1>(get_rotated_vals<fType, 18>(v));
	case 54: return get_mc_vol_case<14, 0>(get_rotated_vals<fType, 10>(v));
	case 55: return get_mc_vol_case<17, 0>(get_rotated_vals<fType, 16>(v));
	case 56: return get_mc_vol_case<6, 0>(get_rotated_vals<fType, 5>(v));
	case 57: return get_mc_vol_case<11, 0>(get_rotated_vals<fType, 13>(v));
	case 58: return get_mc_vol_case<12, 1>(get_rotated_vals<fType, 10>(v));
	case 59: return get_mc_vol_case<17, 0>(get_rotated_vals<fType, 3>(v));
	case 60: return get_mc_vol_case<10, 0>(get_rotated_vals<fType, 21>(v));
	case 61: return get_mc_vol_case<16, 0>(get_rotated_vals<fType, 2>(v));
	case 62: return get_mc_vol_case<16, 0>(get_rotated_vals<fType, 20>(v));
	case 63: return get_mc_vol_case<20, 0>(get_rotated_vals<fType, 20>(v));
	case 64: return get_mc_vol_case<1, 0>(get_rotated_vals<fType, 20>(v));
	case 65: return get_mc_vol_case<4, 0>(get_rotated_vals<fType, 0>(v));
	case 66: return get_mc_vol_case<3, 0>(get_rotated_vals<fType, 23>(v));
	case 67: return get_mc_vol_case<6, 0>(get_rotated_vals<fType, 0>(v));
	case 68: return get_mc_vol_case<2, 0>(get_rotated_vals<fType, 21>(v));
	case 69: return get_mc_vol_case<6, 0>(get_rotated_vals<fType, 21>(v));
	case 70: return get_mc_vol_case<5, 0>(get_rotated_vals<fType, 14>(v));
	case 71: return get_mc_vol_case<14, 0>(get_rotated_vals<fType, 14>(v));
	case 72: return get_mc_vol_case<3, 0>(get_rotated_vals<fType, 21>(v));
	case 73: return get_mc_vol_case<6, 0>(get_rotated_vals<fType, 11>(v));
	case 74: return get_mc_vol_case<7, 0>(get_rotated_vals<fType, 13>(v));
	case 75: return get_mc_vol_case<12, 1>(get_rotated_vals<fType, 8>(v));
	case 76: return get_mc_vol_case<5, 0>(get_rotated_vals<fType, 15>(v));
	case 77: return get_mc_vol_case<11, 0>(get_rotated_vals<fType, 0>(v));
	case 78: return get_mc_vol_case<9, 1>(get_rotated_vals<fType, 20>(v));
	case 79: return get_mc_vol_case<17, 0>(get_rotated_vals<fType, 23>(v));
	case 80: return get_mc_vol_case<3, 0>(get_rotated_vals<fType, 1>(v));
	case 81: return get_mc_vol_case<6, 0>(get_rotated_vals<fType, 10>(v));
	case 82: return get_mc_vol_case<7, 0>(get_rotated_vals<fType, 0>(v));
	case 83: return get_mc_vol_case<12, 0>(get_rotated_vals<fType, 13>(v));
	case 84: return get_mc_vol_case<6, 0>(get_rotated_vals<fType, 16>(v));
	case 85: return get_mc_vol_case<10, 0>(get_rotated_vals<fType, 7>(v));
	case 86: return get_mc_vol_case<12, 1>(get_rotated_vals<fType, 14>(v));
	case 87: return get_mc_vol_case<16, 0>(get_rotated_vals<fType, 19>(v));
	case 88: return get_mc_vol_case<7, 0>(get_rotated_vals<fType, 8>(v));
	case 89: return get_mc_vol_case<12, 0>(get_rotated_vals<fType, 17>(v));
	case 90: return get_mc_vol_case<13, 0>(get_rotated_vals<fType, 0>(v));
	case 91: return get_mc_vol_case<15, 0>(get_rotated_vals<fType, 1>(v));
	case 92: return get_mc_vol_case<12, 0>(get_rotated_vals<fType, 15>(v));
	case 93: return get_mc_vol_case<16, 0>(get_rotated_vals<fType, 6>(v));
	case 94: return get_mc_vol_case<15, 0>(get_rotated_vals<fType, 7>(v));
	case 95: return get_mc_vol_case<19, 0>(get_rotated_vals<fType, 17>(v));
	case 96: return get_mc_vol_case<2, 0>(get_rotated_vals<fType, 14>(v));
	case 97: return get_mc_vol_case<6, 0>(get_rotated_vals<fType, 23>(v));
	case 98: return get_mc_vol_case<5, 0>(get_rotated_vals<fType, 12>(v));
	case 99: return get_mc_vol_case<11, 0>(get_rotated_vals<fType, 10>(v));
	case 100: return get_mc_vol_case<5, 0>(get_rotated_vals<fType, 6>(v));
	case 101: return get_mc_vol_case<12, 1>(get_rotated_vals<fType, 6>(v));
	case 102: return get_mc_vol_case<8, 0>(get_rotated_vals<fType, 12>(v));
	case 103: return get_mc_vol_case<17, 0>(get_rotated_vals<fType, 11>(v));
	case 104: return get_mc_vol_case<6, 0>(get_rotated_vals<fType, 14>(v));
	case 105: return get_mc_vol_case<10, 0>(get_rotated_vals<fType, 1>(v));
	case 106: return get_mc_vol_case<12, 0>(get_rotated_vals<fType, 12>(v));
	case 107: return get_mc_vol_case<16, 0>(get_rotated_vals<fType, 22>(v));
	case 108: return get_mc_vol_case<14, 1>(get_rotated_vals<fType, 6>(v));
	case 109: return get_mc_vol_case<16, 0>(get_rotated_vals<fType, 17>(v));
	case 110: return get_mc_vol_case<17, 0>(get_rotated_vals<fType, 19>(v));
	case 111: return get_mc_vol_case<20, 0>(get_rotated_vals<fType, 22>(v));
	case 112: return get_mc_vol_case<5, 0>(get_rotated_vals<fType, 2>(v));
	case 113: return get_mc_vol_case<14, 1>(get_rotated_vals<fType, 18>(v));
	case 114: return get_mc_vol_case<9, 1>(get_rotated_vals<fType, 22>(v));
	case 115: return get_mc_vol_case<17, 0>(get_rotated_vals<fType, 20>(v));
	case 116: return get_mc_vol_case<11, 0>(get_rotated_vals<fType, 12>(v));
	case 117: return get_mc_vol_case<16, 0>(get_rotated_vals<fType, 15>(v));
	case 118: return get_mc_vol_case<17, 0>(get_rotated_vals<fType, 4>(v));
	case 119: return get_mc_vol_case<20, 0>(get_rotated_vals<fType, 19>(v));
	case 120: return get_mc_vol_case<12, 0>(get_rotated_vals<fType, 2>(v));
	case 121: return get_mc_vol_case<16, 0>(get_rotated_vals<fType, 9>(v));
	case 122: return get_mc_vol_case<15, 0>(get_rotated_vals<fType, 21>(v));
	case 123: return get_mc_vol_case<19, 0>(get_rotated_vals<fType, 2>(v));
	case 124: return get_mc_vol_case<16, 0>(get_rotated_vals<fType, 18>(v));
	case 125: return get_mc_vol_case<18, 0>(get_rotated_vals<fType, 6>(v));
	case 126: return get_mc_vol_case<19, 0>(get_rotated_vals<fType, 15>(v));
	case 127: return get_mc_vol_case<21, 0>(get_rotated_vals<fType, 2>(v));
	case 128: return get_mc_vol_case<1, 0>(get_rotated_vals<fType, 2>(v));
	case 129: return get_mc_vol_case<3, 0>(get_rotated_vals<fType, 15>(v));
	case 130: return get_mc_vol_case<4, 0>(get_rotated_vals<fType, 6>(v));
	case 131: return get_mc_vol_case<6, 0>(get_rotated_vals<fType, 18>(v));
	case 132: return get_mc_vol_case<3, 0>(get_rotated_vals<fType, 2>(v));
	case 133: return get_mc_vol_case<7, 0>(get_rotated_vals<fType, 21>(v));
	case 134: return get_mc_vol_case<6, 0>(get_rotated_vals<fType, 9>(v));
	case 135: return get_mc_vol_case<12, 1>(get_rotated_vals<fType, 7>(v));
	case 136: return get_mc_vol_case<2, 0>(get_rotated_vals<fType, 19>(v));
	case 137: return get_mc_vol_case<5, 0>(get_rotated_vals<fType, 4>(v));
	case 138: return get_mc_vol_case<6, 0>(get_rotated_vals<fType, 15>(v));
	case 139: return get_mc_vol_case<11, 1>(get_rotated_vals<fType, 9>(v));
	case 140: return get_mc_vol_case<5, 0>(get_rotated_vals<fType, 20>(v));
	case 141: return get_mc_vol_case<9, 0>(get_rotated_vals<fType, 0>(v));
	case 142: return get_mc_vol_case<14, 0>(get_rotated_vals<fType, 0>(v));
	case 143: return get_mc_vol_case<17, 0>(get_rotated_vals<fType, 2>(v));
	case 144: return get_mc_vol_case<2, 0>(get_rotated_vals<fType, 22>(v));
	case 145: return get_mc_vol_case<5, 0>(get_rotated_vals<fType, 19>(v));
	case 146: return get_mc_vol_case<6, 0>(get_rotated_vals<fType, 17>(v));
	case 147: return get_mc_vol_case<14, 0>(get_rotated_vals<fType, 13>(v));
	case 148: return get_mc_vol_case<6, 0>(get_rotated_vals<fType, 22>(v));
	case 149: return get_mc_vol_case<12, 1>(get_rotated_vals<fType, 19>(v));
	case 150: return get_mc_vol_case<10, 0>(get_rotated_vals<fType, 13>(v));
	case 151: return get_mc_vol_case<16, 0>(get_rotated_vals<fType, 14>(v));
	case 152: return get_mc_vol_case<5, 0>(get_rotated_vals<fType, 11>(v));
	case 153: return get_mc_vol_case<8, 0>(get_rotated_vals<fType, 4>(v));
	case 154: return get_mc_vol_case<12, 0>(get_rotated_vals<fType, 11>(v));
	case 155: return get_mc_vol_case<17, 0>(get_rotated_vals<fType, 6>(v));
	case 156: return get_mc_vol_case<11, 1>(get_rotated_vals<fType, 4>(v));
	case 157: return get_mc_vol_case<17, 0>(get_rotated_vals<fType, 12>(v));
	case 158: return get_mc_vol_case<16, 0>(get_rotated_vals<fType, 23>(v));
	case 159: return get_mc_vol_case<20, 0>(get_rotated_vals<fType, 14>(v));
	case 160: return get_mc_vol_case<3, 0>(get_rotated_vals<fType, 17>(v));
	case 161: return get_mc_vol_case<7, 0>(get_rotated_vals<fType, 7>(v));
	case 162: return get_mc_vol_case<6, 0>(get_rotated_vals<fType, 6>(v));
	case 163: return get_mc_vol_case<12, 1>(get_rotated_vals<fType, 1>(v));
	case 164: return get_mc_vol_case<7, 0>(get_rotated_vals<fType, 1>(v));
	case 165: return get_mc_vol_case<13, 0>(get_rotated_vals<fType, 19>(v));
	case 166: return get_mc_vol_case<12, 1>(get_rotated_vals<fType, 21>(v));
	case 167: return get_mc_vol_case<15, 0>(get_rotated_vals<fType, 8>(v));
	case 168: return get_mc_vol_case<6, 0>(get_rotated_vals<fType, 19>(v));
	case 169: return get_mc_vol_case<12, 0>(get_rotated_vals<fType, 4>(v));
	case 170: return get_mc_vol_case<10, 0>(get_rotated_vals<fType, 0>(v));
	case 171: return get_mc_vol_case<16, 0>(get_rotated_vals<fType, 16>(v));
	case 172: return get_mc_vol_case<12, 1>(get_rotated_vals<fType, 20>(v));
	case 173: return get_mc_vol_case<15, 0>(get_rotated_vals<fType, 0>(v));
	case 174: return get_mc_vol_case<16, 0>(get_rotated_vals<fType, 10>(v));
	case 175: return get_mc_vol_case<19, 0>(get_rotated_vals<fType, 1>(v));
	case 176: return get_mc_vol_case<5, 0>(get_rotated_vals<fType, 23>(v));
	case 177: return get_mc_vol_case<9, 0>(get_rotated_vals<fType, 13>(v));
	case 178: return get_mc_vol_case<11, 1>(get_rotated_vals<fType, 18>(v));
	case 179: return get_mc_vol_case<17, 0>(get_rotated_vals<fType, 15>(v));
	case 180: return get_mc_vol_case<12, 0>(get_rotated_vals<fType, 23>(v));
	case 181: return get_mc_vol_case<15, 0>(get_rotated_vals<fType, 13>(v));
	case 182: return get_mc_vol_case<16, 0>(get_rotated_vals<fType, 11>(v));
	case 183: return get_mc_vol_case<19, 0>(get_rotated_vals<fType, 21>(v));
	case 184: return get_mc_vol_case<14, 1>(get_rotated_vals<fType, 23>(v));
	case 185: return get_mc_vol_case<17, 0>(get_rotated_vals<fType, 14>(v));
	case 186: return get_mc_vol_case<16, 0>(get_rotated_vals<fType, 21>(v));
	case 187: return get_mc_vol_case<20, 0>(get_rotated_vals<fType, 21>(v));
	case 188: return get_mc_vol_case<16, 0>(get_rotated_vals<fType, 0>(v));
	case 189: return get_mc_vol_case<19, 0>(get_rotated_vals<fType, 23>(v));
	case 190: return get_mc_vol_case<18, 0>(get_rotated_vals<fType, 0>(v));
	case 191: return get_mc_vol_case<21, 0>(get_rotated_vals<fType, 20>(v));
	case 192: return get_mc_vol_case<2, 0>(get_rotated_vals<fType, 20>(v));
	case 193: return get_mc_vol_case<6, 0>(get_rotated_vals<fType, 20>(v));
	case 194: return get_mc_vol_case<6, 0>(get_rotated_vals<fType, 2>(v));
	case 195: return get_mc_vol_case<10, 0>(get_rotated_vals<fType, 12>(v));
	case 196: return get_mc_vol_case<5, 0>(get_rotated_vals<fType, 3>(v));
	case 197: return get_mc_vol_case<12, 0>(get_rotated_vals<fType, 3>(v));
	case 198: return get_mc_vol_case<11, 1>(get_rotated_vals<fType, 6>(v));
	case 199: return get_mc_vol_case<16, 0>(get_rotated_vals<fType, 5>(v));
	case 200: return get_mc_vol_case<5, 0>(get_rotated_vals<fType, 16>(v));
	case 201: return get_mc_vol_case<14, 1>(get_rotated_vals<fType, 4>(v));
	case 202: return get_mc_vol_case<12, 0>(get_rotated_vals<fType, 16>(v));
	case 203: return get_mc_vol_case<16, 0>(get_rotated_vals<fType, 1>(v));
	case 204: return get_mc_vol_case<8, 0>(get_rotated_vals<fType, 3>(v));
	case 205: return get_mc_vol_case<17, 0>(get_rotated_vals<fType, 10>(v));
	case 206: return get_mc_vol_case<17, 0>(get_rotated_vals<fType, 18>(v));
	case 207: return get_mc_vol_case<20, 0>(get_rotated_vals<fType, 1>(v));
	case 208: return get_mc_vol_case<5, 0>(get_rotated_vals<fType, 5>(v));
	case 209: return get_mc_vol_case<11, 1>(get_rotated_vals<fType, 23>(v));
	case 210: return get_mc_vol_case<12, 0>(get_rotated_vals<fType, 5>(v));
	case 211: return get_mc_vol_case<16, 0>(get_rotated_vals<fType, 3>(v));
	case 212: return get_mc_vol_case<14, 1>(get_rotated_vals<fType, 3>(v));
	case 213: return get_mc_vol_case<16, 0>(get_rotated_vals<fType, 13>(v));
	case 214: return get_mc_vol_case<16, 0>(get_rotated_vals<fType, 7>(v));
	case 215: return get_mc_vol_case<18, 0>(get_rotated_vals<fType, 5>(v));
	case 216: return get_mc_vol_case<9, 0>(get_rotated_vals<fType, 4>(v));
	case 217: return get_mc_vol_case<17, 0>(get_rotated_vals<fType, 21>(v));
	case 218: return get_mc_vol_case<15, 0>(get_rotated_vals<fType, 4>(v));
	case 219: return get_mc_vol_case<19, 0>(get_rotated_vals<fType, 13>(v));
	case 220: return get_mc_vol_case<17, 0>(get_rotated_vals<fType, 1>(v));
	case 221: return get_mc_vol_case<20, 0>(get_rotated_vals<fType, 13>(v));
	case 222: return get_mc_vol_case<19, 0>(get_rotated_vals<fType, 0>(v));
	case 223: return get_mc_vol_case<21, 0>(get_rotated_vals<fType, 5>(v));
	case 224: return get_mc_vol_case<5, 0>(get_rotated_vals<fType, 22>(v));
	case 225: return get_mc_vol_case<12, 1>(get_rotated_vals<fType, 22>(v));
	case 226: return get_mc_vol_case<14, 0>(get_rotated_vals<fType, 12>(v));
	case 227: return get_mc_vol_case<16, 0>(get_rotated_vals<fType, 8>(v));
	case 228: return get_mc_vol_case<9, 0>(get_rotated_vals<fType, 12>(v));
	case 229: return get_mc_vol_case<15, 0>(get_rotated_vals<fType, 12>(v));
	case 230: return get_mc_vol_case<17, 0>(get_rotated_vals<fType, 17>(v));
	case 231: return get_mc_vol_case<19, 0>(get_rotated_vals<fType, 7>(v));
	case 232: return get_mc_vol_case<11, 1>(get_rotated_vals<fType, 3>(v));
	case 233: return get_mc_vol_case<16, 0>(get_rotated_vals<fType, 12>(v));
	case 234: return get_mc_vol_case<16, 0>(get_rotated_vals<fType, 4>(v));
	case 235: return get_mc_vol_case<18, 0>(get_rotated_vals<fType, 8>(v));
	case 236: return get_mc_vol_case<17, 0>(get_rotated_vals<fType, 13>(v));
	case 237: return get_mc_vol_case<19, 0>(get_rotated_vals<fType, 4>(v));
	case 238: return get_mc_vol_case<20, 0>(get_rotated_vals<fType, 10>(v));
	case 239: return get_mc_vol_case<21, 0>(get_rotated_vals<fType, 1>(v));
	case 240: return get_mc_vol_case<8, 0>(get_rotated_vals<fType, 23>(v));
	case 241: return get_mc_vol_case<17, 0>(get_rotated_vals<fType, 0>(v));
	case 242: return get_mc_vol_case<17, 0>(get_rotated_vals<fType, 9>(v));
	case 243: return get_mc_vol_case<20, 0>(get_rotated_vals<fType, 8>(v));
	case 244: return get_mc_vol_case<17, 0>(get_rotated_vals<fType, 8>(v));
	case 245: return get_mc_vol_case<19, 0>(get_rotated_vals<fType, 3>(v));
	case 246: return get_mc_vol_case<20, 0>(get_rotated_vals<fType, 7>(v));
	case 247: return get_mc_vol_case<21, 0>(get_rotated_vals<fType, 3>(v));
	case 248: return get_mc_vol_case<17, 0>(get_rotated_vals<fType, 7>(v));
	case 249: return get_mc_vol_case<20, 0>(get_rotated_vals<fType, 12>(v));
	case 250: return get_mc_vol_case<19, 0>(get_rotated_vals<fType, 12>(v));
	case 251: return get_mc_vol_case<21, 0>(get_rotated_vals<fType, 8>(v));
	case 252: return get_mc_vol_case<20, 0>(get_rotated_vals<fType, 0>(v));
	case 253: return get_mc_vol_case<21, 0>(get_rotated_vals<fType, 6>(v));
	case 254: return get_mc_vol_case<21, 0>(get_rotated_vals<fType, 0>(v));
	case 255: return get_mc_vol_case<22, 0>(get_rotated_vals<fType, 0>(v));
	default: return 0.0;//dummy
	}
}

constexpr fType get_mc_area(const S8& v)
{
	switch (const iType table_index = get_mc_table_index(v); table_index) {
	case 0: return get_mc_area_case<0, 0>(get_rotated_vals<fType, 0>(v));
	case 1: return get_mc_area_case<1, 0>(get_rotated_vals<fType, 0>(v));
	case 2: return get_mc_area_case<1, 0>(get_rotated_vals<fType, 6>(v));
	case 3: return get_mc_area_case<2, 0>(get_rotated_vals<fType, 0>(v));
	case 4: return get_mc_area_case<1, 0>(get_rotated_vals<fType, 8>(v));
	case 5: return get_mc_area_case<3, 0>(get_rotated_vals<fType, 12>(v));
	case 6: return get_mc_area_case<2, 0>(get_rotated_vals<fType, 12>(v));
	case 7: return get_mc_area_case<5, 0>(get_rotated_vals<fType, 7>(v));
	case 8: return get_mc_area_case<1, 0>(get_rotated_vals<fType, 3>(v));
	case 9: return get_mc_area_case<2, 0>(get_rotated_vals<fType, 7>(v));
	case 10: return get_mc_area_case<3, 0>(get_rotated_vals<fType, 3>(v));
	case 11: return get_mc_area_case<5, 0>(get_rotated_vals<fType, 8>(v));
	case 12: return get_mc_area_case<2, 0>(get_rotated_vals<fType, 8>(v));
	case 13: return get_mc_area_case<5, 0>(get_rotated_vals<fType, 9>(v));
	case 14: return get_mc_area_case<5, 0>(get_rotated_vals<fType, 0>(v));
	case 15: return get_mc_area_case<8, 0>(get_rotated_vals<fType, 0>(v));
	case 16: return get_mc_area_case<1, 0>(get_rotated_vals<fType, 1>(v));
	case 17: return get_mc_area_case<2, 0>(get_rotated_vals<fType, 10>(v));
	case 18: return get_mc_area_case<3, 0>(get_rotated_vals<fType, 4>(v));
	case 19: return get_mc_area_case<5, 0>(get_rotated_vals<fType, 13>(v));
	case 20: return get_mc_area_case<4, 0>(get_rotated_vals<fType, 8>(v));
	case 21: return get_mc_area_case<6, 0>(get_rotated_vals<fType, 4>(v));
	case 22: return get_mc_area_case<6, 0>(get_rotated_vals<fType, 12>(v));
	case 23: return get_mc_area_case<11, 0>(get_rotated_vals<fType, 8>(v));
	case 24: return get_mc_area_case<3, 0>(get_rotated_vals<fType, 7>(v));
	case 25: return get_mc_area_case<5, 0>(get_rotated_vals<fType, 17>(v));
	case 26: return get_mc_area_case<7, 0>(get_rotated_vals<fType, 12>(v));
	case 27: return get_mc_area_case<9, 1>(get_rotated_vals<fType, 19>(v));
	case 28: return get_mc_area_case<6, 0>(get_rotated_vals<fType, 8>(v));
	case 29: return get_mc_area_case<14, 1>(get_rotated_vals<fType, 9>(v));
	case 30: return get_mc_area_case<12, 0>(get_rotated_vals<fType, 0>(v));
	case 31: return get_mc_area_case<17, 0>(get_rotated_vals<fType, 22>(v));
	case 32: return get_mc_area_case<1, 0>(get_rotated_vals<fType, 5>(v));
	case 33: return get_mc_area_case<3, 0>(get_rotated_vals<fType, 0>(v));
	case 34: return get_mc_area_case<2, 0>(get_rotated_vals<fType, 13>(v));
	case 35: return get_mc_area_case<5, 0>(get_rotated_vals<fType, 1>(v));
	case 36: return get_mc_area_case<3, 0>(get_rotated_vals<fType, 13>(v));
	case 37: return get_mc_area_case<7, 0>(get_rotated_vals<fType, 4>(v));
	case 38: return get_mc_area_case<5, 0>(get_rotated_vals<fType, 21>(v));
	case 39: return get_mc_area_case<9, 1>(get_rotated_vals<fType, 14>(v));
	case 40: return get_mc_area_case<4, 0>(get_rotated_vals<fType, 5>(v));
	case 41: return get_mc_area_case<6, 0>(get_rotated_vals<fType, 7>(v));
	case 42: return get_mc_area_case<6, 0>(get_rotated_vals<fType, 13>(v));
	case 43: return get_mc_area_case<14, 0>(get_rotated_vals<fType, 8>(v));
	case 44: return get_mc_area_case<6, 0>(get_rotated_vals<fType, 3>(v));
	case 45: return get_mc_area_case<12, 1>(get_rotated_vals<fType, 9>(v));
	case 46: return get_mc_area_case<11, 0>(get_rotated_vals<fType, 14>(v));
	case 47: return get_mc_area_case<17, 0>(get_rotated_vals<fType, 5>(v));
	case 48: return get_mc_area_case<2, 0>(get_rotated_vals<fType, 1>(v));
	case 49: return get_mc_area_case<5, 0>(get_rotated_vals<fType, 18>(v));
	case 50: return get_mc_area_case<5, 0>(get_rotated_vals<fType, 10>(v));
	case 51: return get_mc_area_case<8, 0>(get_rotated_vals<fType, 13>(v));
	case 52: return get_mc_area_case<6, 0>(get_rotated_vals<fType, 1>(v));
	case 53: return get_mc_area_case<12, 1>(get_rotated_vals<fType, 18>(v));
	case 54: return get_mc_area_case<14, 0>(get_rotated_vals<fType, 10>(v));
	case 55: return get_mc_area_case<17, 0>(get_rotated_vals<fType, 16>(v));
	case 56: return get_mc_area_case<6, 0>(get_rotated_vals<fType, 5>(v));
	case 57: return get_mc_area_case<11, 0>(get_rotated_vals<fType, 13>(v));
	case 58: return get_mc_area_case<12, 1>(get_rotated_vals<fType, 10>(v));
	case 59: return get_mc_area_case<17, 0>(get_rotated_vals<fType, 3>(v));
	case 60: return get_mc_area_case<10, 0>(get_rotated_vals<fType, 21>(v));
	case 61: return get_mc_area_case<16, 0>(get_rotated_vals<fType, 2>(v));
	case 62: return get_mc_area_case<16, 0>(get_rotated_vals<fType, 20>(v));
	case 63: return get_mc_area_case<20, 0>(get_rotated_vals<fType, 20>(v));
	case 64: return get_mc_area_case<1, 0>(get_rotated_vals<fType, 20>(v));
	case 65: return get_mc_area_case<4, 0>(get_rotated_vals<fType, 0>(v));
	case 66: return get_mc_area_case<3, 0>(get_rotated_vals<fType, 23>(v));
	case 67: return get_mc_area_case<6, 0>(get_rotated_vals<fType, 0>(v));
	case 68: return get_mc_area_case<2, 0>(get_rotated_vals<fType, 21>(v));
	case 69: return get_mc_area_case<6, 0>(get_rotated_vals<fType, 21>(v));
	case 70: return get_mc_area_case<5, 0>(get_rotated_vals<fType, 14>(v));
	case 71: return get_mc_area_case<14, 0>(get_rotated_vals<fType, 14>(v));
	case 72: return get_mc_area_case<3, 0>(get_rotated_vals<fType, 21>(v));
	case 73: return get_mc_area_case<6, 0>(get_rotated_vals<fType, 11>(v));
	case 74: return get_mc_area_case<7, 0>(get_rotated_vals<fType, 13>(v));
	case 75: return get_mc_area_case<12, 1>(get_rotated_vals<fType, 8>(v));
	case 76: return get_mc_area_case<5, 0>(get_rotated_vals<fType, 15>(v));
	case 77: return get_mc_area_case<11, 0>(get_rotated_vals<fType, 0>(v));
	case 78: return get_mc_area_case<9, 1>(get_rotated_vals<fType, 20>(v));
	case 79: return get_mc_area_case<17, 0>(get_rotated_vals<fType, 23>(v));
	case 80: return get_mc_area_case<3, 0>(get_rotated_vals<fType, 1>(v));
	case 81: return get_mc_area_case<6, 0>(get_rotated_vals<fType, 10>(v));
	case 82: return get_mc_area_case<7, 0>(get_rotated_vals<fType, 0>(v));
	case 83: return get_mc_area_case<12, 0>(get_rotated_vals<fType, 13>(v));
	case 84: return get_mc_area_case<6, 0>(get_rotated_vals<fType, 16>(v));
	case 85: return get_mc_area_case<10, 0>(get_rotated_vals<fType, 7>(v));
	case 86: return get_mc_area_case<12, 1>(get_rotated_vals<fType, 14>(v));
	case 87: return get_mc_area_case<16, 0>(get_rotated_vals<fType, 19>(v));
	case 88: return get_mc_area_case<7, 0>(get_rotated_vals<fType, 8>(v));
	case 89: return get_mc_area_case<12, 0>(get_rotated_vals<fType, 17>(v));
	case 90: return get_mc_area_case<13, 0>(get_rotated_vals<fType, 0>(v));
	case 91: return get_mc_area_case<15, 0>(get_rotated_vals<fType, 1>(v));
	case 92: return get_mc_area_case<12, 0>(get_rotated_vals<fType, 15>(v));
	case 93: return get_mc_area_case<16, 0>(get_rotated_vals<fType, 6>(v));
	case 94: return get_mc_area_case<15, 0>(get_rotated_vals<fType, 7>(v));
	case 95: return get_mc_area_case<19, 0>(get_rotated_vals<fType, 17>(v));
	case 96: return get_mc_area_case<2, 0>(get_rotated_vals<fType, 14>(v));
	case 97: return get_mc_area_case<6, 0>(get_rotated_vals<fType, 23>(v));
	case 98: return get_mc_area_case<5, 0>(get_rotated_vals<fType, 12>(v));
	case 99: return get_mc_area_case<11, 0>(get_rotated_vals<fType, 10>(v));
	case 100: return get_mc_area_case<5, 0>(get_rotated_vals<fType, 6>(v));
	case 101: return get_mc_area_case<12, 1>(get_rotated_vals<fType, 6>(v));
	case 102: return get_mc_area_case<8, 0>(get_rotated_vals<fType, 12>(v));
	case 103: return get_mc_area_case<17, 0>(get_rotated_vals<fType, 11>(v));
	case 104: return get_mc_area_case<6, 0>(get_rotated_vals<fType, 14>(v));
	case 105: return get_mc_area_case<10, 0>(get_rotated_vals<fType, 1>(v));
	case 106: return get_mc_area_case<12, 0>(get_rotated_vals<fType, 12>(v));
	case 107: return get_mc_area_case<16, 0>(get_rotated_vals<fType, 22>(v));
	case 108: return get_mc_area_case<14, 1>(get_rotated_vals<fType, 6>(v));
	case 109: return get_mc_area_case<16, 0>(get_rotated_vals<fType, 17>(v));
	case 110: return get_mc_area_case<17, 0>(get_rotated_vals<fType, 19>(v));
	case 111: return get_mc_area_case<20, 0>(get_rotated_vals<fType, 22>(v));
	case 112: return get_mc_area_case<5, 0>(get_rotated_vals<fType, 2>(v));
	case 113: return get_mc_area_case<14, 1>(get_rotated_vals<fType, 18>(v));
	case 114: return get_mc_area_case<9, 1>(get_rotated_vals<fType, 22>(v));
	case 115: return get_mc_area_case<17, 0>(get_rotated_vals<fType, 20>(v));
	case 116: return get_mc_area_case<11, 0>(get_rotated_vals<fType, 12>(v));
	case 117: return get_mc_area_case<16, 0>(get_rotated_vals<fType, 15>(v));
	case 118: return get_mc_area_case<17, 0>(get_rotated_vals<fType, 4>(v));
	case 119: return get_mc_area_case<20, 0>(get_rotated_vals<fType, 19>(v));
	case 120: return get_mc_area_case<12, 0>(get_rotated_vals<fType, 2>(v));
	case 121: return get_mc_area_case<16, 0>(get_rotated_vals<fType, 9>(v));
	case 122: return get_mc_area_case<15, 0>(get_rotated_vals<fType, 21>(v));
	case 123: return get_mc_area_case<19, 0>(get_rotated_vals<fType, 2>(v));
	case 124: return get_mc_area_case<16, 0>(get_rotated_vals<fType, 18>(v));
	case 125: return get_mc_area_case<18, 0>(get_rotated_vals<fType, 6>(v));
	case 126: return get_mc_area_case<19, 0>(get_rotated_vals<fType, 15>(v));
	case 127: return get_mc_area_case<21, 0>(get_rotated_vals<fType, 2>(v));
	case 128: return get_mc_area_case<1, 0>(get_rotated_vals<fType, 2>(v));
	case 129: return get_mc_area_case<3, 0>(get_rotated_vals<fType, 15>(v));
	case 130: return get_mc_area_case<4, 0>(get_rotated_vals<fType, 6>(v));
	case 131: return get_mc_area_case<6, 0>(get_rotated_vals<fType, 18>(v));
	case 132: return get_mc_area_case<3, 0>(get_rotated_vals<fType, 2>(v));
	case 133: return get_mc_area_case<7, 0>(get_rotated_vals<fType, 21>(v));
	case 134: return get_mc_area_case<6, 0>(get_rotated_vals<fType, 9>(v));
	case 135: return get_mc_area_case<12, 1>(get_rotated_vals<fType, 7>(v));
	case 136: return get_mc_area_case<2, 0>(get_rotated_vals<fType, 19>(v));
	case 137: return get_mc_area_case<5, 0>(get_rotated_vals<fType, 4>(v));
	case 138: return get_mc_area_case<6, 0>(get_rotated_vals<fType, 15>(v));
	case 139: return get_mc_area_case<11, 1>(get_rotated_vals<fType, 9>(v));
	case 140: return get_mc_area_case<5, 0>(get_rotated_vals<fType, 20>(v));
	case 141: return get_mc_area_case<9, 0>(get_rotated_vals<fType, 0>(v));
	case 142: return get_mc_area_case<14, 0>(get_rotated_vals<fType, 0>(v));
	case 143: return get_mc_area_case<17, 0>(get_rotated_vals<fType, 2>(v));
	case 144: return get_mc_area_case<2, 0>(get_rotated_vals<fType, 22>(v));
	case 145: return get_mc_area_case<5, 0>(get_rotated_vals<fType, 19>(v));
	case 146: return get_mc_area_case<6, 0>(get_rotated_vals<fType, 17>(v));
	case 147: return get_mc_area_case<14, 0>(get_rotated_vals<fType, 13>(v));
	case 148: return get_mc_area_case<6, 0>(get_rotated_vals<fType, 22>(v));
	case 149: return get_mc_area_case<12, 1>(get_rotated_vals<fType, 19>(v));
	case 150: return get_mc_area_case<10, 0>(get_rotated_vals<fType, 13>(v));
	case 151: return get_mc_area_case<16, 0>(get_rotated_vals<fType, 14>(v));
	case 152: return get_mc_area_case<5, 0>(get_rotated_vals<fType, 11>(v));
	case 153: return get_mc_area_case<8, 0>(get_rotated_vals<fType, 4>(v));
	case 154: return get_mc_area_case<12, 0>(get_rotated_vals<fType, 11>(v));
	case 155: return get_mc_area_case<17, 0>(get_rotated_vals<fType, 6>(v));
	case 156: return get_mc_area_case<11, 1>(get_rotated_vals<fType, 4>(v));
	case 157: return get_mc_area_case<17, 0>(get_rotated_vals<fType, 12>(v));
	case 158: return get_mc_area_case<16, 0>(get_rotated_vals<fType, 23>(v));
	case 159: return get_mc_area_case<20, 0>(get_rotated_vals<fType, 14>(v));
	case 160: return get_mc_area_case<3, 0>(get_rotated_vals<fType, 17>(v));
	case 161: return get_mc_area_case<7, 0>(get_rotated_vals<fType, 7>(v));
	case 162: return get_mc_area_case<6, 0>(get_rotated_vals<fType, 6>(v));
	case 163: return get_mc_area_case<12, 1>(get_rotated_vals<fType, 1>(v));
	case 164: return get_mc_area_case<7, 0>(get_rotated_vals<fType, 1>(v));
	case 165: return get_mc_area_case<13, 0>(get_rotated_vals<fType, 19>(v));
	case 166: return get_mc_area_case<12, 1>(get_rotated_vals<fType, 21>(v));
	case 167: return get_mc_area_case<15, 0>(get_rotated_vals<fType, 8>(v));
	case 168: return get_mc_area_case<6, 0>(get_rotated_vals<fType, 19>(v));
	case 169: return get_mc_area_case<12, 0>(get_rotated_vals<fType, 4>(v));
	case 170: return get_mc_area_case<10, 0>(get_rotated_vals<fType, 0>(v));
	case 171: return get_mc_area_case<16, 0>(get_rotated_vals<fType, 16>(v));
	case 172: return get_mc_area_case<12, 1>(get_rotated_vals<fType, 20>(v));
	case 173: return get_mc_area_case<15, 0>(get_rotated_vals<fType, 0>(v));
	case 174: return get_mc_area_case<16, 0>(get_rotated_vals<fType, 10>(v));
	case 175: return get_mc_area_case<19, 0>(get_rotated_vals<fType, 1>(v));
	case 176: return get_mc_area_case<5, 0>(get_rotated_vals<fType, 23>(v));
	case 177: return get_mc_area_case<9, 0>(get_rotated_vals<fType, 13>(v));
	case 178: return get_mc_area_case<11, 1>(get_rotated_vals<fType, 18>(v));
	case 179: return get_mc_area_case<17, 0>(get_rotated_vals<fType, 15>(v));
	case 180: return get_mc_area_case<12, 0>(get_rotated_vals<fType, 23>(v));
	case 181: return get_mc_area_case<15, 0>(get_rotated_vals<fType, 13>(v));
	case 182: return get_mc_area_case<16, 0>(get_rotated_vals<fType, 11>(v));
	case 183: return get_mc_area_case<19, 0>(get_rotated_vals<fType, 21>(v));
	case 184: return get_mc_area_case<14, 1>(get_rotated_vals<fType, 23>(v));
	case 185: return get_mc_area_case<17, 0>(get_rotated_vals<fType, 14>(v));
	case 186: return get_mc_area_case<16, 0>(get_rotated_vals<fType, 21>(v));
	case 187: return get_mc_area_case<20, 0>(get_rotated_vals<fType, 21>(v));
	case 188: return get_mc_area_case<16, 0>(get_rotated_vals<fType, 0>(v));
	case 189: return get_mc_area_case<19, 0>(get_rotated_vals<fType, 23>(v));
	case 190: return get_mc_area_case<18, 0>(get_rotated_vals<fType, 0>(v));
	case 191: return get_mc_area_case<21, 0>(get_rotated_vals<fType, 20>(v));
	case 192: return get_mc_area_case<2, 0>(get_rotated_vals<fType, 20>(v));
	case 193: return get_mc_area_case<6, 0>(get_rotated_vals<fType, 20>(v));
	case 194: return get_mc_area_case<6, 0>(get_rotated_vals<fType, 2>(v));
	case 195: return get_mc_area_case<10, 0>(get_rotated_vals<fType, 12>(v));
	case 196: return get_mc_area_case<5, 0>(get_rotated_vals<fType, 3>(v));
	case 197: return get_mc_area_case<12, 0>(get_rotated_vals<fType, 3>(v));
	case 198: return get_mc_area_case<11, 1>(get_rotated_vals<fType, 6>(v));
	case 199: return get_mc_area_case<16, 0>(get_rotated_vals<fType, 5>(v));
	case 200: return get_mc_area_case<5, 0>(get_rotated_vals<fType, 16>(v));
	case 201: return get_mc_area_case<14, 1>(get_rotated_vals<fType, 4>(v));
	case 202: return get_mc_area_case<12, 0>(get_rotated_vals<fType, 16>(v));
	case 203: return get_mc_area_case<16, 0>(get_rotated_vals<fType, 1>(v));
	case 204: return get_mc_area_case<8, 0>(get_rotated_vals<fType, 3>(v));
	case 205: return get_mc_area_case<17, 0>(get_rotated_vals<fType, 10>(v));
	case 206: return get_mc_area_case<17, 0>(get_rotated_vals<fType, 18>(v));
	case 207: return get_mc_area_case<20, 0>(get_rotated_vals<fType, 1>(v));
	case 208: return get_mc_area_case<5, 0>(get_rotated_vals<fType, 5>(v));
	case 209: return get_mc_area_case<11, 1>(get_rotated_vals<fType, 23>(v));
	case 210: return get_mc_area_case<12, 0>(get_rotated_vals<fType, 5>(v));
	case 211: return get_mc_area_case<16, 0>(get_rotated_vals<fType, 3>(v));
	case 212: return get_mc_area_case<14, 1>(get_rotated_vals<fType, 3>(v));
	case 213: return get_mc_area_case<16, 0>(get_rotated_vals<fType, 13>(v));
	case 214: return get_mc_area_case<16, 0>(get_rotated_vals<fType, 7>(v));
	case 215: return get_mc_area_case<18, 0>(get_rotated_vals<fType, 5>(v));
	case 216: return get_mc_area_case<9, 0>(get_rotated_vals<fType, 4>(v));
	case 217: return get_mc_area_case<17, 0>(get_rotated_vals<fType, 21>(v));
	case 218: return get_mc_area_case<15, 0>(get_rotated_vals<fType, 4>(v));
	case 219: return get_mc_area_case<19, 0>(get_rotated_vals<fType, 13>(v));
	case 220: return get_mc_area_case<17, 0>(get_rotated_vals<fType, 1>(v));
	case 221: return get_mc_area_case<20, 0>(get_rotated_vals<fType, 13>(v));
	case 222: return get_mc_area_case<19, 0>(get_rotated_vals<fType, 0>(v));
	case 223: return get_mc_area_case<21, 0>(get_rotated_vals<fType, 5>(v));
	case 224: return get_mc_area_case<5, 0>(get_rotated_vals<fType, 22>(v));
	case 225: return get_mc_area_case<12, 1>(get_rotated_vals<fType, 22>(v));
	case 226: return get_mc_area_case<14, 0>(get_rotated_vals<fType, 12>(v));
	case 227: return get_mc_area_case<16, 0>(get_rotated_vals<fType, 8>(v));
	case 228: return get_mc_area_case<9, 0>(get_rotated_vals<fType, 12>(v));
	case 229: return get_mc_area_case<15, 0>(get_rotated_vals<fType, 12>(v));
	case 230: return get_mc_area_case<17, 0>(get_rotated_vals<fType, 17>(v));
	case 231: return get_mc_area_case<19, 0>(get_rotated_vals<fType, 7>(v));
	case 232: return get_mc_area_case<11, 1>(get_rotated_vals<fType, 3>(v));
	case 233: return get_mc_area_case<16, 0>(get_rotated_vals<fType, 12>(v));
	case 234: return get_mc_area_case<16, 0>(get_rotated_vals<fType, 4>(v));
	case 235: return get_mc_area_case<18, 0>(get_rotated_vals<fType, 8>(v));
	case 236: return get_mc_area_case<17, 0>(get_rotated_vals<fType, 13>(v));
	case 237: return get_mc_area_case<19, 0>(get_rotated_vals<fType, 4>(v));
	case 238: return get_mc_area_case<20, 0>(get_rotated_vals<fType, 10>(v));
	case 239: return get_mc_area_case<21, 0>(get_rotated_vals<fType, 1>(v));
	case 240: return get_mc_area_case<8, 0>(get_rotated_vals<fType, 23>(v));
	case 241: return get_mc_area_case<17, 0>(get_rotated_vals<fType, 0>(v));
	case 242: return get_mc_area_case<17, 0>(get_rotated_vals<fType, 9>(v));
	case 243: return get_mc_area_case<20, 0>(get_rotated_vals<fType, 8>(v));
	case 244: return get_mc_area_case<17, 0>(get_rotated_vals<fType, 8>(v));
	case 245: return get_mc_area_case<19, 0>(get_rotated_vals<fType, 3>(v));
	case 246: return get_mc_area_case<20, 0>(get_rotated_vals<fType, 7>(v));
	case 247: return get_mc_area_case<21, 0>(get_rotated_vals<fType, 3>(v));
	case 248: return get_mc_area_case<17, 0>(get_rotated_vals<fType, 7>(v));
	case 249: return get_mc_area_case<20, 0>(get_rotated_vals<fType, 12>(v));
	case 250: return get_mc_area_case<19, 0>(get_rotated_vals<fType, 12>(v));
	case 251: return get_mc_area_case<21, 0>(get_rotated_vals<fType, 8>(v));
	case 252: return get_mc_area_case<20, 0>(get_rotated_vals<fType, 0>(v));
	case 253: return get_mc_area_case<21, 0>(get_rotated_vals<fType, 6>(v));
	case 254: return get_mc_area_case<21, 0>(get_rotated_vals<fType, 0>(v));
	case 255: return get_mc_area_case<22, 0>(get_rotated_vals<fType, 0>(v));
	default: return 0.0;//dummy
	}
}

//=================================================================================================
//		
//=================================================================================================

}

#endif
