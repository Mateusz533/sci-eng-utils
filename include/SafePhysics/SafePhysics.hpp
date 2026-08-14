#pragma once
//
#include <concepts>
#include <cstddef>
#include <cstdint>

namespace Physics
{
	using usize = std::size_t;
	using isize = std::ptrdiff_t;
	using u8 = std::uint8_t;
	using u16 = std::uint16_t;
	using u32 = std::uint32_t;
	using u64 = std::uint64_t;
	using i8 = std::int8_t;
	using i16 = std::int16_t;
	using i32 = std::int32_t;
	using i64 = std::int64_t;
	using f32 = float;
	using f64 = double;
	using f128 = long double;

	static_assert(sizeof(f32) == 4);
	static_assert(sizeof(f64) == 8);
	static_assert(sizeof(f128) == 16);

	template<typename T>
	concept Arithmetic = std::integral<T> || std::floating_point<T>;
}