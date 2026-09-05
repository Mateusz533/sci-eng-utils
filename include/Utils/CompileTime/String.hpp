#pragma once
//
#include <array>
#include <string_view>
//
#include "Utils/Types.hpp"

namespace Utils::CompileTime
{
	using namespace Utils::Types;

	template<usize CAPACITY>
		requires(CAPACITY > 0)
	class String
	{
		template<usize N>
			requires(N > 0)
		friend class String;

	public:
		consteval String() = default;
		consteval String(const std::array<char, CAPACITY - 1>& text) {
			for(usize i = 0; i < Size(); ++i) {
				data[i] = text[i];
			}
			data[Size()] = '\0';
		}
		consteval String(const char (&s)[CAPACITY]) {
			for(usize i = 0; i < CAPACITY; ++i) {
				data[i] = s[i];
			}
		}

		template<usize N>
		consteval auto Concat(const String<N>& other) const {
			String<Size() + N> result{};
			for(usize i = 0; i < Size(); ++i) {
				result.data[i] = data[i];
			}
			for(usize i = 0; i < N; ++i) {
				result.data[Size() + i] = other.data[i];
			}
			return result;
		}

		template<usize N, usize... Ns>
		consteval auto Join(const String<N>& first, const String<Ns>&... parts) const {
			if constexpr(sizeof...(parts) == 0) {
				return Concat(first);
			} else {
				return Concat(first).Join(parts...);
			}
		}

		template<usize OFFSET, usize LENGTH>
			requires((OFFSET + LENGTH < CAPACITY))
		consteval auto SubString() const {
			std::array<char, LENGTH> result = {};
			for(usize i = 0; i < result.size(); ++i) {
				result[i] = data[i + OFFSET];
			}
			return String<result.size() + 1>(result);
		}

		consteval std::string_view View() const { return {data.data(), Size()}; }
		consteval const char* CString() const { return data.data(); }
		static consteval usize Size() { return CAPACITY - 1; }

	private:
		std::array<char, CAPACITY> data{};	// necessarily null-terminated
	};

	template<usize VALUE, usize BASE = 10>
		requires(BASE > 1)
	consteval usize DigitCount() {
		usize counter = 1;
		usize remainingValue = VALUE;
		while(remainingValue >= BASE) {
			remainingValue /= BASE;
			++counter;
		}
		return counter;
	}

	template<i64 VALUE, i64 BASE = 10>
		requires(BASE > 1)
	consteval auto IntegerToString() {
		constexpr usize N = DigitCount<((VALUE >= 0) ? VALUE : -VALUE)>();

		std::array<char, (VALUE >= 0) ? N : N + 1> result{};
		i64 remainingValue = (VALUE >= 0) ? VALUE : -VALUE;
		for(i64 i = result.size() - 1; i >= 0; --i) {
			result[i] = static_cast<char>('0' + (remainingValue % BASE));
			remainingValue /= BASE;
		}
		if constexpr(VALUE < 0) {
			result[0] = '-';
		}
		return String<result.size() + 1>(result);
	}
}