#include <array>
#include <concepts>
#include <cstddef>
#include <cstdint>
#include <string_view>
#include <type_traits>

namespace Utils::CompileTime
{
	template<std::size_t CAPACITY>
		requires(CAPACITY > 0)
	class String
	{
		template<std::size_t N>
			requires(N > 0)
		friend class String;

	public:
		consteval String() = default;
		consteval String(const std::array<char, CAPACITY - 1>& text) {
			for(std::size_t i = 0; i < Size(); ++i) {
				data[i] = text[i];
			}
			data[Size()] = '\0';
		}
		consteval String(const char (&s)[CAPACITY]) {
			for(std::size_t i = 0; i < CAPACITY; ++i) {
				data[i] = s[i];
			}
		}

		template<std::size_t N>
		consteval auto Concat(const String<N>& other) const {
			String<Size() + N> result{};
			for(std::size_t i = 0; i < Size(); ++i) {
				result.data[i] = data[i];
			}
			for(std::size_t i = 0; i < N; ++i) {
				result.data[Size() + i] = other.data[i];
			}
			return result;
		}

		template<std::size_t N, std::size_t... Ns>
		consteval auto Join(const String<N>& first, const String<Ns>&... parts) const {
			if constexpr(sizeof...(parts) == 0) {
				return Concat(first);
			} else {
				return Concat(first).Join(parts...);
			}
		}

		template<std::size_t OFFSET, std::size_t LENGTH>
			requires((OFFSET + LENGTH < CAPACITY))
		consteval auto SubString() const {
			std::array<char, LENGTH> result = {};
			for(std::size_t i = 0; i < result.size(); ++i) {
				result[i] = data[i + OFFSET];
			}
			return String<result.size() + 1>(result);
		}

		consteval std::string_view View() const { return {data.data(), Size()}; }
		consteval const char* CString() const { return data.data(); }
		static consteval std::size_t Size() { return CAPACITY - 1; }

	private:
		std::array<char, CAPACITY> data{};	// necessarily null-terminated
	};

	template<std::size_t VALUE, std::size_t BASE = 10>
		requires(BASE > 1)
	consteval std::size_t DigitCount() {
		std::size_t counter = 1;
		std::size_t remainingValue = VALUE;
		while(remainingValue >= BASE) {
			remainingValue /= BASE;
			++counter;
		}
		return counter;
	}

	template<long VALUE, long BASE = 10>
		requires(BASE > 1)
	consteval auto IntegerToString() {
		constexpr std::size_t N = DigitCount<((VALUE >= 0) ? VALUE : -VALUE)>();

		std::array<char, (VALUE >= 0) ? N : N + 1> result{};
		long remainingValue = (VALUE >= 0) ? VALUE : -VALUE;
		for(long i = result.size() - 1; i >= 0; --i) {
			result[i] = static_cast<char>('0' + (remainingValue % BASE));
			remainingValue /= BASE;
		}
		if constexpr(VALUE < 0) {
			result[0] = '-';
		}
		return String<result.size() + 1>(result);
	}

	template<typename T, std::uint8_t N>
		requires(std::is_integral_v<T> || std::is_floating_point_v<T>)
	static constexpr T Pow(T base) noexcept {
		if constexpr(N > 0) {
			return base * Pow<T, N - 1U>(base);
		} else {
			return T{1};
		}
	}

	static consteval std::int64_t Exp10(std::uint8_t exp) noexcept {
		return (exp > 0) ? 10L * Exp10(exp - 1) : 1L;
	}

	template<std::floating_point T, std::uint8_t N>
		requires(N > 0)
	static constexpr T RootNewtonRaphson(T x, T curr, T prev) noexcept {
		return (curr == prev) ? curr : RootNewtonRaphson<T, N>(x, (T{N - 1U} * curr + x / Pow<T, N - 1U>(curr)) / T{N}, curr);
	}

	template<std::floating_point T>
	static constexpr T SqrtNewtonRaphson(T x, T curr, T prev) noexcept {
		return RootNewtonRaphson<T, 2>(x, curr, prev);
	}

	template<std::floating_point T>
	static constexpr T CbrtNewtonRaphson(T x, T curr, T prev) noexcept {
		return RootNewtonRaphson<T, 3>(x, curr, prev);
	}
}