#pragma once
//
#include <cmath>
#include <numbers>
#include <type_traits>

namespace GenericMath
{
	template<typename T>
		requires(std::is_floating_point_v<T>)
	constexpr T DegreesToRadians(T angleDeg) { return angleDeg * T(std::numbers::pi / 180.0); }

	template<typename T>
		requires(std::is_floating_point_v<T>)
	constexpr T RadiansToDegrees(T angleRad) { return angleRad * T(180.0 / std::numbers::pi); }

	inline double NormalizeRad(double angle) {
		angle = std::fmod(angle, 2.0 * std::numbers::pi);
		angle = (angle < 0.0) ? angle + 2.0 * std::numbers::pi : angle;
		return (angle >= std::numbers::pi) ? angle - 2.0 * std::numbers::pi : angle;
	}
	inline double NormalizeRadPositive(double angle) {
		angle = std::fmod(angle, 2.0 * std::numbers::pi);
		return angle < 0 ? angle + 2.0 * std::numbers::pi : angle;
	}

	inline double NormalizeDeg(double angle) {
		angle = std::fmod(angle, 360.0);
		angle = (angle < 0.0) ? angle + 360.0 : angle;
		return (angle >= 180.0) ? angle - 360.0 : angle;
	}
	inline double NormalizeDegPositive(double angle) {
		angle = std::fmod(angle, 360.0);
		return angle < 0 ? angle + 360.0 : angle;
	}
}