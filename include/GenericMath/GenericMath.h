#pragma once
//
#include <cmath>
#include <type_traits>

namespace GenericMath
{
	template<typename T>
		requires(std::is_floating_point_v<T>)
	constexpr T DegreesToRadians(T angleDeg) { return angleDeg * T(M_PI / 180.0); }

	template<typename T>
		requires(std::is_floating_point_v<T>)
	constexpr T RadiansToDegrees(T angleRad) { return angleRad * T(180.0 / M_PI); }

	inline double NormalizeRad(double angle) {
		angle = std::fmod(angle, 2.0 * M_PI);
		angle = (angle < 0.0) ? angle + 2.0 * M_PI : angle;
		return (angle >= M_PI) ? angle - 2.0 * M_PI : angle;
	}
	inline double NormalizeRadPositive(double angle) {
		angle = std::fmod(angle, 2.0 * M_PI);
		return angle < 0 ? angle + 2.0 * M_PI : angle;
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