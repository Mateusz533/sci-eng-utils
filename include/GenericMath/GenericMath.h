#pragma once

#include <cmath>

namespace GenericMath
{
	template<typename T>
	constexpr inline T DegreesToRadians(T angleDeg) { return angleDeg * M_PI / 180.0; }

	template<typename T>
	constexpr inline T RadiansToDegrees(T angleRad) { return angleRad * 180.0 / M_PI; }

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