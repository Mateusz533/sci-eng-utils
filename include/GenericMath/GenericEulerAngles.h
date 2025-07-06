#pragma once

#include "GenericMath.h"
#include "GenericQuaternion.h"

namespace GenericMath
{
	enum class RotationAxis : std::uint8_t {
		X = 0,
		Y = 1,
		Z = 2,
	};

	template<typename T, RotationAxis A0, RotationAxis A1, RotationAxis A2>
		requires(std::is_floating_point_v<T>)
	class EulerAngles
	{
		static_assert(A0 != A1 && A1 != A2, "Invalid axis order!");
		static_assert(A0 != A2, "Not implemented yet!");

	public:
		constexpr EulerAngles() = delete;
		constexpr EulerAngles(T rad0, T rad1, T rad2) : mAnglesRad{rad0, rad1, rad2} {}
		constexpr EulerAngles(const std::array<T, 3>& anglesRad) : mAnglesRad{anglesRad} {}

		constexpr T& operator()(size_t index) { return mAnglesRad[index]; }
		constexpr const T& operator()(size_t index) const { return mAnglesRad[index]; }

		static constexpr EulerAngles FromMatrix(const Matrix3<T>& mat) {
			const auto angle1 = std::asin(-mat((int)A0, (int)A2));
			const bool noSingularity = std::cos(angle1) > SINGULARITY_SCOPE_RAD;
			const auto angle0 = noSingularity ? std::atan2(+mat((int)A1, (int)A2), mat((int)A2, (int)A2))
											  : std::atan2(-mat((int)A2, (int)A1), mat((int)A1, (int)A1));	// or `0.0`
			const auto angle2 = noSingularity ? std::atan2(+mat((int)A0, (int)A1), mat((int)A0, (int)A0))
											  : 0.0;  // or `std::atan2(-mat((int)A1, (int)A0), mat((int)A1, (int)A1))`

			return EulerAngles{DIR * angle0, DIR * angle1, DIR * angle2};
		}
		static constexpr EulerAngles FromQuaternion(const Quaternion<T>& q) {
			const auto& v = q.Vector();
			const auto &v0 = v((int)A0), &v1 = v((int)A1), &v2 = v((int)A2);

			const auto angle1 = std::asin(std::clamp(2.0 * (q.w() * v1 - DIR * v0 * v2), -1.0, 1.0));
			const bool noSingularity = std::cos(angle1) > SINGULARITY_SCOPE_RAD;
			const auto angle0 = noSingularity ? std::atan2(2.0 * (q.w() * v0 + DIR * v1 * v2), 1.0 - 2.0 * (v1 * v1 + v0 * v0))
											  : 2.0 * std::atan2(v0, q.w());  // or `0.0`
			const auto angle2 = noSingularity ? std::atan2(2.0 * (q.w() * v2 + DIR * v1 * v0), 1.0 - 2.0 * (v1 * v1 + v2 * v2))
											  : 0.0;  // or `2.0 * std::atan2(v2, q.w())`

			return EulerAngles{angle0, angle1, angle2};
		}

	public:
		constexpr Matrix3<T> ToMatrix() const {
			const auto s0 = DIR * std::sin(mAnglesRad[0]), c0 = std::cos(mAnglesRad[0]);
			const auto s1 = DIR * std::sin(mAnglesRad[1]), c1 = std::cos(mAnglesRad[1]);
			const auto s2 = DIR * std::sin(mAnglesRad[2]), c2 = std::cos(mAnglesRad[2]);

			Matrix3<T> result;

			result((int)A0, (int)A0) = c1 * c2;
			result((int)A0, (int)A1) = c1 * s2;
			result((int)A0, (int)A2) = -s1;

			result((int)A1, (int)A0) = s0 * s1 * c2 - c0 * s2;
			result((int)A1, (int)A1) = s0 * s1 * s2 + c0 * c2;
			result((int)A1, (int)A2) = s0 * c1;

			result((int)A2, (int)A0) = c0 * s1 * c2 + s0 * s2;
			result((int)A2, (int)A1) = c0 * s1 * s2 - s0 * c2;
			result((int)A2, (int)A2) = c0 * c1;

			return result;
		}
		constexpr Quaternion<T> ToQuaternion() const {
			const auto angle0 = NormalizeRad(mAnglesRad[0]) / 2;
			const auto angle1 = NormalizeRad(mAnglesRad[1]) / 2;
			const auto angle2 = NormalizeRad(mAnglesRad[2]) / 2;

			const auto s0 = std::sin(angle0), c0 = std::cos(angle0);
			const auto s1 = std::sin(angle1), c1 = std::cos(angle1);
			const auto s2 = std::sin(angle2), c2 = std::cos(angle2);

			Quaternion<T> result;

			result.Scalar() /*    */ = c0 * c1 * c2 + DIR * s0 * s1 * s2;
			result.Vector()((int)A0) = s0 * c1 * c2 - DIR * c0 * s1 * s2;
			result.Vector()((int)A1) = c0 * s1 * c2 + DIR * s0 * c1 * s2;
			result.Vector()((int)A2) = c0 * c1 * s2 - DIR * s0 * s1 * c2;

			return result;
		}

		constexpr bool AreForwardAngles() const {
			return std::cos(mAnglesRad[1]) >= 0;
		}
		constexpr bool AreBackwardAngles() const {
			return !AreForwardAngles();
		}
		constexpr EulerAngles CalcForwardAngles() const {
			return AreBackwardAngles() ? CalcDualAngles() : *this;
		}
		constexpr EulerAngles CalcBackwardAngles() const {
			return AreForwardAngles() ? CalcDualAngles() : *this;
		}
		constexpr EulerAngles CalcDualAngles() const {
			return EulerAngles{mAnglesRad[0] - M_PI, M_PI - mAnglesRad[1], mAnglesRad[2] - M_PI};
		}
		constexpr EulerAngles Normalize() const {
			return EulerAngles{NormalizeRad(mAnglesRad[0]), NormalizeRad(mAnglesRad[1]), NormalizeRad(mAnglesRad[2])};
		}

	private:
		static inline constexpr T DIR = ((3u + (uint)A0 - (uint)A1) % 3u == 1u) ? T(1) : T(-1);
		static inline constexpr T SINGULARITY_SCOPE_RAD = 1e-6;
		std::array<T, 3> mAnglesRad;
	};

	template<typename T>
	using EulerAnglesLocalXYZ = EulerAngles<T, RotationAxis::X, RotationAxis::Y, RotationAxis::Z>;
	template<typename T>
	using EulerAnglesLocalXZY = EulerAngles<T, RotationAxis::X, RotationAxis::Z, RotationAxis::Y>;
	template<typename T>
	using EulerAnglesLocalYZX = EulerAngles<T, RotationAxis::Y, RotationAxis::Z, RotationAxis::X>;
	template<typename T>
	using EulerAnglesLocalYXZ = EulerAngles<T, RotationAxis::Y, RotationAxis::X, RotationAxis::Z>;
	template<typename T>
	using EulerAnglesLocalZXY = EulerAngles<T, RotationAxis::Z, RotationAxis::X, RotationAxis::Y>;
	template<typename T>
	using EulerAnglesLocalZYX = EulerAngles<T, RotationAxis::Z, RotationAxis::Y, RotationAxis::X>;

	template<typename T>
	using EulerAnglesGlobalXYZ = EulerAnglesLocalZYX<T>;
	template<typename T>
	using EulerAnglesGlobalXZY = EulerAnglesLocalYZX<T>;
	template<typename T>
	using EulerAnglesGlobalYZX = EulerAnglesLocalXZY<T>;
	template<typename T>
	using EulerAnglesGlobalYXZ = EulerAnglesLocalZXY<T>;
	template<typename T>
	using EulerAnglesGlobalZXY = EulerAnglesLocalYXZ<T>;
	template<typename T>
	using EulerAnglesGlobalZYX = EulerAnglesLocalXYZ<T>;

	template<typename T>
	class StdEulerAngles : public EulerAnglesLocalZYX<T>
	{
	public:
		constexpr StdEulerAngles() = delete;
		constexpr StdEulerAngles(const EulerAnglesLocalZYX<T>& base) : EulerAnglesLocalZYX<T>{base} {}

		static constexpr StdEulerAngles FromMatrix(const Matrix3<T>& mat) {
			return EulerAnglesLocalZYX<T>::FromMatrix(mat);
		}
		static constexpr StdEulerAngles FromQuaternion(const Quaternion<T>& q) {
			return EulerAnglesLocalZYX<T>::FromQuaternion(q);
		}
		static constexpr StdEulerAngles FromVectorRadRPY(const Vector3<T>& vec) {
			return EulerAnglesLocalZYX<T>{vec(2), vec(1), vec(0)};
		}
		static constexpr StdEulerAngles FromVectorDegRPY(const Vector3<T>& vec) {
			return EulerAnglesLocalZYX<T>{DegreesToRadians(vec(2)), DegreesToRadians(vec(1)), DegreesToRadians(vec(0))};
		}

	public:
		constexpr T& Roll() { return (*this)(2); }
		constexpr T& Pitch() { return (*this)(1); }
		constexpr T& Yaw() { return (*this)(0); }

		constexpr const T& Roll() const { return (*this)(2); }
		constexpr const T& Pitch() const { return (*this)(1); }
		constexpr const T& Yaw() const { return (*this)(0); }

		constexpr Vector3<T> GetRadRPY() const {
			return Vector3<T>{Roll(), Pitch(), Yaw()};
		}
		constexpr Vector3<T> ToDegRPY() const {
			return Vector3<T>{RadiansToDegrees(Roll()), RadiansToDegrees(Pitch()), RadiansToDegrees(Yaw())};
		}
	};

	using StdEulerAnglesd = StdEulerAngles<double>;
	using StdEulerAnglesf = StdEulerAngles<float>;
}