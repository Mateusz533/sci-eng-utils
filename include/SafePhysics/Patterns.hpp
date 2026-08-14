#pragma once
//
#include <numbers>
//
#include "UnitsSI.hpp"

namespace Physics
{
	namespace Constants
	{
		using namespace Units::SI;

		inline constexpr Radians<> RADIAN_PER_REVOLUTION{2.0 * std::numbers::pi};
		inline constexpr MetersPerSecond<> SPEED_OF_LIGHT{2.99792458e8};
		inline constexpr JouleSeconds<> PLANCK_CONSTANT{6.62607015e-34};
		inline constexpr Coulombs<> ELEMENTARY_CHARGE{1.602176634e-19};
		inline constexpr JoulesPerKelvin<> BOLTZMANN_CONSTANT{1.380649e-23};
		inline constexpr PartsPerMole<> AVOGADRO_CONSTANT{6.02214076e23};

		inline constexpr JouleSecondsPerRadian<> DIRAC_CONSTANT{PLANCK_CONSTANT / RADIAN_PER_REVOLUTION};
		inline constexpr HenriesPerMeter<> MAGNETIC_CONSTANT{4e-7 * std::numbers::pi};
		inline constexpr FaradsPerMeter<> ELECTRIC_CONSTANT = Scale<>{1} / (MAGNETIC_CONSTANT * SPEED_OF_LIGHT.Power<2>());
	}

	namespace Calculate
	{
		using namespace Units::SI;

		/* Mechanics */;

		template<Arithmetic T>
		constexpr MetersPerSecond<T> TangentialVelocity(RadiansPerSecond<T> angularVelocity, RadialMeters<T> radius) noexcept {
			return angularVelocity * radius;
		}
		template<Arithmetic T>
		constexpr MetersPerSecondSquare<T> TangentialAcceleration(RadiansPerSecondSquare<T> angularAcceleration, RadialMeters<T> radius) noexcept {
			return angularAcceleration * radius;
		}
		template<Arithmetic T>
		constexpr NewtonMeterSeconds<T> AngularMomentum(NewtonSeconds<T> momentum, RadialMeters<T> radius) noexcept {
			return momentum * radius;
		}
		template<Arithmetic T>
		constexpr Joules<T> Energy(KiloGrams<T> mass, MetersPerSecond<T> velocity) noexcept {
			return mass * velocity * velocity / Scale<T>(2);
		}
		template<Arithmetic T>
		constexpr Joules<T> Energy(KiloGramMetersSquare<T> momentOfInertia, RadiansPerSecond<T> angularVelocity) noexcept {
			return momentOfInertia * angularVelocity * angularVelocity / Scale<T>(2);
		}
		template<Arithmetic T>
		constexpr Watts<T> Power(Newtons<T> force, MetersPerSecond<T> velocity) noexcept {
			return force * velocity;
		}
		template<Arithmetic T>
		constexpr Watts<T> Power(NewtonMeters<T> torque, RadiansPerSecond<T> angularVelocity) noexcept {
			return torque * angularVelocity;
		}
		template<Arithmetic T>
		constexpr KiloGramMetersSquare<T> Inertia(KiloGrams<T> mass, RadialMeters<T> radius, Scale<T> factor) noexcept {
			return factor * mass * radius * radius;
		}

		/* Electronics */;

		template<Arithmetic T>
		constexpr Amperes<T> Current(Coulombs<T> charge, Seconds<T> time) noexcept {
			return charge / time;
		}
		template<Arithmetic T>
		constexpr Coulombs<T> Charge(Amperes<T> current, Seconds<T> time) noexcept {
			return current * time;
		}
		template<Arithmetic T>
		constexpr Ohms<T> Resistance(Volts<T> voltage, Amperes<T> current) noexcept {
			return voltage / current;
		}
		template<Arithmetic T>
		constexpr Coulombs<T> Charge(Farads<T> capacity, Volts<T> voltage) noexcept {
			return capacity * voltage;
		}
		template<Arithmetic T>
		constexpr Webers<T> MagneticFlux(Henries<T> inductance, Amperes<T> current) noexcept {
			return inductance * current;
		}

		/* Oscillations and waves */;

		template<Arithmetic T>
		constexpr RadiansPerMeter<T> WaveNumber(Meters<T> waveLength) noexcept {
			return Constants::RADIAN_PER_REVOLUTION / waveLength;
		}
		template<Arithmetic T>
		constexpr Seconds<T> Period(Hertzes<T> frequency) noexcept {
			return Scale<T>(1) / frequency;
		}

		/* Relativistic mechanics */;

		template<Arithmetic T>
		constexpr Scale<T> TimeDilation(MetersPerSecond<T> velocity) noexcept {
			return Scale<T>(1) / (Scale<T>(1) - (velocity / Constants::SPEED_OF_LIGHT).template Power<2>()).Sqrt();
		}

		// TODO: Add more patterns
	}
}