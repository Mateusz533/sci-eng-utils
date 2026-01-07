#pragma once
//
#include "UnitsSI.hpp"

namespace Physics
{
	namespace Constants
	{
		using namespace Units::SI;

		inline constexpr Radians<> RADIAN_PER_REVOLUTION{2 * M_PI};
		inline constexpr MetersPerSecond<> SPEED_OF_LIGHT{2.99792458e8};
		inline constexpr JouleSeconds<> PLANCK_CONSTANT{6.62607015e-34};
		inline constexpr Coulombs<> ELEMENTARY_CHARGE{1.602176634e-19};
		inline constexpr JoulesPerKelvin<> BOLTZMANN_CONSTANT{1.380649e-23};
		inline constexpr PartsPerMole<> AVOGADRO_CONSTANT{6.02214076e23};

		inline constexpr JouleSecondsPerRadian<> DIRAC_CONSTANT{PLANCK_CONSTANT / RADIAN_PER_REVOLUTION};
		inline constexpr HenriesPerMeter<> MAGNETIC_CONSTANT{4e-7 * M_PI};
		inline constexpr FaradsPerMeter<> ELECTRIC_CONSTANT = Scale<>{1} / (MAGNETIC_CONSTANT * SPEED_OF_LIGHT.Power<2>());
	}

	namespace Calculate
	{
		using namespace Units::SI;

		/* Mechanics */;

		template<typename T>
		constexpr MetersPerSecond<T> TangentialVelocity(RadiansPerSecond<T> angularVelocity, RadialMeters<T> radius) {
			return angularVelocity * radius;
		}
		template<typename T>
		constexpr MetersPerSecondSquare<T> TangentialAcceleration(RadiansPerSecondSquare<T> angularAcceleration, RadialMeters<T> radius) {
			return angularAcceleration * radius;
		}
		template<typename T>
		constexpr NewtonMeterSeconds<T> AngularMomentum(NewtonSeconds<T> momentum, RadialMeters<T> radius) {
			return momentum * radius;
		}
		template<typename T>
		constexpr Joules<T> Energy(KiloGrams<T> mass, MetersPerSecond<T> velocity) {
			return mass * velocity * velocity / Scale<T>(2);
		}
		template<typename T>
		constexpr Joules<T> Energy(KiloGramMetersSquare<T> momentOfInertia, RadiansPerSecond<T> angularVelocity) {
			return momentOfInertia * angularVelocity * angularVelocity / Scale<T>(2);
		}
		template<typename T>
		constexpr Watts<T> Power(Newtons<T> force, MetersPerSecond<T> velocity) {
			return force * velocity;
		}
		template<typename T>
		constexpr Watts<T> Power(NewtonMeters<T> torque, RadiansPerSecond<T> angularVelocity) {
			return torque * angularVelocity;
		}
		template<typename T>
		constexpr KiloGramMetersSquare<T> Inertia(KiloGrams<T> mass, RadialMeters<T> radius, Scale<T> factor) {
			return factor * mass * radius * radius;
		}

		/* Electronics */;

		template<typename T>
		constexpr Amperes<T> Current(Coulombs<T> charge, Seconds<T> time) {
			return charge / time;
		}
		template<typename T>
		constexpr Coulombs<T> Charge(Amperes<T> current, Seconds<T> time) {
			return current * time;
		}
		template<typename T>
		constexpr Ohms<T> Resistance(Volts<T> voltage, Amperes<T> current) {
			return voltage / current;
		}
		template<typename T>
		constexpr Coulombs<T> Charge(Farads<T> capacity, Volts<T> voltage) {
			return capacity * voltage;
		}
		template<typename T>
		constexpr Webers<T> MagneticFlux(Henries<T> inductance, Amperes<T> current) {
			return inductance * current;
		}

		/* Oscilations and waves */;

		template<typename T>
		constexpr RadiansPerMeter<T> WaveNumber(Meters<T> waveLength) {
			return Constants::RADIAN_PER_REVOLUTION / waveLength;
		}
		template<typename T>
		constexpr Seconds<T> Period(Hertzes<T> frequency) {
			return Scale<T>(1) / frequency;
		}

		/* Relativistic mechanics */;

		template<typename T>
		constexpr Scale<T> TimeDilation(MetersPerSecond<T> velocity) {
			return Scale<T>(1) / (Scale<T>(1) - (velocity / Constants::SPEED_OF_LIGHT).template Power<2>()).Sqrt();
		}

		// TODO: Add more patterns
	}
}