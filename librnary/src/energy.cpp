//
// Created by max on 8/14/16.
//

#include "energy.hpp"

#include <cmath>

using namespace std;

librnary::energy_t librnary::KCalToEnergy(librnary::kcalmol_t kcalsmol) {
	return static_cast<librnary::energy_t>(round(kcalsmol / librnary::EnergyIota));
}


librnary::kcalmol_t librnary::EnergyToKCal(librnary::energy_t e) {
	return e * librnary::EnergyIota;
}