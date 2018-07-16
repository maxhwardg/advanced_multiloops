//
// Created by max on 8/14/16.
// Contains types and functions for working with units of free energy.


#ifndef RNARK_ENERGY_HPP
#define RNARK_ENERGY_HPP

namespace librnary {

/// A value in Kcals/mol.
typedef double kcalmol_t;

/// The increment size of an energy value. Needed since energy is actually an integral value.
/// Currently desgiend to be the same as RNAstructure for compatibility.
const kcalmol_t EnergyIota = 0.1; // TODO: Change this when RNAstructure is no longer relied on.

/// Standard energy representation. Needs to be an integer type, and stores energy as a multiple of EnergyIota.
typedef int energy_t;

/**
 * @param kcalsmol Energy value is kcals/mol.
 * @return Energy value in increments.
 */
energy_t KCalToEnergy(kcalmol_t kcalsmol);

/**
 * Converts energy increments to kcal/mol.
 * @param e Energy in increments.
 * @return The value in kcal/mol.
 */
kcalmol_t EnergyToKCal(energy_t e);
}

#endif //RNARK_ENERGY_HPP
