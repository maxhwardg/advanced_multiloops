//
// Created by max on 5/4/16.
//

/**
 * Contains functions and types related to RNA primary structure.
 * The primary structure of an RNA is the nucleotide sequence.
 */

#ifndef RNARK_PRIMARY_STRUCTURE_HPP
#define RNARK_PRIMARY_STRUCTURE_HPP

#include <vector>
#include <string>
#include <random>

namespace librnary {
/// Represents a base on a nucleotide.
enum Base { A = 0, U, G, C, NUMBASES };
/// The primary structure of an RNA. The sequence of bases itself.
typedef std::vector<Base> PrimeStructure;
/// Converts a given character into the corresponding Base.
Base CharToBase(char c);
/// Converts a base into the corresponding character.
char BaseToChar(Base b);
/// Parses a string into a primary structure.
PrimeStructure StringToPrimary(const std::string &seq);
/// Represents a primary structure as a string.
std::string PrimaryToString(const PrimeStructure &primary);
/// Generates a random primary structure of particular length.
PrimeStructure RandomPrimary(std::default_random_engine &re, unsigned length);
}

#endif //RNARK_PRIMARY_STRUCTURE_HPP
