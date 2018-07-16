//
// Created by max on 6/20/16.
//

#ifndef RNARK_READ_CTS_HPP
#define RNARK_READ_CTS_HPP

#include <string>
#include <istream>
#include <vector>

#include "primary_structure.hpp"
#include "secondary_structure.hpp"

namespace librnary {

struct CTData {
	std::string name;
	PrimeStructure primary;
	Matching match;
};

/// Reads a single file give a path to the CT file.
CTData ReadFile(const std::string &file, const std::string &path);

/**
 * Takes the base path to a ct folder, and an input stream.
 * Reads data from in, assuming CT Set format, and returns every data set as a vector.
 * CT Set format has the following format:
 * ct_file_name 1
 * ct_file_name 2
 * ...
 * end set
 * ct_file_name 1
 * ...
 * end
 */
std::vector<std::vector<CTData>> ReadFilesInCTSetFormat(const std::string &base_path, std::istream &in);

/**
 * Reads all the CT files given in a stream. Assumes CT Set format (see ReadFilesInCTSetFormat).
 * @param base_path The base prefix to any ct file path.
 * @param in An input stream. Will read every line from this stream and assumes they constitute a CT Set format file.
 * @return A list of all the CTs.
 */
std::vector<librnary::CTData> ReadAllCTs(const std::string &base_path, std::istream &in);
}

#endif //RNARK_READ_CTS_HPP
