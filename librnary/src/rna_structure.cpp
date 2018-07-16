//
// Created by max on 5/30/16.
//

#include "rna_structure.hpp"

#include <cstring>

using namespace std;

unique_ptr<datatable> librnary::LoadDatatable(const string &path) {
	unique_ptr<datatable> dt(new datatable());
	vector<string> names = {"loop", "stack", "tstackh", "tstacki",
							"tloop", "miscloop", "dangle", "int22", "int21", "coaxial",
							"tstackcoax", "coaxstack", "tstack", "tstackm", "triloop",
							"int11", "hexaloop", "tstacki23", "tstacki1n"};

	char *paths[19];

	for (unsigned i = 0; i < names.size(); ++i) {
		string file_path = path + "/" + names[i] + ".dat";
		paths[i] = new char[file_path.size() + 1];
		strcpy(paths[i], file_path.c_str());
	}

	int rv = opendat(paths[0], paths[1], paths[2], paths[3], paths[4], paths[5],
					 paths[6], paths[7], paths[8], paths[9], paths[10], paths[11],
					 paths[12], paths[13], paths[14], paths[15], paths[16],
					 paths[17], paths[18], dt.get());
	if (rv == 0) {
		cerr << "An error occured when loading the data tables" << endl;
	}
	for (char* p : paths) {
		delete[] p;
	}
	return dt;
}

unique_ptr<structure> librnary::LoadStructure(const librnary::PrimeStructure &primary) {
	unique_ptr<structure> seq(new structure());
	seq->allocate(static_cast<int>(primary.size()));
	for (int i = 1; i <= seq->GetSequenceLength(); ++i) {
		librnary::Base nuc = primary[i - 1];
		if (nuc == librnary::A) seq->numseq[i] = 1;
		else if (nuc == librnary::C) seq->numseq[i] = 2;
		else if (nuc == librnary::G) seq->numseq[i] = 3;
		else if (nuc == librnary::U) seq->numseq[i] = 4;
		else seq->numseq[i] = 0;
		seq->nucs[i] = librnary::BaseToChar(nuc);
		seq->hnumber[i] = static_cast<short>(i);
	}
	return seq;
}

std::unique_ptr<structure> librnary::LoadStructure(const librnary::PrimeStructure &primary_seq,
												   const librnary::Matching &match) {
	assert(primary_seq.size() == match.size()); // Ensure there are no indexing issues later.
	unique_ptr<structure> struc_ptr = LoadStructure(primary_seq);
	struc_ptr->AddStructure(); // Create a single structure for the structure (nomeclature is weird).
	for (int i = 0; i < static_cast<int>(match.size()); ++i) {
		// Make sure the matching contains only valid base pairs.
		assert(match[i] == i || librnary::ValidPair(primary_seq[i], primary_seq[match[i]]));
		if (i < match[i]) // Don't add a pair twice.
			struc_ptr->SetPair(i + 1, match[i] + 1); // structure uses 1 indexing, so +1 everywhere.
	}
	return struc_ptr;
}

unique_ptr<structure> librnary::LoadStructure(const string &ct_file) {
	unique_ptr<structure> seq(new structure());
	long error_code = seq->openct(ct_file.c_str());
	if (error_code != 0) {
		cerr << "Could not read " << ct_file << endl;
		throw seq->openct(ct_file.c_str());
	}
	return seq;
}

unique_ptr<structure> librnary::LoadSeqFile(const std::string &seq_file) {
	unique_ptr<structure> seq(new structure());
	long error_code = seq->openseq(seq_file.c_str());
	if (error_code != 1)
		throw seq->openseq(seq_file.c_str());
	return seq;
}

librnary::Matching librnary::StructureToMatching(structure &struc) {
	auto match = EmptyMatching(static_cast<unsigned>(struc.GetSequenceLength()));
	for (int i = 0; i < static_cast<int>(match.size()); ++i) {
		int p = struc.GetPair(i + 1);
		match[i] = p == 0 ? i : p - 1;
	}
	return match;
}

librnary::PrimeStructure librnary::StructureToPrimary(structure &struc) {
	string prim;
	for (int i = 0; i < struc.GetSequenceLength(); ++i)
		prim.push_back(struc.nucs[i + 1]);
	return librnary::StringToPrimary(prim);
}

librnary::energy_t librnary::RunEFN2(datatable &dt, structure &struc) {
	efn2(&dt, &struc, 1);
	return struc.GetEnergy(1);
}

librnary::energy_t librnary::RunEFN2WithSimpleMulti(datatable &dt, structure &struc) {
	efn2(&dt, &struc, 1, true);
	return struc.GetEnergy(1);
}

librnary::energy_t librnary::RunMFEFold(datatable &dt, structure &struc, int two_loop_max_nts) {
	dynamic(&struc, &dt, 1, 0, 0, nullptr, false, nullptr, two_loop_max_nts, false, false);
	return struc.GetEnergy(1);
}
