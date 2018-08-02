//
// Created by max on 8/2/17.
// Program to train parameters for the linear/affine multi-loop model.

#include <string>
#include <sstream>

#include "cxxopts.hpp"

#include "training/IBF_multiloop.hpp"
#include "models/nn_affine_model.hpp"
#include "folders/nn_affine_folder.hpp"
#include "scorers/nn_scorer.hpp"

using namespace std;

class LinearParameterSet {
protected:
	librnary::energy_t init, branch, unpaired;
public:
	LinearParameterSet(librnary::energy_t ini,
					   librnary::energy_t br,
					   librnary::energy_t up) : init(ini), branch(br), unpaired(up) {}
	string to_string() const {
		stringstream ss;
		ss << "init = " << init << " branch = " << branch << " unpaired = " << unpaired;
		return ss.str();
	}
	bool operator==(const LinearParameterSet &rhs) const {
		return init == rhs.init && branch == rhs.branch && unpaired == rhs.unpaired;
	}
	void LoadInto(librnary::NNAffineModel &model) const {
		model.SetMLParams(init, branch, unpaired);
	}
	class MultiInfo {
	protected:
		int branches{}, unpaired{};
	public:
		MultiInfo() = default;
		MultiInfo(const librnary::Surface &surf) {
			auto lr = librnary::LoopRegion(surf);
			branches = librnary::ExtractBranches(lr);
			unpaired = librnary::ExtractUnpaired(lr);
		}
		librnary::energy_t MLClosure(const librnary::NNAffineModel &model) const {
			return model.MLClosure(branches, unpaired);
		}
	};
};

int main(int argc, char **argv) {
	cxxopts::Options
		options("Train linear model using IBF", "Trains the parameters of the linear model using IBF.");

	options.add_options()
		("d,data_path", "Path to data_tables", cxxopts::value<string>()->default_value("data_tables/"))
		("c,ct_path", "Path to the folder of CTs", cxxopts::value<string>()->default_value("data_set/ct_files/"))
		("t,threads",
		 "Number of threads to use",
		 cxxopts::value<int>()->default_value(std::to_string(std::thread::hardware_concurrency())))
		("h,help", "Print help");

	string data_tables, ct_path;
	size_t threads;

	try {
		options.parse(argc, argv);
		data_tables = options["data_path"].as<string>();
		ct_path = options["ct_path"].as<string>();
		threads = static_cast<size_t>(options["threads"].as<int>());
		if (options.count("help") == 1) {
			cout << options.help({"", "Group"}) << endl;
			return 0;
		}

	} catch (const cxxopts::OptionException &e) {
		cout << "Argument parsing error: " << e.what() << endl;
		return 1;
	}

	// Read CTs.
	auto cts = librnary::ReadFilesInCTSetFormat(ct_path, cin);

	// Generate parameter list.
	vector<LinearParameterSet> params;

	for (librnary::energy_t init = 30; init <= 200; ++init) {
		for (librnary::energy_t br = -60; br <= 30; ++br) {
			for (librnary::energy_t up = -60; up <= 30; ++up) {
				params.emplace_back(init, br, up);
			}
		}
	}

	// Set up needed instances.
	librnary::NNAffineModel model(data_tables);
	librnary::NNAffineFolder folder(model);
	folder.SetMaxTwoLoop(30);
	folder.SetLonelyPairs(false);

	// Give multi-loops no cost.
	model.SetMLParams(0, 0, 0);

	// Make the trainer and train!
	librnary::IBFMultiLoop<LinearParameterSet,
						   librnary::NNAffineModel,
						   librnary::NNScorer<librnary::NNAffineModel>,
						   librnary::NNAffineFolder> trainer(model, cts, cout);
	trainer.SetNumStructureSeeds(5);
	trainer.SetThreads(threads);
	auto best_params = trainer.Train(params, params.front(), folder, 50);

	cout << "Best parameters: " << best_params.to_string() << endl;


	return 0;
}
