//
// Created by max on 7/27/17.
// Program to train parameters for the logarithmic multi-loop model.

#include <string>
#include <sstream>

#include "cxxopts.hpp"

#include "training/IBF_multiloop.hpp"
#include "models/nn_unpaired_model.hpp"
#include "folders/nn_unpaired_folder.hpp"
#include "scorers/nn_scorer.hpp"

using namespace std;

class LogarithmicParameterSet {
protected:
	librnary::energy_t init, branch, unpaired;
	librnary::kcalmol_t log_mult;
	int pivot;
public:
	LogarithmicParameterSet(librnary::energy_t ini,
							librnary::energy_t br,
							librnary::energy_t up,
							librnary::kcalmol_t log,
							int piv) : init(ini), branch(br), unpaired(up), log_mult(log), pivot(piv) {}
	string to_string() const {
		stringstream ss;
		ss << "init = " << init << " branch = " << branch << " unpaired = "
		   << unpaired << " log_mult = " << log_mult << " pivot = " << pivot;
		return ss.str();
	}
	bool operator==(const LogarithmicParameterSet &rhs) const {
		return init == rhs.init && branch == rhs.branch
			&& unpaired == rhs.unpaired && log_mult == rhs.log_mult && pivot == rhs.pivot;
	}
	void LoadInto(librnary::NNUnpairedModel &model) const {
		model.SetMLParams(init, branch, unpaired, log_mult, pivot);
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
		librnary::energy_t MLClosure(const librnary::NNUnpairedModel &model) const {
			return model.MLClosure(unpaired, branches);
		}
	};
};

int main(int argc, char **argv) {
	cxxopts::Options
		options("Train logarithmic model using IBF", "Trains the parameters of the logarithmic model using IBF.");

	options.add_options()
		("d,data_path", "Path to data_tables", cxxopts::value<string>()->default_value("data_tables/"))
		("c,ct_path", "Path to the folder of CTs", cxxopts::value<string>()->default_value("RNAs/ArchiveIII/"))
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
	vector<LogarithmicParameterSet> params;


	for (librnary::energy_t init = 30; init <= 200; ++init) {
		for (librnary::energy_t br = -30; br <= 30; ++br) {
			for (librnary::energy_t up = -30; up <= 30; ++up) {
				for (librnary::energy_t log = -30; log <= 30; ++log) {
					for (int pivot = 1; pivot <= 8; pivot += 1) {
						params.emplace_back(init, br, up, librnary::EnergyToKCal(log), pivot);
					}
				}
			}
		}
	}

	// Set up needed instances.
	librnary::NNUnpairedModel model(data_tables);
	librnary::NNUnpairedFolder folder(model);
	folder.SetMaxTwoLoop(30);
	folder.SetLonelyPairs(false);

	// Give multi-loops no cost.
	model.SetMLParams(0, 0, 0, 0, 999999);

	// Make the trainer and train!
	librnary::IBFMultiLoop<LogarithmicParameterSet,
						   librnary::NNUnpairedModel,
						   librnary::NNScorer<librnary::NNUnpairedModel>,
						   librnary::NNUnpairedFolder> trainer(model, cts, clog);
	trainer.SetNumStructureSeeds(5);
	trainer.SetThreads(threads);
	auto best_params = trainer.Train(params, params.front(), folder, 50);

	cout << "Best parameters: " << best_params.to_string() << endl;

	return 0;
}