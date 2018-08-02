//
// Created by max on 8/29/17.
// Program to train parameters for the Aalberts + Nandagopal multi-loop model.

#include <string>
#include <sstream>
#include <folders/aalberts_folder.hpp>

#include "cxxopts.hpp"

#include "training/IBF_multiloop_aalberts.hpp"

using namespace std;

class AalbertsParameterSet {
protected:
	double a, b;
	librnary::kcalmol_t c;
public:
	AalbertsParameterSet(double _a, double _b, librnary::kcalmol_t _c) : a(_a), b(_b), c(_c) {}
	string to_string() const {
		stringstream ss;
		ss << "a = " << a << " b = " << b << " C = " << c;
		return ss.str();
	}
	bool operator==(const AalbertsParameterSet &rhs) const {
		return a == rhs.a && b == rhs.b && c == rhs.c;
	}
	void LoadInto(librnary::AalbertsModel &model) const {
		model.SetAdditiveConstant(c);
		model.SetNCoeffBase(a);
		model.SetMCoeffBase(b);
	}
	class MultiInfo {
	protected:
		int Alength{}, Blength{};
	public:
		MultiInfo() = default;
		/**
		 * This exists to avoid a compiler issue with templated inheritance.
		 */
		MultiInfo(const librnary::Surface &surf) {
			assert(false);
		}
		MultiInfo(const librnary::AalbertsScorer &scorer, const librnary::Surface &surf) {
			auto trace = scorer.TraceMLConfig(surf);
			Alength = get<1>(trace)["N"];
			Blength = get<1>(trace)["M"];
		}
		librnary::energy_t MLClosure(const librnary::AalbertsModel &model) const {
			return model.MLInit(Alength, Blength);
		}
	};
};

int main(int argc, char **argv) {
	cxxopts::Options
		options("Traing AN model using IBF", "Trains the parameters of the AN model using IBF.");

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
	vector<AalbertsParameterSet> params;

	librnary::AalbertsModel model(data_tables);

	for (int a = -100; a <= 50; ++a) {
		for (int b = -50; b <= 50; ++b) {
			for (librnary::energy_t C = -50; C <= 50; ++C) {
				params.emplace_back(model.NCoeffBase() + a * 0.05, model.MCoeffBase() + b * 0.05, C*0.1);
			}
		}
	}

	// Set up needed instances.
	librnary::AalbertsFolder folder(model);
	folder.SetMaxTwoLoop(30);
	folder.SetLonelyPairs(false);

	// Give multi-loops no cost.
	model.SetNCoeffBase(0);
	model.SetMCoeffBase(0);
	model.SetAdditiveConstant(0);

	// Make the trainer and train!
	librnary::IBFMultiLoopAalberts<AalbertsParameterSet> trainer(model, cts, clog);
	trainer.SetNumStructureSeeds(5);
	trainer.SetThreads(threads);
	auto best_params = trainer.Train(params, params.back(), folder, 50);

	cout << "Best parameters: " << best_params.to_string() << endl;


	return 0;
}

