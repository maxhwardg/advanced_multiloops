//
// Created by max on 6/20/16.
//


#include <iostream>
#include <models/nn_unpaired_model.hpp>
#include <scorers/nn_scorer.hpp>

using namespace std;

/**
 * This is an energy calculator program that assumes the logarithmic multi-loop model with NNDB parameters.
 * It expects the first argument to be the path to the data_tables/ directory.
 * The second argument is optional, and if it is "-describe" then a breakdown of free energy will printed.
 * The program expects a stream of input on stdin.
 * Every case is a pair of whitepsace separated strings. The first is expected to be the primary sequence.
 * The second should be a dot-bracket representation of the secondary structure to score.
 * The output comprises 3 lines, optionally followed by an energy breakdown.
 * The first line is the parsed version of the primary sequence. The second is the parsed secondary structure.
 * The third is the free energy.
 */
int main(int argc, char **argv) {
    if (argc < 2) {
        cerr << "Expected at least one argument, the relative path to the data tables folder." << endl;
        return 1;
    }

    bool describe = argc >= 3 && string(argv[2]) == "-describe";


    string data_tables_path = argv[1];
    librnary::NNUnpairedModel model(data_tables_path);

    librnary::NNScorer<librnary::NNUnpairedModel> scorer(model);

    while (cin.good()) {
        string primary, db;
        cin >> primary >> db;

        auto parsed_primary = librnary::StringToPrimary(primary);
        auto matching = librnary::DotBracketToMatching(db);

        scorer.SetRNA(parsed_primary);

        librnary::SSTree sstree(matching);

        cout << librnary::PrimaryToString(parsed_primary) << endl;
        cout << librnary::MatchingToDotBracket(matching) << endl;
        cout << librnary::EnergyToKCal(scorer.ScoreExterior(sstree.RootSurface())) << " kcal/mol" << endl;

        if (describe) {
            cout << scorer.TraceExterior(sstree.RootSurface()).Describe(' ', 0) << endl;
        }

        cin >> ws;
    }

}
