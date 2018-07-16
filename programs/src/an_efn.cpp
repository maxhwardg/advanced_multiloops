//
// Created by max on 6/20/16.
//


#include <iostream>
#include <scorers/aalberts_scorer.hpp>
#include <scorers/nn_scorer.hpp>

using namespace std;

/**
 * This is an energy calculator program that assumes the Aalberts & Nandagopal model of multi-loop energy.
 * It expects the first argument to be the path to the data_tables/ directory.
 * Currently an energy breakdown is not implemented for this model, though the free energy calculation works correctly.
 * The program expects a stream of input on stdin.
 * Every case is a pair of whitepsace separated strings. The first is expected to be the primary sequence.
 * The second should be a dot-bracket representation of the secondary structure to score.
 * The output comprises 3 lines.
 * The first line is the parsed version of the primary sequence. The second is the parsed secondary structure.
 * The third is the free energy.
 */
int main(int argc, char **argv) {

    if (argc != 2) {
        cerr << "Expected exactly one argument, the relative path to the data tables folder." << endl;
        return 1;
    }


    string data_tables_path = argv[1];
    librnary::AalbertsModel model(data_tables_path);

    librnary::AalbertsScorer scorer(model);

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


        cin >> ws;
    }

}
