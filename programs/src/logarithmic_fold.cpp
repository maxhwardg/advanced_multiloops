//
// Created by max on 6/20/16.
//


#include <iostream>
#include <folders/nn_unpaired_folder.hpp>

using namespace std;

/**
 * This program allows basic folding to be done under the logarithmic model using the default valued published on the
 * RNA Nearest Neighbor Database (NNDB).
 * The algorithm is O(n^4) and can be slow for very large RNA.
 * It expects a single argument on command line. The relative path to the data_tables/ folder, which is an
 * RNAstructure format database of free energy parameters for RNA.
 * It takes whitespace separated RNA primary sequences as input on stdin.
 * These should be words from the alphabet {A,U,G,C} (note capitalization).
 * Any letters not from this alphabet are defined as invalid, and will be left unpaired.
 * For each sequence, three lines will be output. The parsed version of the primary sequence,
 * the minimum free energy of any for that primary sequence structure, and a minimum free energy structure.
 *
 * A sample usage:
 * ./the_build_directory/logarithmic_fold data-tables/
 * GAUCGAUCGACJJJJJJJXXXXauguagcuacgCCCCCC
 * > GAUCGAUCGACXXXXXXXXXXXXXXXXXXXXXXCCCCCC
 * > 0 kcal/mol
 * > .......................................
 * GGGGAUAUGCUAGCUGAUGCUACGUAGCUGAcgugcuagcuagcugaucguGCUACGUAGCUACGGCGCGCCCCC
 * > GGGGAUAUGCUAGCUGAUGCUACGUAGCUGAXXXXXXXXXXXXXXXXXXXXGCUACGUAGCUACGGCGCGCCCCC
 * > -27.5 kcal/mol
 * > ((((...(((..((((..((((((((((.......................))))))))))..))))))))))).
 * GGGGGGGGGGGGGGGGGGCCCCCCCCCCCCCCCCCCCCCC
 * > GGGGGGGGGGGGGGGGGGCCCCCCCCCCCCCCCCCCCCCC
 * > -50 kcal/mol
 * > ((((((((((((((((((...)))))))))))))))))).
 * GGGUGAUGUGGGUACGUCGGGGGGAAAAAACCCCCCCCCCCCCCAAAAAAAAAGGGGGGGGGUUUUUUUUCCCCCCCCCCCCCC
 * > GGGUGAUGUGGGUACGUCGGGGGGAAAAAACCCCCCCCCCCCCCAAAAAAAAAGGGGGGGGGUUUUUUUUCCCCCCCCCCCCCC
 * > -48.2 kcal/mol
 * > (((.(((((....)))))((((((((((((..(((((((((............)))))))))..))))))))))))....))).

 */
int main(int argc, char **argv) {
    if (argc != 2) {
        cerr << "Expected exactly one argument, the relative path to the data tables folder." << endl;
        return 1;
    }

    string data_tables_path = argv[1];
    librnary::NNUnpairedModel model(data_tables_path);

    librnary::NNUnpairedFolder folder(model);
    folder.SetLonelyPairs(false);

    while (cin.good()) {
        string primary;
        cin >> primary;

        auto parsed_primary = librnary::StringToPrimary(primary);

        librnary::energy_t mfe = folder.Fold(parsed_primary);
        cout << librnary::PrimaryToString(parsed_primary) << endl;
        cout << librnary::EnergyToKCal(mfe) << " kcal/mol" << endl;
        cout << librnary::MatchingToDotBracket(folder.Traceback()) << endl;

        cin >> ws;
    }

}
