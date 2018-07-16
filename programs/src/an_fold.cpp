//
// Created by max on 6/20/16.
//


#include <iostream>
#include <folders/aalberts_folder.hpp>

using namespace std;

/**
 * This program allows basic folding to be done under the Aalberts & Nandagopal model of multi-loops.
 * See "A two-length-scale polymer theory for RNA loop free energies and helix stacking" (2010).
 * The algorithm is O(n^5) and can be slow for moderately large RNA.
 * Also it can use a lot of memory, since the memory requirement is O(n^4).
 * It expects a single argument on command line. The relative path to the data_tables/ folder, which is an
 * RNAstructure format database of free energy parameters for RNA.
 * It takes whitespace separated RNA primary sequences as input on stdin.
 * These should be words from the alphabet {A,U,G,C} (note capitalization).
 * Any letters not from this alphabet are defined as invalid, and will be left unpaired.
 * For each sequence, three lines will be output. The parsed version of the primary sequence,
 * the minimum free energy of any for that primary sequence structure, and a minimum free energy structure.
 *
 * A sample usage:
 * ./the_build_directory/an_fold data-tables/
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
 * > -49.8 kcal/mol
 * > (((.(((((....)))))((((((((((((..(((((((((............)))))))))..))))))))))))..)))...

 */
int main(int argc, char **argv) {
    if (argc != 2) {
        cerr << "Expected exactly one argument, the relative path to the data tables folder." << endl;
        return 1;
    }

    string data_tables_path = argv[1];
    librnary::AalbertsModel model(data_tables_path);

    librnary::AalbertsFolder folder(model);
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
