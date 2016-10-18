# Algorithms for Advanced Multi-Loops
The repository contains a code library and associated programs for minimum free energy (MFE) Ribonucleic acid (RNA) folding algorithms, and sequence-structure free energy calculation programs. These programs correspond to the paper "Advanced Multi-Loop Algorithms for RNA Secondary Structure Prediction Reveal that the Simplest Model is Best" by Ward et al. The folding and free energy calculation programs in particular encompass the logarithmic model, which can be found under the Turner 1999 rules in the [Nearest Neighbor Database](http://rna.urmc.rochester.edu/NNDB/turner99/mb.html); the model of Aalberts & Nandagopal (see "A two-length-scale polymer theory for RNA loop free energies and helix stacking"); and the typical, linear model found in most MFE folding algorithms. In addition, there are a few undocumented goodies, not least among them is an algorithm for finding the MFE folding under an average multi-loop asymmetry model (beware, it is O(n^7)), an affine multi-loop asymmetry model folding algorithm, an algorithm for any non piecewise function taking both the number of branches and unpaired in a multi-loop, and a quite efficient brute force folding algorithm.

## Paper Results
The results from the paper can be found (in CSV format) in the "results" directory. The database of known structures is available upon request from the Mathews group. (See "Exact calculation of loop formation probability identifies folding motifs in RNA secondary structures" by Sloma & Mathews.)

## Compilation & Requirements
 To compile the programs you will need a few things. The first is a C++ compiler that supports the C++11 standard, and the second is CMake (version >= 3.5.x), which is the build system used. The code has only been tested using the default GCC (G++), and CMake installation under Ubuntu 16.04, though anything meeting the previous requirements should work. With this in mind, compilation can be done the usual CMake way. For example, let's assume we are using GNU Linux system, and we're in the root directory for this project:

> mkdir bin
> cd bin
> cmake ../ -DCMAKE_BUILD_TYPE=Release
> make all
> ./programs/logarithmic_fold  ../data_tables/

## Programs
### Input & Output Format
All programs take input and give output on the standard input and output streams. For RNA sequences, a string of characters from the alphabet {'A', 'U', 'G', 'C'} is expected. Any character that is non-whitepsace, and not from this alphabet is assumed to be an "unknown" nucleotide, and thus will not be paired in any folding algorithm. Secondary structures are always in dot-bracket format. No pseudoknots are allowed, so any characters other than '(' and ')' are assumed to be unpaired nucleotides.

### Folding Programs
All the folding programs expect a single argument on command line. This is the relative path to the data_tables directory (you will find this in the root directory of this project). This is a collection of files which contain the free energy parameters in RNAstructure format. We use this since RNAstructure is used to provide some basic free energy functions.

#### Logarithmic Fold Program
This program (logarithmic_fold) does MFE folding under the Turner 1999 logarithmic model of and parameters for multi-loops (see [Nearest Neighbor Database](http://rna.urmc.rochester.edu/NNDB/turner99/mb.html)). It expects a series of whitespace separated strings as input. Each string should be a primary sequence to be folded. For each primary sequence, the output will be 3 lines. The first will be the parsed version of the primary sequence (where 'X' denotes an unknown nucleotide), the second will be the MFE value, and the third will be a dot-bracket MFE structure. Here is an example usage:

> ./the_build_directory/logarithmic_fold data-tables/
> GAUCGAUCGACJJJJJJJXXXXauguagcuacgCCCCCC
> GAUCGAUCGACXXXXXXXXXXXXXXXXXXXXXXCCCCCC
> 0 kcal/mol
> .......................................
> GGGGAUAUGCUAGCUGAUGCUACGUAGCUGAcgugcuagcuagcugaucguGCUACGUAGCUACGGCGCGCCCCC
> GGGGAUAUGCUAGCUGAUGCUACGUAGCUGAXXXXXXXXXXXXXXXXXXXXGCUACGUAGCUACGGCGCGCCCCC
> -27.5 kcal/mol
> ((((...(((..((((..((((((((((.......................))))))))))..))))))))))).
> GGGGGGGGGGGGGGGGGGCCCCCCCCCCCCCCCCCCCCCC
> GGGGGGGGGGGGGGGGGGCCCCCCCCCCCCCCCCCCCCCC
> -50 kcal/mol
> ((((((((((((((((((...)))))))))))))))))).
> GGGUGAUGUGGGUACGUCGGGGGGAAAAAACCCCCCCCCCCCCCAAAAAAAAAGGGGGGGGGUUUUUUUUCCCCCCCCCCCCCC
> GGGUGAUGUGGGUACGUCGGGGGGAAAAAACCCCCCCCCCCCCCAAAAAAAAAGGGGGGGGGUUUUUUUUCCCCCCCCCCCCCC
> -48.2 kcal/mol
> (((.(((((....)))))((((((((((((..(((((((((............)))))))))..))))))))))))....))).

The time and space complexities are O(n^4) and O(n^3) respectively, so be careful about accidentally allocating all your RAM folding a 10k nt RNA.


#### Aalberts & Nandagopal Fold Program
This program (an_fold) does MFE folding under the Aalberts & Nandagopal multi-loop model (see [Nearest Neighbor Database](http://rna.urmc.rochester.edu/NNDB/turner99/mb.html)). It expects a series of whitespace separated strings as input. Each string should be a primary sequence to be folded. For each primary sequence, the output will be 3 lines. The first will be the parsed version of the primary sequence (where 'X' denotes an unknown nucleotide), the second will be the MFE value, and the third will be a dot-bracket MFE structure. Here is an example usage:

> ./the_build_directory/an_fold data-tables/
> GAUCGAUCGACJJJJJJJXXXXauguagcuacgCCCCCC
> GAUCGAUCGACXXXXXXXXXXXXXXXXXXXXXXCCCCCC
> 0 kcal/mol
> .......................................
> GGGGAUAUGCUAGCUGAUGCUACGUAGCUGAcgugcuagcuagcugaucguGCUACGUAGCUACGGCGCGCCCCC
> GGGGAUAUGCUAGCUGAUGCUACGUAGCUGAXXXXXXXXXXXXXXXXXXXXGCUACGUAGCUACGGCGCGCCCCC
> -27.5 kcal/mol
> ((((...(((..((((..((((((((((.......................))))))))))..))))))))))).
> GGGGGGGGGGGGGGGGGGCCCCCCCCCCCCCCCCCCCCCC
> GGGGGGGGGGGGGGGGGGCCCCCCCCCCCCCCCCCCCCCC
> -50 kcal/mol
> ((((((((((((((((((...)))))))))))))))))).
> GGGUGAUGUGGGUACGUCGGGGGGAAAAAACCCCCCCCCCCCCCAAAAAAAAAGGGGGGGGGUUUUUUUUCCCCCCCCCCCCCC
> GGGUGAUGUGGGUACGUCGGGGGGAAAAAACCCCCCCCCCCCCCAAAAAAAAAGGGGGGGGGUUUUUUUUCCCCCCCCCCCCCC
> -49.8 kcal/mol
> (((.(((((....)))))((((((((((((..(((((((((............)))))))))..))))))))))))..)))...

The time and space complexities are O(n^5) and O(n^4) respectively. This means that folding RNAs of more than about 200 nts can take minutes (or conceivably hours on a slow machine), and require gigabytes of RAM, so be careful.

#### Linear Fold Program
This program (linear_fold) does MFE folding under the typical linear model with the Turner 2004 parameters, which are those used by RNAstructure 5.8. It expects a series of whitespace separated strings as input. Each string should be a primary sequence to be folded. For each primary sequence, the output will be 3 lines. The first will be the parsed version of the primary sequence (where 'X' denotes an unknown nucleotide), the second will be the MFE value, and the third will be a dot-bracket MFE structure. Here is an example usage:

> ./the_build_directory/linear_fold data-tables/
> GAUCGAUCGACJJJJJJJXXXXauguagcuacgCCCCCC
> GAUCGAUCGACXXXXXXXXXXXXXXXXXXXXXXCCCCCC
> 0 kcal/mol
> .......................................
> GGGGAUAUGCUAGCUGAUGCUACGUAGCUGAcgugcuagcuagcugaucguGCUACGUAGCUACGGCGCGCCCCC
> GGGGAUAUGCUAGCUGAUGCUACGUAGCUGAXXXXXXXXXXXXXXXXXXXXGCUACGUAGCUACGGCGCGCCCCC
> -27.5 kcal/mol
> ((((...(((..((((..((((((((((.......................))))))))))..))))))))))).
> GGGGGGGGGGGGGGGGGGCCCCCCCCCCCCCCCCCCCCCC
> GGGGGGGGGGGGGGGGGGCCCCCCCCCCCCCCCCCCCCCC
> -50 kcal/mol
> ((((((((((((((((((...)))))))))))))))))).
> GGGUGAUGUGGGUACGUCGGGGGGAAAAAACCCCCCCCCCCCCCAAAAAAAAAGGGGGGGGGUUUUUUUUCCCCCCCCCCCCCC
> GGGUGAUGUGGGUACGUCGGGGGGAAAAAACCCCCCCCCCCCCCAAAAAAAAAGGGGGGGGGUUUUUUUUCCCCCCCCCCCCCC
> -48.4 kcal/mol
> (((.(((((....)))))((((((((((((..(((((((((............)))))))))..))))))))))))..)))...

This program can be a bit slower than expected, as it is done by emulation via the logarithmic folding algorithm.

### Energy Calculation Programs
These programs take an RNA sequence and a structure and calculate the expected free energy of that pair using a particular energy model. All these programs expect a sequence of whitespace separated sequence structure pairs as input. Each pair should be a sequence and structure separated by whitespace. These programs will all output 3 lines for each input pair. The parsed primary sequence, the parsed secondary structure, and the free energy. Every energy calculation program expects the first argument on command line to be the path to the data_tables folder. This should be the relative path to the data_tables directory (you will find this in the root directory of this project). In addition, the logarithmic and linear (but not Aalberts & Nandagoapl) programs can take a second command line argument. If it is "-describe", then they will also print out a breakdown of the free energy for every loop in the structure.

#### Logarithmic Energy Calculator
This program (logarithmic_efn) scores energy using the logarithmic model (see above). Here is an example usage:
> GGGUGAUGUGGGUACGUCGGGGGGAAAAAACCCCCCCCCCCCCCAAAAAAAAAGGGGGGGGGUUUUUUUUCCCCCCCCCCCCCC
> (((.(((((....)))))((((((((((((..(((((((((............)))))))))..))))))))))))....))).
> GGGUGAUGUGGGUACGUCGGGGGGAAAAAACCCCCCCCCCCCCCAAAAAAAAAGGGGGGGGGUUUUUUUUCCCCCCCCCCCCCC
> (((.(((((....)))))((((((((((((..(((((((((............)))))))))..))))))))))))....))).
> -48.2 kcal/mol

#### Aalberts & Nandagopal Energy Calculator
This program (an_efn) scores energy using the Aalberts & Nandagopal model (see above). Here is an example usage:
> GGGUGAUGUGGGUACGUCGGGGGGAAAAAACCCCCCCCCCCCCCAAAAAAAAAGGGGGGGGGUUUUUUUUCCCCCCCCCCCCCC
> (((.(((((....)))))((((((((((((..(((((((((............)))))))))..))))))))))))..)))...
> GGGUGAUGUGGGUACGUCGGGGGGAAAAAACCCCCCCCCCCCCCAAAAAAAAAGGGGGGGGGUUUUUUUUCCCCCCCCCCCCCC
> (((.(((((....)))))((((((((((((..(((((((((............)))))))))..))))))))))))..)))...
> -49.8 kcal/mol

#### Linear Energy Calculator
This program (linear_efn) scores energy using the linear model (see above). Here is an example usage:
> GGGUGAUGUGGGUACGUCGGGGGGAAAAAACCCCCCCCCCCCCCAAAAAAAAAGGGGGGGGGUUUUUUUUCCCCCCCCCCCCCC
> (((.(((((....)))))((((((((((((..(((((((((............)))))))))..))))))))))))..)))...
> GGGUGAUGUGGGUACGUCGGGGGGAAAAAACCCCCCCCCCCCCCAAAAAAAAAGGGGGGGGGUUUUUUUUCCCCCCCCCCCCCC
> (((.(((((....)))))((((((((((((..(((((((((............)))))))))..))))))))))))..)))...
> -48.4 kcal/mol

### A Note on data_tables
A common usage mistake is incorrectly providing the relative path to the data_tables directory to the executables. If you are getting strange looking structures with odd free energies, this is the most likely cause.
