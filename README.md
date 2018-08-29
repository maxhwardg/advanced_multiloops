# Algorithms for Advanced Multi-Loops
The repository contains a code library and associated programs for two things. First, minimum free energy (MFE) Ribonucleic acid (RNA) folding algorithms, and sequence-structure free energy calculation programs. Second, algorithms for training model parameters, particularly multi-loop models.

The MFE folding and energy calculation algorithms correspond to the paper "Advanced Multi-Loop Algorithms for RNA Secondary Structure Prediction Reveal that the Simplest Model is Best" by Ward et al. The folding and free energy calculation programs in particular encompass the logarithmic model, which can be found under the Turner 1999 rules in the [Nearest Neighbor Database](http://rna.urmc.rochester.edu/NNDB/turner99/mb.html); the model of Aalberts & Nandagopal (see "A two-length-scale polymer theory for RNA loop free energies and helix stacking"); and the typical, linear model found in most MFE folding algorithms. In addition, there are a few undocumented goodies, not least among them is an algorithm for finding the MFE folding under an average multi-loop asymmetry model (beware, it is O(n^7)), an affine multi-loop asymmetry model folding algorithm, an algorithm for any non piecewise function taking both the number of branches and unpaired in a multi-loop, and a quite efficient brute force folding algorithm.

The parameter training algorithms correspond to a paper called "Determining Parameters for Non-Linear Models of Multi-Loop Free Energy Change". Programs are included that can train parameters for the linear, logarithmic, and Aalberts & Nandagopal models. Further, the complete training set used in the paper can be found.

Instructions for using these programs are now included.

## Compilation & Requirements
 To compile the programs you will need a few things. The first is a C++ compiler that supports the C++11 standard, and the second is CMake (version >= 3.5.x), which is the build system used. The code has only been tested using the default GCC (G++), and CMake installation under Ubuntu 16.04, though anything meeting the previous requirements should work. With this in mind, compilation can be done the usual CMake way. For example, let's assume we are using GNU Linux system, and we're in the root directory for this project:

```
mkdir bin
cd bin
cmake ../ -DCMAKE_BUILD_TYPE=Release
make all
./programs/logarithmic_fold  ../data_tables/
```

## Folding & Energy Calculation Programs
### Input & Output Format
All programs take input and give output on the standard input and output streams. For RNA sequences, a string of characters from the alphabet {'A', 'U', 'G', 'C'} is expected. Any character that is non-whitepsace, and not from this alphabet is assumed to be an "unknown" nucleotide, and thus will not be paired in any folding algorithm. Secondary structures are always in dot-bracket format. No pseudoknots are allowed, so any characters other than '(' and ')' are assumed to be unpaired nucleotides.

### Folding Programs
All the folding programs expect a single argument on command line. This is the relative path to the data_tables directory (you will find this in the root directory of this project). This is a collection of files which contain the free energy parameters in RNAstructure format. We use this since RNAstructure is used to provide some basic free energy functions.

#### Logarithmic Fold Program
This program (logarithmic_fold) does MFE folding under the Turner 1999 logarithmic model and parameters for multi-loops (see [Nearest Neighbor Database](http://rna.urmc.rochester.edu/NNDB/turner99/mb.html)). It expects a series of whitespace separated strings as input. Each string should be a primary sequence to be folded. For each primary sequence, the output will be 3 lines. The first will be the parsed version of the primary sequence (where 'X' denotes an unknown nucleotide), the second will be the MFE value, and the third will be a dot-bracket MFE structure. Here is an example usage:

```
./the_build_directory/logarithmic_fold data-tables/
GAUCGAUCGACJJJJJJJXXXXauguagcuacgCCCCCC
> GAUCGAUCGACXXXXXXXXXXXXXXXXXXXXXXCCCCCC
> 0 kcal/mol
> .......................................
GGGGAUAUGCUAGCUGAUGCUACGUAGCUGAcgugcuagcuagcugaucguGCUACGUAGCUACGGCGCGCCCCC
> GGGGAUAUGCUAGCUGAUGCUACGUAGCUGAXXXXXXXXXXXXXXXXXXXXGCUACGUAGCUACGGCGCGCCCCC
> -27.5 kcal/mol
> ((((...(((..((((..((((((((((.......................))))))))))..))))))))))).
GGGGGGGGGGGGGGGGGGCCCCCCCCCCCCCCCCCCCCCC
> GGGGGGGGGGGGGGGGGGCCCCCCCCCCCCCCCCCCCCCC
> -50 kcal/mol
> ((((((((((((((((((...)))))))))))))))))).
GGGUGAUGUGGGUACGUCGGGGGGAAAAAACCCCCCCCCCCCCCAAAAAAAAAGGGGGGGGGUUUUUUUUCCCCCCCCCCCCCC
> GGGUGAUGUGGGUACGUCGGGGGGAAAAAACCCCCCCCCCCCCCAAAAAAAAAGGGGGGGGGUUUUUUUUCCCCCCCCCCCCCC
> -48.2 kcal/mol
> (((.(((((....)))))((((((((((((..(((((((((............)))))))))..))))))))))))....))).
```

The time and space complexities are O(n^4) and O(n^3) respectively, so be careful about accidentally allocating all your RAM folding a 10k nt RNA.


#### Aalberts & Nandagopal Fold Program
This program (an_fold) does MFE folding under the Aalberts & Nandagopal multi-loop model (see [Nearest Neighbor Database](http://rna.urmc.rochester.edu/NNDB/turner99/mb.html)). It expects a series of whitespace separated strings as input. Each string should be a primary sequence to be folded. For each primary sequence, the output will be 3 lines. The first will be the parsed version of the primary sequence (where 'X' denotes an unknown nucleotide), the second will be the MFE value, and the third will be a dot-bracket MFE structure. Here is an example usage:

```
./the_build_directory/an_fold data-tables/
GAUCGAUCGACJJJJJJJXXXXauguagcuacgCCCCCC
> GAUCGAUCGACXXXXXXXXXXXXXXXXXXXXXXCCCCCC
> 0 kcal/mol
> .......................................
GGGGAUAUGCUAGCUGAUGCUACGUAGCUGAcgugcuagcuagcugaucguGCUACGUAGCUACGGCGCGCCCCC
> GGGGAUAUGCUAGCUGAUGCUACGUAGCUGAXXXXXXXXXXXXXXXXXXXXGCUACGUAGCUACGGCGCGCCCCC
> -27.5 kcal/mol
> ((((...(((..((((..((((((((((.......................))))))))))..))))))))))).
GGGGGGGGGGGGGGGGGGCCCCCCCCCCCCCCCCCCCCCC
> GGGGGGGGGGGGGGGGGGCCCCCCCCCCCCCCCCCCCCCC
> -50 kcal/mol
> ((((((((((((((((((...)))))))))))))))))).
GGGUGAUGUGGGUACGUCGGGGGGAAAAAACCCCCCCCCCCCCCAAAAAAAAAGGGGGGGGGUUUUUUUUCCCCCCCCCCCCCC
> GGGUGAUGUGGGUACGUCGGGGGGAAAAAACCCCCCCCCCCCCCAAAAAAAAAGGGGGGGGGUUUUUUUUCCCCCCCCCCCCCC
> -49.8 kcal/mol
> (((.(((((....)))))((((((((((((..(((((((((............)))))))))..))))))))))))..)))...
```

The time and space complexities are O(n^5) and O(n^4) respectively. This means that folding RNAs of more than about 200 nts can take minutes (or conceivably hours on a slow machine), and require gigabytes of RAM, so be careful.

#### Linear Fold Program
This program (linear_fold) does MFE folding under the typical linear model with the Turner 2004 parameters, which are those used by RNAstructure 5.8. It expects a series of whitespace separated strings as input. Each string should be a primary sequence to be folded. For each primary sequence, the output will be 3 lines. The first will be the parsed version of the primary sequence (where 'X' denotes an unknown nucleotide), the second will be the MFE value, and the third will be a dot-bracket MFE structure. Here is an example usage:

```
./the_build_directory/linear_fold data-tables/
GAUCGAUCGACJJJJJJJXXXXauguagcuacgCCCCCC
> GAUCGAUCGACXXXXXXXXXXXXXXXXXXXXXXCCCCCC
> 0 kcal/mol
> .......................................
GGGGAUAUGCUAGCUGAUGCUACGUAGCUGAcgugcuagcuagcugaucguGCUACGUAGCUACGGCGCGCCCCC
> GGGGAUAUGCUAGCUGAUGCUACGUAGCUGAXXXXXXXXXXXXXXXXXXXXGCUACGUAGCUACGGCGCGCCCCC
> -27.5 kcal/mol
> ((((...(((..((((..((((((((((.......................))))))))))..))))))))))).
GGGGGGGGGGGGGGGGGGCCCCCCCCCCCCCCCCCCCCCC
> GGGGGGGGGGGGGGGGGGCCCCCCCCCCCCCCCCCCCCCC
> -50 kcal/mol
> ((((((((((((((((((...)))))))))))))))))).
GGGUGAUGUGGGUACGUCGGGGGGAAAAAACCCCCCCCCCCCCCAAAAAAAAAGGGGGGGGGUUUUUUUUCCCCCCCCCCCCCC
> GGGUGAUGUGGGUACGUCGGGGGGAAAAAACCCCCCCCCCCCCCAAAAAAAAAGGGGGGGGGUUUUUUUUCCCCCCCCCCCCCC
> -48.4 kcal/mol
> (((.(((((....)))))((((((((((((..(((((((((............)))))))))..))))))))))))..)))...
```

This program can be a bit slower than expected, as it is done by emulation via the logarithmic folding algorithm.

### Energy Calculation Programs
These programs take an RNA sequence and a structure and calculate the expected free energy of that pair using a particular energy model. All these programs expect a sequence of whitespace separated sequence structure pairs as input. Each pair should be a sequence and structure separated by whitespace. These programs will all output 3 lines for each input pair. The parsed primary sequence, the parsed secondary structure, and the free energy. Every energy calculation program expects the first argument on command line to be the path to the data_tables folder. This should be the relative path to the data_tables directory (you will find this in the root directory of this project). In addition, the logarithmic and linear (but not Aalberts & Nandagoapl) programs can take a second command line argument. If it is "-describe", then they will also print out a breakdown of the free energy for every loop in the structure.

#### Logarithmic Energy Calculator
This program (logarithmic_efn) scores energy using the logarithmic model (see above). Here is an example usage:

```
./the_build_directory/logarithmic_efn data-tables/
GGGUGAUGUGGGUACGUCGGGGGGAAAAAACCCCCCCCCCCCCCAAAAAAAAAGGGGGGGGGUUUUUUUUCCCCCCCCCCCCCC
(((.(((((....)))))((((((((((((..(((((((((............)))))))))..))))))))))))....))).
> GGGUGAUGUGGGUACGUCGGGGGGAAAAAACCCCCCCCCCCCCCAAAAAAAAAGGGGGGGGGUUUUUUUUCCCCCCCCCCCCCC
> (((.(((((....)))))((((((((((((..(((((((((............)))))))))..))))))))))))....))).
> -48.2 kcal/mol
```

#### Aalberts & Nandagopal Energy Calculator
This program (an_efn) scores energy using the Aalberts & Nandagopal model (see above). Here is an example usage:

```
./the_build_directory/an_efn data-tables/
GGGUGAUGUGGGUACGUCGGGGGGAAAAAACCCCCCCCCCCCCCAAAAAAAAAGGGGGGGGGUUUUUUUUCCCCCCCCCCCCCC
(((.(((((....)))))((((((((((((..(((((((((............)))))))))..))))))))))))..)))...
> GGGUGAUGUGGGUACGUCGGGGGGAAAAAACCCCCCCCCCCCCCAAAAAAAAAGGGGGGGGGUUUUUUUUCCCCCCCCCCCCCC
> (((.(((((....)))))((((((((((((..(((((((((............)))))))))..))))))))))))..)))...
> -49.8 kcal/mol
```

#### Linear Energy Calculator
This program (linear_efn) scores energy using the linear model (see above). Here is an example usage:

```
./the_build_directory/linear_efn data-tables/
GGGUGAUGUGGGUACGUCGGGGGGAAAAAACCCCCCCCCCCCCCAAAAAAAAAGGGGGGGGGUUUUUUUUCCCCCCCCCCCCCC
(((.(((((....)))))((((((((((((..(((((((((............)))))))))..))))))))))))..)))...
> GGGUGAUGUGGGUACGUCGGGGGGAAAAAACCCCCCCCCCCCCCAAAAAAAAAGGGGGGGGGUUUUUUUUCCCCCCCCCCCCCC
> (((.(((((....)))))((((((((((((..(((((((((............)))))))))..))))))))))))..)))...
> -48.4 kcal/mol
```

### A Note on data_tables
A common usage mistake is incorrectly providing the relative path to the data_tables directory to the executables. If you are getting strange looking structures with odd free energies, this is the most likely cause.

### Samples From the Paper
In the paper, two sample RNA folds are given (Figures 1 and 2). We shall reproduce them here to complete the tutorial. Obviously, with the full database, our full set of results could easily be reproduced the same way.


#### Figure 1
Consider the tRNA found in Figure 1 of the paper (Sprinzl ID RA1661):

```
GGGGGCAUAGCUCAGCUGGGAGAGCGCCUGCUUUGCACGCAGGAGGUCUGCGGUUCGAUCCCGCGCGCUCCCACCA
((((((...((((........)))).(((((.......))))).....(((((.......))))).))))))...
```

Here is how to reproduce the linear folding:

```
./the_build_directory/linear_fold data_tables/
GGGGGCAUAGCUCAGCUGGGAGAGCGCCUGCUUUGCACGCAGGAGGUCUGCGGUUCGAUCCCGCGCGCUCCCACCA
> GGGGGCAUAGCUCAGCUGGGAGAGCGCCUGCUUUGCACGCAGGAGGUCUGCGGUUCGAUCCCGCGCGCUCCCACCA
> -30.9 kcal/mol
> ((((((...((((........)))).(((((.......))))).....(((((.......))))).))))))....
```

And the logarithmic folding:

```
./the_build_directory/logarithmic_fold data_tables/
GGGGGCAUAGCUCAGCUGGGAGAGCGCCUGCUUUGCACGCAGGAGGUCUGCGGUUCGAUCCCGCGCGCUCCCACCA
> GGGGGCAUAGCUCAGCUGGGAGAGCGCCUGCUUUGCACGCAGGAGGUCUGCGGUUCGAUCCCGCGCGCUCCCACCA
> -30.2 kcal/mol
> ..((((...))))...((((((.((((.......((.(((((.....)))))))........)))).))))))...
```

Out of interest, lets get an energy breakdown of the linear model's result, but using the logarithmic model's energy calculator. This will score any multi-loops in the linear prediction as though they had a logarithmic energy function.

```
./the_build_directory/logarithmic_efn data_tables/ -describe
GGGGGCAUAGCUCAGCUGGGAGAGCGCCUGCUUUGCACGCAGGAGGUCUGCGGUUCGAUCCCGCGCGCUCCCACCA
((((((...((((........)))).(((((.......))))).....(((((.......))))).))))))....
> GGGGGCAUAGCUCAGCUGGGAGAGCGCCUGCUUUGCACGCAGGAGGUCUGCGGUUCGAUCCCGCGCGCUCCCACCA
> ((((((...((((........)))).(((((.......))))).....(((((.......))))).))))))....
> -30.1 kcal/mol
> External-loop: -17/-301
> AU/GU penalties: 0
> Stacking: -17
> Stacking(0, 71, 3' dangle) 
>  Two-loop closed by (0, 71) and (1, 70): -33/-284
>   Two-loop closed by (1, 70) and (2, 69): -33/-251
>    Two-loop closed by (2, 69) and (3, 68): -15/-218
>     Two-loop closed by (3, 68) and (4, 67): -21/-203
>      Two-loop closed by (4, 67) and (5, 66): -34/-182
>       Multi-loop closed by (5, 66): 9/-148
>       AU/GU penalties: 5
>       Multi-loop closure: 77
>       Multi-loop closure featues: branches=4 unpaired=10 
>      Stacking: -73
>       Stacking(9, 24, mismatch-mediated coax (26,42)) Stacking(5, 66, mismatch-mediated coax (48,64)) 
>       Two-loop closed by (9, 24) and (10, 23): -34/-39
>        Two-loop closed by (10, 23) and (11, 22): -21/-5
>         Two-loop closed by (11, 22) and (12, 21): -24/16
>          One-loop closed by (12, 21): 40/40
>        Two-loop closed by (26, 42) and (27, 41): -33/-63
>         Two-loop closed by (27, 41) and (28, 40): -21/-30
>          Two-loop closed by (28, 40) and (29, 39): -21/-9
>           Two-loop closed by (29, 39) and (30, 38): -34/12
>            One-loop closed by (30, 38): 46/46
>        Two-loop closed by (48, 64) and (49, 63): -14/-55
>         Two-loop closed by (49, 63) and (50, 62): -34/-41
>          Two-loop closed by (50, 62) and (51, 61): -24/-7
>           Two-loop closed by (51, 61) and (52, 60): -33/17
>            One-loop closed by (52, 60): 50/50
```

Notice that in the breakdown, every loop is listed with recursive indentation. The free energy contribution of the loop, followed by the recursive energy of the closed substructure, is given after each loop. Note, a minor oddity of the energy breakdown is that energies are reported as multiples of 0.1 kcal/mol.

#### Figure 2

Now, consider this signal recognition particle RNA from Figure 2:

```
GGUCCCCUCGCAACGAUCAGCCGUGAACCCGGUCAGGCCCGGAAGGGAGCAGCCGCAGCGGUGACAUUGUGUGCCGGGGUGUGGCUGGUAG
((((.((((((.(((((..(((((.....((((....(((....)))....))))..)))))...))))))...)))))...)))).....
```

Here is the linear fold invocation:

```
./the_build_directory/linear_fold data_tables/
GGUCCCCUCGCAACGAUCAGCCGUGAACCCGGUCAGGCCCGGAAGGGAGCAGCCGCAGCGGUGACAUUGUGUGCCGGGGUGUGGCUGGUAG
> GGUCCCCUCGCAACGAUCAGCCGUGAACCCGGUCAGGCCCGGAAGGGAGCAGCCGCAGCGGUGACAUUGUGUGCCGGGGUGUGGCUGGUAG
> -32.9 kcal/mol
> ((((.((((((((((((..(((((.....((((....(((....)))....))))..)))))...))))).)).)))))...)))).....
```

And now the Aalberts & Nandagopal fold invocation (this can take a while):
```
./the_build_directory/an_fold data_tables/
GGUCCCCUCGCAACGAUCAGCCGUGAACCCGGUCAGGCCCGGAAGGGAGCAGCCGCAGCGGUGACAUUGUGUGCCGGGGUGUGGCUGGUAG
> GGUCCCCUCGCAACGAUCAGCCGUGAACCCGGUCAGGCCCGGAAGGGAGCAGCCGCAGCGGUGACAUUGUGUGCCGGGGUGUGGCUGGUAG
> -33.4 kcal/mol
> ..((((.((((...........))))..((((......))))..)))).((((((((.((((.(((...)))))))...))))))))....
```

Now let's check what the energy calculator for the linear model says for the Aalberts & Nandagopal model's prediction:

```
./the_build_directory/linear_efn data_tables/
GGUCCCCUCGCAACGAUCAGCCGUGAACCCGGUCAGGCCCGGAAGGGAGCAGCCGCAGCGGUGACAUUGUGUGCCGGGGUGUGGCUGGUAG
..((((.((((...........))))..((((......))))..)))).((((((((.((((.(((...)))))))...))))))))....
> GGUCCCCUCGCAACGAUCAGCCGUGAACCCGGUCAGGCCCGGAAGGGAGCAGCCGCAGCGGUGACAUUGUGUGCCGGGGUGUGGCUGGUAG
> ..((((.((((...........))))..((((......))))..)))).((((((((.((((.(((...)))))))...))))))))....
> -32.3 kcal/mol
```

### Using Alternate Parameters
The efn and folding programs can be modified to use non-standard multi-loop parameters. Both the AalbertsModel and NNUnpairedModel classes (see the source code in the programs/ folder for example usages of these classes) have a method called 'SetMLParams' that can be used to change the multi-loop energy model parameters. For example, the programs for the linear model use this method to simulate the linear model using a reparameterized logarithmic model.

## Parameter Training Algorithms
The parameter training programs are train_linear, train_logarithmic, and train_an. They all have similar input requirements. Instructions for flags can be found by calling a program with the flag "-h".

### Basic Usage
All parameter training algorithms need three pieces of input to function. The path to the data tables, which are the parameters for the entire nearest neighbor model as per RNAstructure; the path to the data set, which is a folder containing all the .ct files for RNAs used in training; a list of .ct file names to use for training. If run from the root directory of this repository, the default paths for the data tables and data set will work out of the box. To change them, use the flag "-h" to see how. The list of .ct files to use should be given via standard input. The folder "data_sets/" contains several .ctset files. These are the lists of .ct files used for the paper. The folder "data_sets/ct_files" contains the archive of training data used in the paper.

The range of parameters to optimize over is hard coded in each program. They are those used in the paper. To change the parameter ranges, the code for a program must be changed.

### An Example
Let us consider an example run of parameter optimization. For this example, I will assume that we are in the root directory of the project. Also, I will assume that "path_to_programs/" is the path to the folder the programs are in. The following should bring up usage information for the linear parameter training program.

```
./path_to_programs/train_linear -h
```

It should look a bit like this:

```
Trains the parameters of the linear model using IBF.
Usage:
  Train linear model using IBF [OPTION...]

  -d, --data_path arg  Path to data_tables (default: data_tables/)
  -c, --ct_path arg    Path to the folder of CTs (default:
                       data_set/ct_files/)
  -t, --threads arg    Number of threads to use (default: 16)
  -h, --help           Print help
```

Now, let us run the program using the small training set. By default, this will use all the cores on your machine, so be wary! The following invocation should starting training:

```
./bin/programs/train_linear < data_set/training_small.ctset
```

It can take a while. You should see something like the following as a result:

```
Seed #1 init = 32 branch = -58 unpaired = 3: 0.101219
Seed #2 init = 112 branch = -38 unpaired = -4: 0.472647
Seed #3 init = 183 branch = -9 unpaired = -16: 0.096862
Seed #4 init = 125 branch = -44 unpaired = 2: 0.463108
Seed #5 init = 193 branch = 5 unpaired = -5: 0.464895
Epoch #0: 
	Average F-Score = 0.0342677
	Best score = 0.768073
	init = 109 branch = 30 unpaired = 11
Epoch #1: 
	Average F-Score = 0.417572
	Best score = 0.630142
	init = 75 branch = -7 unpaired = 0
Epoch #2: 
	Average F-Score = 0.617731
	Best score = 0.636124
	init = 163 branch = -25 unpaired = -2
Epoch #3: 
	Average F-Score = 0.590204
	Best score = 0.629207
	init = 115 branch = -15 unpaired = 0
Epoch #4: 
	Average F-Score = 0.631208
	Best score = 0.633857
	init = 111 branch = -18 unpaired = 1
Epoch #5: 
	Average F-Score = 0.625122
	Best score = 0.631451
	init = 110 branch = -10 unpaired = -1
Epoch #6: 
	Average F-Score = 0.632576
	Best score = 0.635146
	init = 124 branch = -14 unpaired = -1
Epoch #7: 
	Average F-Score = 0.631371
	Best score = 0.633531
	init = 109 branch = -13 unpaired = 0
Epoch #8: 
	Average F-Score = 0.633684
	Best score = 0.634308
	init = 109 branch = -13 unpaired = 0
Best parameters: init = 109 branch = -13 unpaired = 0
```

Each epoch is a training set. We can see the F-Score steadily increases until we converge on a parameter set, which is given at the end. Note that the parameters found are different from those in the paper, since we training on the small data rather than the large.

The same training process will work for any of the parameter optimization programs. Please keep in mind that parameter training uses all cores and can therefore use a lot of memory!

