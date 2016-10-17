
#if !defined(ALGORITHM_H)
#define ALGORITHM_H

#include "arrayclass.h"
#include "forceclass.h"
#include "dotarray.h"
#include "rna_library.h"

#ifdef _WINDOWS_GUI
#include "../RNAstructure_windows_interface/TProgressDialog.h"
#else

#ifdef _JAVA_GUI
#include "../java_interface/SWIG/TProgressDialog.h"
#else
#include "TProgressDialog.h"
#endif // JAVA GUI

#endif //WINDOWS

//***********************************Structures:

/////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////

void de_allocate (int **v,int i);//deallocates memory for a 2d array
void de_allocate (bool **v,int i);//alternative form of de_allocate
void de_allocate (short int **v,int i);//alternative form of de_allocate



//**********************************prototypes:



void getout (char *energyfile);//get name of file to output
										//	energy info

void efn2(datatable *data,structure *ct, int structnum = 0, bool simplemb = false, const char *outputfile=NULL);//energy calculator


void energyout(structure *ct,char *enrgyfile);

//dynamic programming algorithm for secondary structure prediction by free energy minimization
	//this is the dynamic folding algorithm of Zuker
         //cntrl6 = #tracebacks
         //cntrl8 = percent sort
         //cntrl9 = window
		//TProgressDialog is an interface for returning the progress of the calculation
		//Savfile is for creating a file with arrays and parameters for refolding with different 
			//suboptimal tracebacks
		//quickenergy indicates whether to determine the lowest free energy for the sequence without a structure
		//quickstructure is a bool that will generate only the lowest free energy structure.  No savefiles can generated. 
		//maxinter is the maximum number of unpaired nucleotides allowed in an internal loop
	//This returns an error code, where zero is no error and non-zero indicates a traceback error.
int dynamic (structure *ct,datatable *data,int cntrl6,int cntrl8,int cntrl9,
	TProgressDialog* update=0, bool quickenergy = false, char* savfile = 0, int maxinter = 30, bool quickstructure = false, bool simple_iloops = true);


void fill(structure *ct, arrayclass &v, arrayclass &w, arrayclass &wmb, forceclass &fce, int &vmin,bool *lfce, bool *mod,
          integersize *w5, integersize *w3, bool quickenergy,
          datatable *data, arrayclass *w2, arrayclass *wmb2, arrayclass *we,TProgressDialog* update = 0, int maxinter = 30, bool quickstructure = false, bool simple_iloops = true);


//The fill step of the dynamic programming algorithm for free energy minimization:
void fill(structure *ct, arrayclass &v, arrayclass &w, arrayclass &wmb, forceclass &fce, int &vmin,bool *lfce, bool *mod,
		  integersize *w5, integersize *w3, bool qickenergy,
		  datatable *data, arrayclass *w2, arrayclass *wmb2, TProgressDialog* update=0, int maxinter = 30, bool quickstructure = false, bool simple_iloops = true);

//this overloaded dynamic function is used by NAPSS program to generate a special format dotplot
void dynamic (structure *ct,datatable* data,int cntrl6, int cntrl8,int cntrl9,
              arrayclass *v, arrayclass *vmb/*tracks MB loops*/, arrayclass *vext/*tracks exterior loops*/,
              TProgressDialog* update=0, bool quickenergy = false, char* savefile = 0, int maxinter = 30, bool quickstructure = false, bool simple_iloops = true);
//this overloaded fill function is used to NAPSS program to generate a special format dotplot
void fill(structure *ct, arrayclass &v, arrayclass &vmb, arrayclass &vext, arrayclass &w, arrayclass &wmb, forceclass &fce, 
          int &vmin, bool *lfce, bool *mod,integersize *w5, integersize *w3, bool quickenergy,
          datatable *data, arrayclass *w2, arrayclass *wmb2, TProgressDialog* update=0, int maxinter = 30, bool quickstructure = false, bool simple_iloops = true);

void errmsg(int err,int err1);//function for outputting info in case of an error
void update (int i);//function informs user of progress of fill algorithm

//the force... functions are used to initialize the arrays used to apply
//constraints during the folding process
void forcepair(int x,int y,structure *ct,forceclass *v);
void forcesingle(int x,structure* ct,forceclass *v);
void forcedbl(int dbl,structure* ct,forceclass *w,bool *v);
void forceinter(int dbl,structure* ct,forceclass *w);
void forceinterefn(int dbl,structure* ct,forceclass *w);

//filter is used to choose structures to output after efn2
//	this can make the number of structures more reasonable for inspection
//	it takes a structure, ct, which also contains the final output,
//	percent sort, maximum number of structures, and window size
void filter(structure* ct, int percent, int max, int window);

//force is used to prepare arrays for the function dynamic, used during the
//	fill routines - it coordinates the force...() functions above
void force(structure *ct,forceclass *fce, bool *lfce);

//Use the fill information to generate a set of suboptimal structures using the mfold heuristic.
	//This returns an error code, where zero is no error and non-zero indicates a traceback error.
int traceback(structure *ct, datatable *data, arrayclass *v, arrayclass *w, arrayclass *wmb, arrayclass *w2,arrayclass *wmb2, integersize *w3, integersize *w5, forceclass *fce,
	bool *lfce,integersize vmin, int cntrl6, int cntrl8, int cntrl9, bool *mod);

//this function is used to calculate the values of all the dots in a dot plot
void dotefn2(structure *ct, datatable *data, arrayclass *v, arrayclass *w, arrayclass *w2,
	int *w3, int *w5, short int **fce, bool *lfce,int vmin,dotarray *dots,
   TProgressDialog* PD = 0);
void calcpnum(dotarray *dots, int *pnum, int increment, short int numofbases,
	TProgressDialog *PD = 0);
void savefile(int i, std::ofstream* sav);//this function is used to make a save file
											//after the fill algorithm
short int readfile(std::ifstream *read);//this function is used to read save files
void savedot(dotarray *dots,structure *ct, char *filename); //save dot plot info
void readdot(dotarray *dots, structure *ct, char *filename);//read a dot plot file
void dpalign(dotarray *dots1,dotarray *dots2,structure* ct1,structure *ct2,short int *align);
short int getbestdot(dotarray *dots1,dotarray *dots2, structure* ct1,
	structure *ct2, short int i, short int j);//return the best dot for base i
   //in dots1 and j in dots2
//dpalign will align two dot plots and store the info in the array align
void energydump (structure *ct, arrayclass *v, datatable *data, int n,char *filename, int i, int j);
//energydump will spit out the composite free energies for a traceback
void energydump (structure *ct, datatable *data,arrayclass *v, int n,char *filename);
//energydump2 will spit out the composite free energies for a traceback -- with
//the au penalty associated with the correct entity
int checknp(bool lfce1,bool lfce2); 
//this function is used by the fill and trace to check whether nucleotides 
//contained in a dangling end are forced double stranded

void opensav(char* filename, structure* ct, int cntrl6, int cntrl8,int cntrl9);//opens a save file with information filled by
   									//fill algorithm

#ifdef DYNALIGN_II
int trace(structure *ct, datatable *data, int ii, int ji,
          arrayclass *v, arrayclass *w, arrayclass *wmb, arrayclass *w2, arrayclass *wmb2, 
          bool *lfce, forceclass *fce, integersize *w3, integersize *w5,bool *mod,arrayclass *we = NULL,integersize energy = 0,short open = 0, short pair = 0, bool quickstructure = false);
#else
int trace(structure *ct, datatable *data, int ii, int ji,
		arrayclass *v, arrayclass *w, arrayclass *wmb, arrayclass *w2, arrayclass *wmb2, 
          bool *lfce, forceclass *fce, integersize *w3, integersize *w5,bool *mod,bool quickstructure=false);
#endif
void readsav(const char *filename, structure *ct, arrayclass *w2, arrayclass *wmb2, 
			 integersize *w5, integersize *w3, bool *lfce, bool *mod, datatable *data,
			 arrayclass *v, arrayclass *w, arrayclass *wmb, forceclass *fce, int *vmin);
bool notgu(int i, int j, structure *ct);
void writehelixfile(char *filename,structure *ct,int StructureNumber);//write a file with the helices

//reads the save file

void cctout(structure *ct, char *filename);
//outputs a structure in cct format, which is shortened format compared to ct.

#endif
