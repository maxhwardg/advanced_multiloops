/*
 * Dynalign
 *  by David Mathews, copyright 2002, 2003,2004,2005
 *
 * Contributors: Chris Connett and Andrew Yohn, 2006
 *               Michael Sloma, 2015
 */

#include "arrayclass.h"

#include "defines.h"

arrayclass::arrayclass(int size, integersize energy) {
	infinite = energy;
	Size = size;
	int i,j;
	dg = new integersize *[size+1];

	for (i=0;i<=(size);i++){
		dg[i] = new integersize [size+1];
	}
	for (i=0;i<=size;i++){
		for (j=0;j<size+1;j++){
			dg[i][j] = INFINITE_ENERGY;
		}
	}

	//Now move pointers, to facilitate fast access by
    //avoiding arithmetic during array access function
	for (i=0;i<=size;++i){
		dg[i]-=i;
	}
    //columns are i index, rows are j index
    //because i>j, array is now shaped like this:
    //     n
    // |    |
    //  |    |
    //   |    |
    //    |    |
    //         2n
    //j>n means this is an exterior fragment

}

// the destructor deallocates the space used
arrayclass::~arrayclass(){
	for (int i=0;i<=Size;i++){
		//move pointers back before deleting
		dg[i]+=i;

		//now delete
		delete[] dg[i];
	}
	delete[] dg;
}

arrayclassT::arrayclassT(int size, integersize energy) {
	infinite = energy;
	Size = size;
	dg = new integersize *[2*size+1];

    //because this is transpose of arrayclass,
    //i index is rows, j index is columns
    //array is shaped like this:
    // ||
    // | |
    // |  |
    // |   |  n
    //  |   |
    //   |  |
    //    | |
    //     || 2n
    //      n
	for (int i=0;i<=2*size;i++){
        int rowlength = i<=size ? i+1 : 2*size+1-i;
		dg[i] = new integersize [rowlength];
        for (int j=0;j<rowlength;j++){
            dg[i][j] = infinite;
        }
	}
    //move pointer so we don't have to do arithmetic during array access
    for(int i = size+1;i<=2*size;i++){
        dg[i] -= i-size;
    }
}

// the destructor deallocates the space used
arrayclassT::~arrayclassT(){
	for (int i=0;i<=2*Size;i++){
        if (i>Size){
            dg[i] += (i-Size);
        }
		delete[] dg[i];
	}
	delete[] dg;
}
