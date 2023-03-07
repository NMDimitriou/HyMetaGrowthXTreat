#include "declarations.cuh"

/* Reads the model parameters from params.txt file
 * Nikolaos Dimitriou, McGill 2021
*/ 

void readData(InitialData& id, char* argv)
{
#if ( MODEL == FK )
	ifstream f;
    f.open(argv);
    if (f.is_open()){
        f >> id.da;
        f >> id.s;
        f >> id.ka;
        f >> id.param;
        f >> id.dev;
        f >> id.pid;
        f.close();
    }

#elif ( MODEL == KSC || MODEL == TEST)
	ifstream f;
    f.open(argv);
    if (f.is_open()){
        f >> id.da;
        f >> id.s;
        f >> id.ka;
		f >> id.chi;
		f >> id.db;
		f >> id.param;
		f >> id.dev;
		f >> id.pid;
        f.close();
    }
	id.r=0.;
#elif ( MODEL == KSMD )
	ifstream f;
    f.open(argv);
    if (f.is_open()){
        f >> id.da;
        f >> id.s;
        f >> id.ka;
        f >> id.chi_ecm;
        f >> id.dc;
        f >> id.param;
        f >> id.dev;
        f >> id.pid;
        f.close();
    }
#elif ( MODEL == KSCMD )
	ifstream f;
    f.open(argv);
    if (f.is_open()){
        f >> id.da;
        f >> id.s;
        f >> id.ka;
        f >> id.chi;
		f >> id.chi_ecm;
        f >> id.db;
        f >> id.param;
        f >> id.dev;
        f >> id.pid;
        f.close();
    }
#endif
}
