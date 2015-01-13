#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <array>
#include <utility>
#include <cmath>
#include <algorithm>

#include "panache/DFTensor.h"
#include "panache/CHTensor.h"
#include "panache/Output.h"
#include "panache/c_convert.h"
#include "panache/SimpleMatrix.h"

using namespace panache;
using namespace std;

#include "panache/Parallel.h"

#ifdef PANACHE_MPI
#include <mpi.h>
#endif


#define CHOLESKY_DELTA 1e-4


ostream * out;

void PrintUsage(void)
{
    *out << "\n"
         << "Libpanache C++ example\n"
         << "\n"
         << "Usage: runtest [opt] <dir>\n"
         << "\n"
         << "Options:\n"
         << "-v           Verbose printing\n"
         << "-d           Write Q tensors to disk (rather than in core)\n"
         << "-c           Use Cyclops (MPI)\n" 
         << "-h           Print help (you're looking at it\n"
         << "<dir>        Directory holding the test information\n"
         << "\n\n";
}


// Note - exceptions are turned on for the ifstream object
// so that any parsing errors just throw an exeption. Catch those,
// and throw an exception
SharedBasisSet ReadBasisFile(SharedMolecule mol, const string & filename)
{
    ifstream f(filename.c_str());

    if(!f.is_open())
        throw std::runtime_error("Cannot open basis file!");

    f.exceptions(std::ifstream::failbit |
                 std::ifstream::badbit  |
                 std::ifstream::eofbit);

    try
    {
        int ncenters, nshells, nbf, nao, nprim;
        f >> nshells >> nbf >> nao >> nprim >> ncenters;

        unique_ptr<int[]> nshellspercenter(new int[ncenters]);
        unique_ptr<C_ShellInfo[]> shells(new C_ShellInfo[nshells]);

        for(int i = 0; i < ncenters; i++)
            f >> nshellspercenter[i];

        int dum; // holds the shell index. Not needed

        for(int i = 0; i < nshells; i++)
        {
            f >> dum >> shells[i].nprim >> shells[i].am >> shells[i].ispure;
            shells[i].coef = new double[shells[i].nprim];
            shells[i].exp = new double[shells[i].nprim];

            for(int j = 0; j < shells[i].nprim; j++)
                f >> shells[i].coef[j] >> shells[i].exp[j];
        }


        // Create the basis set
        // These are not normalized
        auto ret = BasisSetFromArrays(mol, ncenters, nshellspercenter.get(), shells.get());

        for(int i = 0; i < nshells; i++)
        {
            delete [] shells[i].exp;
            delete [] shells[i].coef;
        }

        return ret;
    }
    catch(...)
    {
        throw std::runtime_error("Error parsing basis file");
    }
}




SharedMolecule ReadMoleculeFile(const string & filename)
{
    ifstream f(filename.c_str());

    if(!f.is_open())
        throw runtime_error("Cannot open file!");

    f.exceptions(std::ifstream::failbit |
                 std::ifstream::badbit  |
                 std::ifstream::eofbit);
    try
    {
        int natoms, nallatoms;
        std::string schoen, fullpg; // fullpg is no longer used

        f >> natoms >> nallatoms;
        f.ignore(); // remove newline

        getline(f, schoen);
        getline(f, fullpg);

        unique_ptr<C_AtomCenter[]> atoms(new C_AtomCenter[nallatoms]);

        double dummy; // for mass and Z-number

        for(int i = 0; i < nallatoms; i++)
        {
            f.get(atoms[i].symbol, 4, ' ');
            atoms[i].symbol[4] = '\0';

            f >> dummy  // Z is not used
              >> atoms[i].center[0]
              >> atoms[i].center[1]
              >> atoms[i].center[2]
              >> dummy;   // mass is not used

            f.ignore(); // ignore the newline
        }

        auto ret = MoleculeFromArrays(nallatoms, atoms.get());

        return ret;
    }
    catch(...)
    {
        throw std::runtime_error("Error parsing molecule file");
    }
}



/*!
 * \brief Reads the C matrix file
 */
std::shared_ptr<SimpleMatrix> ReadCMatrixFile(const string & filename)
{
    ifstream f(filename.c_str());

    if(!f.is_open())
        throw runtime_error("Cannot open matrix file!");

    f.exceptions(std::ifstream::failbit |
                 std::ifstream::badbit  |
                 std::ifstream::eofbit);
    try
    {
        int nrow, ncol;

        f >> nrow >> ncol;

        std::shared_ptr<SimpleMatrix> mat(new SimpleMatrix(nrow, ncol));

        for(int i = 0; i < nrow; i++)
        for(int j = 0; j < ncol; j++)
            f >> (*mat)(i,j);

        return mat;
    }
    catch(...)
    {
        throw runtime_error("Error parsing matrix file");
    }
}


/*!
 * \brief Reads the orbital energy file
 */
std::vector<double> ReadOrbEnFile(const string & filename)
{
    ifstream f(filename.c_str());

    if(!f.is_open())
        throw runtime_error("Cannot open matrix file!");

    f.exceptions(std::ifstream::failbit |
                 std::ifstream::badbit |
                 std::ifstream::eofbit);
    try
    {
        int n;
        f >> n;
        std::vector<double> vec(n);

        for(int i = 0; i < n; i++)
            f >> vec[i];

        return vec;
    }
    catch(...)
    {
        throw runtime_error("Error parsing matrix file");
    }
}


int ReadNocc(const string & filename)
{
    ifstream f(filename.c_str());

    if(!f.is_open())
        throw runtime_error("Cannot open nocc file!");

    f.exceptions(std::ifstream::failbit |
                 std::ifstream::badbit  |
                 std::ifstream::eofbit);
    try
    {
        int nocc;

        f >> nocc;

        return nocc;
    }
    catch(...)
    {
        throw runtime_error("Error parsing nocc file");
    }
}




string GetNextArg(int & i, int argc, char ** argv)
{
    if(i >= argc)
        throw runtime_error("Error - no more arguments!");

    return argv[i++];
}

int GetIArg(int & i, int argc, char ** argv)
{
    string str = GetNextArg(i, argc, argv);
    try {

        return stoi(str);
    }
    catch(...)
    {
        stringstream ss;
        ss << "Cannot convert to int: " << str;
        throw runtime_error(ss.str());
    }
}



int main(int argc, char ** argv)
{

    int ret = 0;

    #ifdef PANACHE_MPI
    //MPI_Init(&argc, &argv);
    int threadprov;
    MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &threadprov);
    //std::cout << "Wanted " << MPI_THREAD_MULTIPLE << " got " << threadprov << "\n";
    #endif
    panache::parallel::Init(&argc, &argv);


    // if not master, dump output to a string stream
    if(panache::parallel::IsMaster())
        out = &cout;
    else
        out = new stringstream;

    try
    {

        if(argc == 1)
        {
            PrintUsage();
            return 0;
        }


        string dir;

        // program options
        bool verbose = false;
        bool cyclops = false;
        bool disk = false;
        bool keepdisk = false;
        bool readdisk = false;

        int i = 1;
        while(i < argc)
        {
            string starg(GetNextArg(i, argc, argv));
            if(starg == "-v")
                verbose = true;
            else if(starg == "-c")
                cyclops = true;
            else if(starg == "-d")
                disk = true;
            else if(starg == "-k")
                keepdisk = true;
            else if(starg == "-r")
                readdisk = true;
            else if(starg == "-h")
            {
                PrintUsage();
                return 0;
            }
            else if(dir.length() == 0)
            {
                // add trailing slash if needed
                if(starg[starg.length()-1] != '/')
                    starg.append("/");
                dir = starg;
            }
            else
            {
                *out << "\n\n";
                *out << "-----------------------------------------------\n";
                *out << "Error: Unknown argument \"" << starg << "\"\n";
                *out << "-----------------------------------------------\n\n";
                return 1;
            }
        }

        // check for incompatible options
        #ifndef PANACHE_CYCLOPS
        if(cyclops)
            throw std::runtime_error("Error - Panache not compiled with cyclops!");
        #else
        if(cyclops && disk)
            throw std::runtime_error("Incompatible options: cyclops and disk");
        #endif

        if(verbose)
            panache::output::SetOutput(&*out);

        // Get the test description
        string desc_path(dir);
        desc_path.append("desc");

        fstream desc(desc_path.c_str());
        if(!desc.is_open())
        {
            *out << "\nFatal Error: Cannot open " << desc_path << "\"\n\n";
            return -1;
        }

        string desc_str;
        getline(desc, desc_str);

        *out << "--------------------------------------\n";
        *out << "Example: MP2 for " << desc_str << "\n";
        *out << "--------------------------------------\n";

        auto mol = ReadMoleculeFile(dir + "geometry");
        auto primary = ReadBasisFile(mol, dir + "basis.primary");
        auto cmat = ReadCMatrixFile(dir + "cmat");
        auto orben = ReadOrbEnFile(dir + "orben");

        int nso = primary->nbf();
        int nocc = ReadNocc(dir + "nocc");
        int nmo = nso;
        int nvir = nmo - nocc;

        DFTensor dft(primary, dir + "basis.aux.gbs", "/tmp/df",
                     DFOPT_COULOMB | DFOPT_EIGINV,
                     BSORDER_PSI4, 0);

        // *** We are only testing Qso from CHTensor               *** //
        // *** But generating them all (to test for memory issues) *** //
        CHTensor cht(primary, CHOLESKY_DELTA, "/tmp/ch", BSORDER_PSI4, 0);

        dft.SetCMatrix(cmat->pointer(), nmo, false);
        cht.SetCMatrix(cmat->pointer(), nmo, false);
        dft.SetNOcc(nocc);
        cht.SetNOcc(nocc);

        int dfqflags = (QGEN_QSO | QGEN_QOV);
        int chqflags = (QGEN_QSO | QGEN_QOV);

        int qstore = 0;

        if(disk)
            qstore |= QSTORAGE_ONDISK;
        if(keepdisk)
            qstore |= QSTORAGE_KEEPDISK;
        if(readdisk)
            qstore |= QSTORAGE_READDISK;
        #ifdef PANACHE_CYCLOPS
        else if(cyclops)
            qstore |= QSTORAGE_CYCLOPS;
        #endif
        else
            qstore |= QSTORAGE_INMEM;

        dft.GenQTensors(dfqflags, qstore);
        cht.GenQTensors(chqflags, qstore);


        int naux, ndim1, ndim2;

        // ignoring return value
        dft.TensorDimensions(QGEN_QOV, naux, ndim1, ndim2);


        // should be occ x vir
        if(ndim1 != nocc || ndim2 != nvir)
          throw std::runtime_error("Inconsistent dimensions!");

        // could use a plain pointer here, but
        // I like smart pointers
        std::unique_ptr<double[]> buf_ia(new double[naux]);
        std::unique_ptr<double[]> buf_jb(new double[naux]);
        std::unique_ptr<double[]> buf_ib(new double[naux]);
        std::unique_ptr<double[]> buf_ja(new double[naux]);
        double * buf_iap = buf_ia.get();
        double * buf_jbp = buf_jb.get();
        double * buf_ibp = buf_ib.get();
        double * buf_jap = buf_ja.get();

        // calculate MP2 with density fitting
        double e2s = 0;
        double e2t = 0;

        for(int i = 0; i < nocc; i++)
        for(int a = 0; a < nvir; a++)
        {
            int ia = dft.CalcIndex(QGEN_QOV, i, a);
            dft.GetBatch(QGEN_QOV, buf_iap, naux, ia);

            for(int j = 0; j < nocc; j++)
            for(int b = 0; b < nvir; b++)
            {
                int jb = dft.CalcIndex(QGEN_QOV, j, b);
                dft.GetBatch(QGEN_QOV, buf_jbp, naux, jb);
    
                int ib = dft.CalcIndex(QGEN_QOV, i, b);
                int ja = dft.CalcIndex(QGEN_QOV, j, a);
    
                dft.GetBatch(QGEN_QOV, buf_ibp, naux, ib);
                dft.GetBatch(QGEN_QOV, buf_jap, naux, ja);
    
                double iajb = 0;
                double ibja = 0;
    
                // contract
                for(int q = 0; q < naux; q++)
                {
                    iajb += buf_iap[q]*buf_jb[q];
                    ibja += buf_ibp[q]*buf_ja[q];
                }
                
                // remember that a,b goes from [0,nvir) so we have to 
                // factor in the occupied orbitals for orben
                double denom = orben[i]+orben[j]-orben[nocc+a]-orben[nocc+b]; 
                e2s += (iajb*iajb)/denom;
                e2t += (iajb*(iajb-ibja))/denom;
            }
        }
       

        *out << "MP2 Correlation Energy\n"
             << "             Same spin: " << e2t << "\n"
             << "         Opposite spin: " << e2s << "\n"
             << "                 TOTAL: " << e2s+e2t << "\n"; 

        // even if not verbose, print the timints
        if(!verbose)
            panache::output::SetOutput(&*out);

        dft.PrintTimings();
        cht.PrintTimings();

    }
    catch(const exception & ex)
    {
        *out << "\n*****************"
             << "\nCAUGHT EXCEPTION!"
             << "\n" << ex.what()
             << "\n*****************"
             << "\n\n";

        ret = -1;
    }


    panache::parallel::Finalize();
    #ifdef PANACHE_MPI
    MPI_Finalize();
    #endif

    if(!panache::parallel::IsMaster())
        delete out;

    return ret;
}





