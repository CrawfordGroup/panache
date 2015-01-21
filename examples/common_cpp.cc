#include <iostream>
#include <fstream>
#include <sstream>

#include "panache/Output.h"
#include "panache/c_convert.h"
#include "panache/SimpleMatrix.h"

using namespace panache;
using namespace std;

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

