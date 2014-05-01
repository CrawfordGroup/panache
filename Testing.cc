#include <stdexcept>
#include <fstream>
#include <sstream>

#include "Testing.h"
#include "c_convert.h"

using std::string;
using std::ifstream;
using std::stringstream;


namespace panache
{
namespace testing
{

int TestInfo::ReadBasisFile(const string & filename,
                            int * &nshellspercenter,
                            C_ShellInfo * &shells,
                            TestInfo::BasisTest & test)
{
    ifstream f(filename.c_str());
    if(!f.is_open())
        throw std::runtime_error("Cannot open testing file!");

    int ncenters, nshells;
    f >> nshells >> ncenters;

    test.total_nshell.set(nshells);
    test.total_ncenters.set(ncenters);

    //std::cout << "Read: nshells = " << nshells << " ncenters = " << ncenters << "\n";

    nshellspercenter = new int[ncenters];
    for(int i = 0; i < ncenters; i++)
    {
        f >> nshellspercenter[i];
        test.center_nshell.push_back(TestResult<int>(nshellspercenter[i]));
        //std::cout << "   nshellspercenter[" << i << "] = " << nshellspercenter[i] << "\n";
    }


    shells = new C_ShellInfo[nshells];
    int dum; // holds the shell index. Not needed
    int totalnprim = 0;

    for(int i = 0; i < nshells; i++)
    {
        f >> dum >> shells[i].nprim >> shells[i].am >> shells[i].ispure;
        test.shell_nprim.push_back(TestResult<int>(shells[i].nprim));

        /*
        std::cout << "Shell " << i << " nprim = " << shells[i].nprim
                                   << " am = " << shells[i].am
                                   << " ispure = " << shells[i].ispure << "\n";
        */

        shells[i].coef = new double[shells[i].nprim];
        shells[i].exp = new double[shells[i].nprim];

        for(int j = 0; j < shells[i].nprim; j++)
        {
            f >> shells[i].coef[j] >> shells[i].exp[j];
            test.exp.push_back(TestResult<double>(shells[i].exp[j], TEST_EXPONENT_THRESHOLD));
            test.coef.push_back(TestResult<double>(shells[i].coef[j], TEST_COEFFICIENT_THRESHOLD));
            totalnprim++;
            //std::cout << "     exp: " << shells[i].exp[j] << " coef " << shells[i].coef[j] << "\n";
        }

    }

    test.total_nprim.set(totalnprim);

    return ncenters;
}

int TestInfo::ReadMolecule(const string & filename, C_AtomCenter * &atoms, MoleculeTest & test)
{
    ifstream f(filename.c_str());
    if(!f.is_open())
        throw std::runtime_error("Cannot open testing file!");

    int ncenters;
    f >> ncenters;

    test.total_ncenters.set(ncenters);

    atoms = new C_AtomCenter[ncenters];

    for(int i = 0; i < ncenters; i++)
    {
        f.ignore(); // ignore the newline
        char * s = new char[8];
        f.get(s, 4, ' ');
        atoms[i].symbol = s;

        f >> atoms[i].Z
          >> atoms[i].center[0]
          >> atoms[i].center[1]
          >> atoms[i].center[2]
          >> atoms[i].mass;

        test.xyz[0].push_back(TestResult<double>(atoms[i].center[0], TEST_MOLECULE_XYZ_THRESHOLD));
        test.xyz[1].push_back(TestResult<double>(atoms[i].center[1], TEST_MOLECULE_XYZ_THRESHOLD));
        test.xyz[2].push_back(TestResult<double>(atoms[i].center[2], TEST_MOLECULE_XYZ_THRESHOLD));
        test.Z.push_back(TestResult<double>(atoms[i].Z, TEST_MOLECULE_Z_THRESHOLD));
        test.mass.push_back(TestResult<double>(atoms[i].mass, TEST_MOLECULE_MASS_THRESHOLD));
    }

    return ncenters;
}


TestInfo::TestInfo(const std::string & testname, const std::string & dir)
    : testname_(testname)
{
    string primary_basis_filename(dir);
    string aux_basis_filename(dir);
    string molecule_filename(dir);

    primary_basis_filename.append("basis.primary");
    aux_basis_filename.append("basis.aux");
    molecule_filename.append("geometry");

    ncenters_ = ReadMolecule(molecule_filename, atoms_, molecule_test_);
    ReadBasisFile(primary_basis_filename, primary_nshellspercenter_, primary_shells_, primary_test_);
    ReadBasisFile(aux_basis_filename, aux_nshellspercenter_, aux_shells_, aux_test_);

    primary_nshells_ = 0;
    aux_nshells_ = 0;
    for(int i = 0; i < ncenters_; i++)
    {
        primary_nshells_ += primary_nshellspercenter_[i];
        aux_nshells_ += aux_nshellspercenter_[i];
    }
}



TestInfo::~TestInfo()
{
    for(int i = 0; i < primary_nshells_; i++)
    {
        delete [] primary_shells_[i].exp;
        delete [] primary_shells_[i].coef;
    }
    for(int i = 0; i < aux_nshells_; i++)
    {
        delete [] aux_shells_[i].exp;
        delete [] aux_shells_[i].coef;
    }

    delete [] primary_shells_;
    delete [] aux_shells_;
    delete [] primary_nshellspercenter_;
    delete [] aux_nshellspercenter_;

    for(int i = 0; i < ncenters_; i++)
        delete [] atoms_[i].symbol;

    delete [] atoms_;

}


void TestInfo::TestBasisConversion(void)
{
}

void TestInfo::TestMoleculeConversion(void)
{
    auto mol = MoleculeFromArrays(ncenters_, atoms_);
    molecule_test_.total_ncenters.set_result(mol->natom());
    
    for(int i = 0; i < mol->natom(); i++)
    {
        // at() will throw exception on out-of-bounds
        molecule_test_.xyz[0].at(i).set_result(mol->fx(i));
        molecule_test_.xyz[1].at(i).set_result(mol->fy(i));
        molecule_test_.xyz[2].at(i).set_result(mol->fz(i));
        molecule_test_.Z.at(i).set_result(mol->Z(i));
        molecule_test_.mass.at(i).set_result(mol->mass(i));
    }
}


bool TestInfo::PrintResults(std::ostream & out, bool verbose)
{
    int success = 1;
    int molsuccess = 1;
    int pbasissuccess = 1;
    int abasissuccess = 1;
    out << "--------------------------------------\n";
    out << "Results for test: " << testname_ << "\n";
    out << "--------------------------------------\n";
    out << "Molecule Conversion Test\n\n";

    PrintRow(out, "Name", "Reference", "This run", "Diff", "Threshold", "Pass/Fail");
    PrintSeparator(out);
    molsuccess *= Test(out, "# of centers", molecule_test_.total_ncenters);

    for(int i = 0; i < ncenters_; i++)
    {
        for(int c = 0; c < 3; c++)
        {
            stringstream ss;
            ss << "Center " << i << " (" << ((c == 0) ? 'x' : ((c == 1) ? 'y' : 'z')) << ")";
            molsuccess *= Test(out, ss.str(), molecule_test_.xyz[c].at(i));
        }

        stringstream ssmass;
        ssmass << "Center " << i << " (mass)";
        stringstream ssZ;
        ssZ << "Center " << i << " (Z number)";

        molsuccess *= Test(out, ssmass.str(), molecule_test_.mass.at(i));
        molsuccess *= Test(out, ssZ.str(), molecule_test_.Z.at(i));
    }

    out << "\n***************************************************\n";
    out << "Molecule conversion result: " << (molsuccess ? "PASS" : "FAIL");
    out << "\n***************************************************\n";
    out << "\n\n";

    success *= molsuccess; 







    out << "Basis Set Conversion Test (Primary)\n\n";

    PrintRow(out, "Name", "Reference", "This run", "Diff", "Threshold", "Pass/Fail");
    PrintSeparator(out);
    pbasissuccess *= Test(out, "# of centers", molecule_test_.total_ncenters);

    for(int i = 0; i < ncenters_; i++)
    {
        for(int c = 0; c < 3; c++)
        {
            stringstream ss;
            ss << "Center " << i << " (" << ((c == 0) ? 'x' : ((c == 1) ? 'y' : 'z')) << ")";
            pbasissuccess *= Test(out, ss.str(), molecule_test_.xyz[c].at(i));
        }

        stringstream ssmass;
        ssmass << "Center " << i << " (mass)";
        stringstream ssZ;
        ssZ << "Center " << i << " (Z number)";

        pbasissuccess *= Test(out, ssmass.str(), molecule_test_.mass.at(i));
        pbasissuccess *= Test(out, ssZ.str(), molecule_test_.Z.at(i));
    }

    out << "\n***************************************************\n";
    out << "Basis Set (primary) conversion result: " << (pbasissuccess ? "PASS" : "FAIL");
    out << "\n***************************************************\n";










    out << "\n\n";
    return success;
}


}
} //close namespace panache::testing


