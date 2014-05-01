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
        test.shell_am.push_back(TestResult<int>(shells[i].am));
        test.shell_ispure.push_back(TestResult<int>(shells[i].ispure));

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

    primary_nshells_ = aux_nshells_ = 0;
    primary_nprim_ = aux_nprim_ = 0;

    for(int i = 0; i < ncenters_; i++)
    {
        primary_nshells_ += primary_nshellspercenter_[i];
        aux_nshells_ += aux_nshellspercenter_[i];
    }

    for(int i = 0; i < primary_nshells_; i++)
        primary_nprim_ += primary_shells_[i].nprim;

    for(int i = 0; i < aux_nshells_; i++)
        aux_nprim_ += aux_shells_[i].nprim;

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


void TestInfo::TestBasisConversion(int nshells,
                                   int * nshellspercenter,
                                   C_ShellInfo * shells,
                                   BasisTest & test)
{
    auto mol = MoleculeFromArrays(ncenters_, atoms_);
    auto basis = BasisSetFromArrays(mol,
                                    ncenters_,
                                    nshellspercenter,
                                    shells);

    test.total_nprim.set_result(basis->nprimitive());
    test.total_nshell.set_result(basis->nshell());

    for(int i = 0; i < ncenters_; i++)
        test.center_nshell.at(i).set_result(basis->nshell_on_center(i));

    int counter = 0;
    for(int i = 0; i < nshells; i++)
    {
        test.shell_nprim.at(i).set_result(basis->shell(i).nprimitive());
        test.shell_am.at(i).set_result(basis->shell(i).am());
        test.shell_ispure.at(i).set_result(basis->shell(i).is_pure());

        for(int j = 0; j < basis->shell(i).nprimitive(); j++)
        {
            test.exp.at(counter).set_result(basis->shell(i).exp(j));
            test.coef.at(counter).set_result(basis->shell(i).original_coef(j));
            counter++;
        }
    }
}

void TestInfo::TestBasisConversion(void)
{
    TestBasisConversion(primary_nshells_,
                        primary_nshellspercenter_,
                        primary_shells_,
                        primary_test_);    
    TestBasisConversion(aux_nshells_,
                        aux_nshellspercenter_,
                        aux_shells_,
                        aux_test_);    

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


int TestInfo::PrintBasisResults(std::ostream & out, const string & type,
                                 int nshells, int nprim,
                                 BasisTest & test,
                                 bool verbose)
{
    int failures = 0;

    out << "Basis Set Conversion Test (" << type << ")\n\n";

    PrintRow(out, "Name", "Reference", "This run", "Diff", "Threshold", "Pass/Fail");
    PrintSeparator(out);
    failures +=  Test(out, "# of primitives", test.total_nprim);
    failures +=  Test(out, "# of shells", test.total_nshell);

    for(int i = 0; i < ncenters_; i++)
    {
        stringstream ss;
        ss << "# of shells / center " << i;
        failures +=  Test(out, ss.str(), test.center_nshell.at(i));
    }


    for(int i = 0; i < nshells; i++)
    {
        stringstream ss;
        ss << "# of primitives / shell " << i;
        failures +=  Test(out, ss.str(), test.shell_nprim.at(i));
    }
    for(int i = 0; i < nshells; i++)
    {
        stringstream ss;
        ss << "Angular momentum / shell " << i;
        failures +=  Test(out, ss.str(), test.shell_am.at(i));
    }
    for(int i = 0; i < nshells; i++)
    {
        stringstream ss;
        ss << "Is pure? / shell " << i;
        failures +=  Test(out, ss.str(), test.shell_ispure.at(i));
    }


    for(int i = 0; i < nprim; i++)
    {
        stringstream ss;
        ss << "Exponent (Primitive " << i << ")";
        failures +=  Test(out, ss.str(), test.exp.at(i));
    } 

    for(int i = 0; i < nprim; i++)
    {
        stringstream ss;
        ss << "Coefficient (Primitive " << i << ")";
        failures +=  Test(out, ss.str(), test.coef.at(i));
    } 


    out << "\n";
    out << "***********************************************************************\n";
    out << "Basis Set (" << type << ") conversion result: " << (failures ? "FAIL" : "PASS");

    if(failures)
        out << " (" << failures << " failures)";

    out << "\n";
    out << "***********************************************************************\n";
   
    return failures; 
}


int TestInfo::PrintResults(std::ostream & out, bool verbose)
{
    int failures = 0;
    int mol_failures = 0;

    out << "--------------------------------------\n";
    out << "Results for test: " << testname_ << "\n";
    out << "--------------------------------------\n";
    out << "Molecule Conversion Test\n\n";

    PrintRow(out, "Name", "Reference", "This run", "Diff", "Threshold", "Pass/Fail");
    PrintSeparator(out);
    mol_failures += Test(out, "# of centers", molecule_test_.total_ncenters);

    for(int i = 0; i < ncenters_; i++)
    {
        for(int c = 0; c < 3; c++)
        {
            stringstream ss;
            ss << "Center " << i << " (" << ((c == 0) ? 'x' : ((c == 1) ? 'y' : 'z')) << ")";
            mol_failures += Test(out, ss.str(), molecule_test_.xyz[c].at(i));
        }

        stringstream ssmass;
        ssmass << "Center " << i << " (mass)";
        stringstream ssZ;
        ssZ << "Center " << i << " (Z number)";

        mol_failures += Test(out, ssmass.str(), molecule_test_.mass.at(i));
        mol_failures += Test(out, ssZ.str(), molecule_test_.Z.at(i));
    }

    out << "\n";
    out << "***********************************************************************\n";
    out << "Molecule conversion result: " << (mol_failures ? "FAIL" : "PASS");

    if(mol_failures)
        out << " (" << mol_failures << " failures)";

    out << "\n";
    out << "***********************************************************************\n";
    failures += mol_failures; 

    out << "\n\n\n";
    
    failures += PrintBasisResults(out, "Primary", primary_nshells_, primary_nprim_,
                                  primary_test_, verbose);
    out << "\n\n\n";

    failures += PrintBasisResults(out, "Auxiliary", aux_nshells_, aux_nprim_,
                                  aux_test_, verbose);


    out << "\n\n";
    out << "=============================================================\n";
    out << "= OVERALL RESULT: " << (failures ? "FAIL" : "PASS");

    if(failures)
        out << " (" << failures << " failures)";

    out << "\n";
    out << "=============================================================\n";
    out << "\n"; 

    return failures;
}


}
} //close namespace panache::testing


