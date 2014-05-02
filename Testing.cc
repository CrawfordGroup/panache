#include <fstream>
#include <sstream>

#include "Testing.h"
#include "c_convert.h"
#include "DFTensor.h"

#ifdef PANACHE_DEVELOPER_GENERATE
#include <utility>
#include <algorithm>
#include <array>
#include <iomanip>
using std::pair;
using std::array;
#endif

using std::string;
using std::ifstream;
using std::stringstream;


namespace panache
{
namespace testing
{


// Note - exceptions are turned on for the ifstream object
// so that any parsing errors just throw an exeption. Catch those,
// and throw a TestingParserException
int TestInfo::ReadBasisFile(const string & filename,
                            int * &nshellspercenter,
                            C_ShellInfo * &shells,
                            TestInfo::BasisTest & test)
{
    ifstream f(filename.c_str());

    if(!f.is_open())
        throw TestingParserException("Cannot open basis file!", filename);

    f.exceptions(std::ifstream::failbit |
                 std::ifstream::badbit  |
                 std::ifstream::eofbit);

    try
    {
        int ncenters, nshells;
        f >> nshells >> ncenters;

        test.nshell.set(nshells);

        nshellspercenter = new int[ncenters];
        for(int i = 0; i < ncenters; i++)
        {
            f >> nshellspercenter[i];
            test.center_nshell.push_back(ITestResult(nshellspercenter[i]));
        }


        shells = new C_ShellInfo[nshells];

        int dum; // holds the shell index. Not needed
        int totalnprim = 0;

        for(int i = 0; i < nshells; i++)
        {
            f >> dum >> shells[i].nprim >> shells[i].am >> shells[i].ispure;
            test.shell_nprim.push_back(ITestResult(shells[i].nprim));
            test.shell_am.push_back(ITestResult(shells[i].am));
            test.shell_ispure.push_back(ITestResult(shells[i].ispure));

            shells[i].coef = new double[shells[i].nprim];
            shells[i].exp = new double[shells[i].nprim];

            for(int j = 0; j < shells[i].nprim; j++)
            {
                f >> shells[i].coef[j] >> shells[i].exp[j];
                test.exp.push_back(DTestResult(shells[i].exp[j], TEST_EXPONENT_THRESHOLD));
                test.coef.push_back(DTestResult(shells[i].coef[j], TEST_COEFFICIENT_THRESHOLD));
                totalnprim++;
            }

        }

        test.nprim.set(totalnprim);

        return ncenters;
    }
    catch(...)
    {
        throw TestingParserException("Error parsing basis file", filename);
    }
}

int TestInfo::ReadMolecule(const string & filename, C_AtomCenter * &atoms, MoleculeTest & test)
{
    ifstream f(filename.c_str());

    if(!f.is_open())
        throw TestingParserException("Cannot open file!", filename);

    f.exceptions(std::ifstream::failbit |
                 std::ifstream::badbit  |
                 std::ifstream::eofbit);
    try
    {

        int ncenters;
        f >> ncenters;

        test.ncenters.set(ncenters);

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

            test.xyz[0].push_back(DTestResult(atoms[i].center[0], TEST_MOLECULE_XYZ_THRESHOLD));
            test.xyz[1].push_back(DTestResult(atoms[i].center[1], TEST_MOLECULE_XYZ_THRESHOLD));
            test.xyz[2].push_back(DTestResult(atoms[i].center[2], TEST_MOLECULE_XYZ_THRESHOLD));
            test.Z.push_back(DTestResult(atoms[i].Z, TEST_MOLECULE_Z_THRESHOLD));
            test.mass.push_back(DTestResult(atoms[i].mass, TEST_MOLECULE_MASS_THRESHOLD));
        }

        return ncenters;
    }
    catch(...)
    {
        throw TestingParserException("Error parsing molecule file", filename);
    }
}


void TestInfo::ReadMatrixInfo(const string & filename,
                              MatrixTest & test,
                              double threshold)
{
    ifstream f(filename.c_str());

    if(!f.is_open())
        throw TestingParserException("Cannot open file!", filename);

    f.exceptions(std::ifstream::failbit |
                 std::ifstream::badbit  |
                 std::ifstream::eofbit);
    try
    {
        int nrow, ncol;
        double sum, checksum;

        f >> nrow >> ncol >> sum >> checksum;

        test.nrow.set(nrow);
        test.ncol.set(ncol);
        test.sum.set(sum);
        test.checksum.set(checksum);

        double val;
        for(int i = 0; i < 50; i++)
        {
            f >> test.elements[i].index >> val;
            test.elements[i].test.set(val, threshold);
        }
    }
    catch(...)
    {
        throw TestingParserException("Error parsing matrix file", filename);
    }
}


TestInfo::TestInfo(const std::string & testname, const std::string & dir)
    : testname_(testname)
{
    string primary_basis_filename(dir);
    string aux_basis_filename(dir);
    string molecule_filename(dir);
    string qso_filename(dir);

    primary_basis_filename.append("basis.primary");
    aux_basis_filename.append("basis.aux");
    molecule_filename.append("geometry");
    qso_filename.append("qso");

    // Read in the molecule
    ncenters_ = ReadMolecule(molecule_filename, atoms_, molecule_test_);

    // Read in the basis set information
    int pcenters = ReadBasisFile(primary_basis_filename, primary_nshellspercenter_, primary_shells_, primary_test_);
    int acenters = ReadBasisFile(aux_basis_filename,     aux_nshellspercenter_,     aux_shells_,     aux_test_);

    


    // Make sure all files agree on the number of centers
    if(pcenters != acenters || pcenters != ncenters_)
        throw TestingSanityException("Error - not all files agree on the number of centers!");


    // We don't really read in the number of shells, but they can be
    // deduced from the number of centers and number of shells per center
    primary_nshells_ = aux_nshells_ = 0;
    primary_nprim_ = aux_nprim_ = 0;

    for(int i = 0; i < ncenters_; i++)
    {
        primary_nshells_ += primary_nshellspercenter_[i];
        aux_nshells_ += aux_nshellspercenter_[i];
    }

    // Same as above, but we can get the number of primitives
    for(int i = 0; i < primary_nshells_; i++)
        primary_nprim_ += primary_shells_[i].nprim;

    for(int i = 0; i < aux_nshells_; i++)
        aux_nprim_ += aux_shells_[i].nprim;


    // Read in Qso tests
    ReadMatrixInfo(qso_filename, qso_test_, TEST_QSO_ELEMENT_THRESHOLD);
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


void TestInfo::TestBasisConversion(int nshells, int * nshellspercenter, C_ShellInfo * shells, BasisTest & test)
{
    auto mol = MoleculeFromArrays(ncenters_, atoms_);
    auto basis = BasisSetFromArrays(mol, ncenters_, nshellspercenter, shells);

    test.nprim.set_thisrun(basis->nprimitive());
    test.nshell.set_thisrun(basis->nshell());

    for(int i = 0; i < ncenters_; i++)
        test.center_nshell.at(i).set_thisrun(basis->nshell_on_center(i));

    int counter = 0;
    for(int i = 0; i < nshells; i++)
    {
        test.shell_nprim.at(i).set_thisrun(basis->shell(i).nprimitive());
        test.shell_am.at(i).set_thisrun(basis->shell(i).am());
        test.shell_ispure.at(i).set_thisrun(basis->shell(i).is_pure());

        for(int j = 0; j < basis->shell(i).nprimitive(); j++)
        {
            test.exp.at(counter).set_thisrun(basis->shell(i).exp(j));
            test.coef.at(counter).set_thisrun(basis->shell(i).original_coef(j));
            counter++;
        }
    }
}


void TestInfo::TestBasisConversion(void)
{
    TestBasisConversion(primary_nshells_, primary_nshellspercenter_, primary_shells_, primary_test_);
    TestBasisConversion(aux_nshells_,     aux_nshellspercenter_,     aux_shells_,     aux_test_);
}


void TestInfo::TestMoleculeConversion(void)
{
    auto mol = MoleculeFromArrays(ncenters_, atoms_);
    molecule_test_.ncenters.set_thisrun(mol->natom());

    for(int i = 0; i < mol->natom(); i++)
    {
        // at() will throw exception on out-of-bounds
        molecule_test_.xyz[0].at(i).set_thisrun(mol->fx(i));
        molecule_test_.xyz[1].at(i).set_thisrun(mol->fy(i));
        molecule_test_.xyz[2].at(i).set_thisrun(mol->fz(i));
        molecule_test_.Z.at(i).set_thisrun(mol->Z(i));
        molecule_test_.mass.at(i).set_thisrun(mol->mass(i));
    }
}

void TestInfo::TestMatrix(SharedMatrix mat, 
                          MatrixTest & test)
{
    //! \todo matrices with > 1 irrep?
    double * p = &(mat->pointer(0)[0][0]);
    size_t size = mat->colspi()[0] * mat->rowspi()[0];

    test.nrow.set_thisrun(mat->rowspi()[0]);
    test.ncol.set_thisrun(mat->colspi()[0]);

    double sum = 0;
    double checksum = 0;

    for(size_t i = 0; i < size; i++)
    {
        sum += p[i];
        checksum += p[i]*static_cast<double>(i);
    }

    test.sum.set_thisrun(sum);
    test.checksum.set_thisrun(checksum);

    // DON'T TEST ELEMENTS IF DIMENSIONS ARE WRONG
    if(test.nrow.check() && test.ncol.check())
    {
        for(int i = 0; i < 50; i++)
            test.elements[i].test.set_thisrun(p[test.elements[i].index]);
    }
} 


void TestInfo::TestQsoMatrix(void)
{
    auto mol = MoleculeFromArrays(ncenters_, atoms_);
    auto primary = BasisSetFromArrays(mol, ncenters_, primary_nshellspercenter_, primary_shells_);
    auto aux = BasisSetFromArrays(mol, ncenters_, aux_nshellspercenter_, aux_shells_);

    DFTensor dft(primary, aux);
    auto mat = dft.Qso();
    TestMatrix(mat, qso_test_);
}

int TestInfo::PrintBasisResults(std::ostream & out, const string & type,
                                int nshells, int nprim,
                                const BasisTest & test,
                                bool verbose)
{
    int failures = 0;

    out << "Basis Set Conversion Test (" << type << ")\n\n";

    PrintRow(out, "Name", "Reference", "This run", "Diff", "Threshold", "Pass/Fail");
    PrintSeparator(out);
    failures +=  Test(out, "# of primitives", test.nprim);
    failures +=  Test(out, "# of shells", test.nshell);

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


int TestInfo::PrintMatrixResults(std::ostream & out, const string & type,
                                 const MatrixTest & test, bool verbose)
{
    int failures = 0;

    out << "Matrix test: " << type << "\n\n";

    PrintRow(out, "Name", "Reference", "This run", "Diff", "Threshold", "Pass/Fail");
    PrintSeparator(out);
    failures +=  Test(out, "# of rows", test.nrow);
    failures +=  Test(out, "# of columns", test.ncol);
    failures +=  Test(out, "sum", test.sum);
    failures +=  Test(out, "checksum", test.checksum);

    for(int i = 0; i < 50; i++)
    {
        stringstream ss;
        ss << "Element " << test.elements[i].index;
        failures +=  Test(out, ss.str(), test.elements[i].test);
    }

    out << "\n";
    out << "***********************************************************************\n";
    out << "Matrix \"" << type << "\" result: " << (failures ? "FAIL" : "PASS");

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
    mol_failures += Test(out, "# of centers", molecule_test_.ncenters);

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



    ////////////////////////////////////////////////////////////////
    // Call the separate function for the two different basis sets
    ////////////////////////////////////////////////////////////////
    out << "\n\n\n";

    failures += PrintBasisResults(out, "Primary", primary_nshells_, primary_nprim_,
                                  primary_test_, verbose);
    out << "\n\n\n";

    failures += PrintBasisResults(out, "Auxiliary", aux_nshells_, aux_nprim_,
                                  aux_test_, verbose);


    ////////////////////////////////////////////////////////////////
    // Matrix testing
    ////////////////////////////////////////////////////////////////
    out << "\n\n\n";
    failures += PrintMatrixResults(out, "QSO", qso_test_);

    ////////////////////////////////////////////////////////////////
    // Overall results
    ////////////////////////////////////////////////////////////////
    out << "\n\n\n";
    out << "=============================================================\n";
    out << "= OVERALL RESULT: " << (failures ? "FAIL" : "PASS");

    if(failures)
        out << " (" << failures << " failures)";

    out << "\n";
    out << "=============================================================\n";
    out << "\n";

    return failures;
}


#ifdef PANACHE_DEVELOPER_GENERATE
void TestInfo::Generate(void)
{
    auto mol = MoleculeFromArrays(ncenters_, atoms_);
    auto primary = BasisSetFromArrays(mol, ncenters_, primary_nshellspercenter_, primary_shells_);
    auto aux = BasisSetFromArrays(mol, ncenters_, aux_nshellspercenter_, aux_shells_);


    DFTensor dft(primary, aux);
    auto mat = dft.Qso();

    double * p = &(mat->pointer(0)[0][0]);
    size_t size = mat->colspi()[0] * mat->rowspi()[0];

    double sum = 0;
    double checksum = 0;
    array<pair<size_t, double>, 50> largest;


    for(size_t i = 0; i < size; i++)
    {
        sum += p[i];
        checksum += p[i]*static_cast<double>(i);

        if(p[i] > largest[0].second)
        {
            largest[0].first = i;
            largest[0].second = p[i];
            std::sort(largest.begin(), largest.end(), [](pair<size_t, double> p1, pair<size_t, double> p2)
            {
                return p1.second < p2.second;
            });
        }
    }

    std::cout << "---------begin QSO data---------\n";
    std::cout << mat->rowspi(0) << "\n" << mat->colspi(0) << "\n";

    // to return the precision back to normal
    auto prec = std::cout.precision();
    std::cout << std::setprecision(20) << sum << "\n";
    std::cout << checksum << "\n";
    for(auto & it : largest)
        std::cout << it.first << it.second << "\n";
    std::cout.precision(prec);

}
#endif


}
} //close namespace panache::testing


