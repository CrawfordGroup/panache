#include <fstream>
#include <sstream>
#include <cassert>
#include <iomanip>

#include "Testing.h"
#include "c_convert.h"
#include "DFTensor.h"
#include "Orderings.h"
#include "BasisFunctionMacros.h"
#include "SlowERIBase.h"
#include "BasisSetParser.h"

#include "ERI.h"

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
                            int_t * &nshellspercenter,
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
        int ncenters, nshells, nbf, nao, nprim;
        f >> nshells >> nbf >> nao >> nprim >> ncenters;

        test.nshell.set(nshells);
        test.nao.set(nao);
        test.nbf.set(nbf);
        test.nprim.set(nprim);

        nshellspercenter = new int_t[ncenters];
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
        int natoms, nallatoms;
        std::string schoen, fullpg; // fullpg is no longer used

        f >> natoms >> nallatoms;
        f.ignore(); // remove newline

        getline(f, schoen);
        getline(f, fullpg);
        test.schoen.set(schoen);

        test.natom.set(natoms);
        test.nallatom.set(nallatoms);

        atoms = new C_AtomCenter[nallatoms];

        double dummy; // for mass and Z-number

        for(int i = 0; i < nallatoms; i++)
        {
            char * s = new char[8];
            f.get(s, 4, ' ');
            atoms[i].symbol = s;

            f >> dummy  // Z is not used
              >> atoms[i].center[0]
              >> atoms[i].center[1]
              >> atoms[i].center[2]
              >> dummy;   // mass is not used

            f.ignore(); // ignore the newline

            test.xyz[0].push_back(DTestResult(atoms[i].center[0], TEST_MOLECULE_XYZ_THRESHOLD));
            test.xyz[1].push_back(DTestResult(atoms[i].center[1], TEST_MOLECULE_XYZ_THRESHOLD));
            test.xyz[2].push_back(DTestResult(atoms[i].center[2], TEST_MOLECULE_XYZ_THRESHOLD));
        }

        return nallatoms;
    }
    catch(...)
    {
        throw TestingParserException("Error parsing molecule file", filename);
    }
}


void TestInfo::ReadMatrixInfo(const string & filename,
                              MatrixTest & test,
                              double element_threshold,
                              double sum_threshold,
                              double checksum_threshold)
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

        test.length.set(nrow*ncol);
        test.sum.set(sum, sum_threshold);
        test.checksum.set(checksum, checksum_threshold);

        double val;
        for(int i = 0; i < 100; i++)
        {
            f >> test.elements[i].index >> val;
            test.elements[i].test.set(val, element_threshold);
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
    string qmo_filename(dir);
    string cmo_filename(dir);

    primary_basis_filename.append("basis.primary");
    aux_basis_filename.append("basis.aux");
    molecule_filename.append("geometry");
    qso_filename.append("qso");
    qmo_filename.append("qmo");
    cmo_filename.append("cmat");

    // Read in the molecule
    ncenters_ = ReadMolecule(molecule_filename, atoms_, molecule_test_);

    // Read in the basis set information
    int_t pcenters = ReadBasisFile(primary_basis_filename, primary_nshellspercenter_, primary_shells_, primary_test_);
    int_t acenters = ReadBasisFile(aux_basis_filename,     aux_nshellspercenter_,     aux_shells_,     aux_test_);


    // Make sure all files agree on the number of centers
    if(pcenters != acenters || pcenters != ncenters_)
        throw TestingSanityException("Error - not all files agree on the number of centers!");


    //! \todo Could be returned from the ReadBasis function (or a pointer passed)
    // We don't really get back the number of shells, but they can be
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
    ReadMatrixInfo(qso_filename, qso_test_, TEST_QSO_ELEMENT_THRESHOLD,
                   TEST_QSO_SUM_THRESHOLD,
                   TEST_QSO_CHECKSUM_THRESHOLD);

    // Read in Qmo tests
    ReadMatrixInfo(qmo_filename, qmo_test_, TEST_QMO_ELEMENT_THRESHOLD,
                   TEST_QMO_SUM_THRESHOLD,
                   TEST_QMO_CHECKSUM_THRESHOLD);

    cmo_matrix_ = ReadCMatrix(cmo_filename);
}

SharedMatrix TestInfo::ReadCMatrix(const string & filename)
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

        f >> nrow >> ncol;

        SharedMatrix mat(new Matrix(nrow, ncol));

        for(int i = 0; i < nrow; i++)
        for(int j = 0; j < ncol; j++)
            f >> mat->pointer(0)[i][j];

        return mat;
    }
    catch(...)
    {
        throw TestingParserException("Error parsing matrix file", filename);
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


void TestInfo::TestBasisConversion(int_t nshells, int_t * nshellspercenter, C_ShellInfo * shells, BasisTest & test)
{
    auto mol = MoleculeFromArrays(ncenters_, atoms_);
    auto basis = BasisSetFromArrays(mol, ncenters_, nshellspercenter, shells, false);

    test.nprim.set_thisrun(basis->nprimitive());
    test.nbf.set_thisrun(basis->nbf());
    test.nao.set_thisrun(basis->nao());
    test.nshell.set_thisrun(basis->nshell());

    // don't keep testing if stuff is wrong
    if(!(test.nprim.check() && test.nbf.check() && test.nao.check() && test.nshell.check()))
       return; 

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
    molecule_test_.natom.set_thisrun(mol->natom());
    molecule_test_.nallatom.set_thisrun(mol->nallatom());

    molecule_test_.schoen.set_thisrun(mol->schoenflies_symbol());

    for(int i = 0; i < mol->natom(); i++)
    {
        // at() will throw exception on out-of-bounds
        molecule_test_.xyz[0].at(i).set_thisrun(mol->fx(i));
        molecule_test_.xyz[1].at(i).set_thisrun(mol->fy(i));
        molecule_test_.xyz[2].at(i).set_thisrun(mol->fz(i));
    }
}

void TestInfo::TestMatrix(double * mat,
                          int size,
                          MatrixTest & test)
{
    double sum = 0;
    double checksum = 0;

    test.length.set_thisrun(size);



    // DON'T TEST ELEMENTS IF DIMENSIONS ARE WRONG
    for(size_t i = 0; i < size; i++)
    {
        sum += mat[i];
        checksum += mat[i]*static_cast<double>(i+1);
    }

    test.sum.set_thisrun(sum);
    test.checksum.set_thisrun(checksum);

    if(test.length.check())
    {
        for(int i = 0; i < 100; i++)
            test.elements[i].test.set_thisrun(mat[test.elements[i].index]);
    }

}


void TestInfo::TestQsoMatrix(void)
{
    auto mol = MoleculeFromArrays(ncenters_, atoms_);
    auto primary = BasisSetFromArrays(mol, ncenters_, primary_nshellspercenter_, primary_shells_, false);
    auto aux = BasisSetFromArrays(mol, ncenters_, aux_nshellspercenter_, aux_shells_, false);

    // to test reading of basis sets from file
    //std::shared_ptr<Gaussian94BasisSetParser> parser(new Gaussian94BasisSetParser);
    //auto primary = BasisSet::construct(parser, mol, "6-31gss.gbs");
    //auto aux = BasisSet::construct(parser, mol, "cc-pvdz-ri.gbs");

    DFTensor dft(primary, aux);
    int naux, nso2;
    size_t matsize = dft.QsoDimensions(naux, nso2);

    dft.GenQso(true);

    double * mat = new double[matsize];

    // Test getting it all at once
    // dft.GetBatch(mat, matsize);

    // test getting in batches
    double * buf = new double[3*nso2];
    int curq = 0;
    int n;
    while(n = dft.GetBatch_Qso(buf, 3*nso2))
    {
        std::copy(buf, buf+n*nso2, mat+curq*nso2);
        curq += n;
    }
    delete [] buf;

    TestMatrix(mat, matsize, qso_test_);
    delete [] mat;
}

void TestInfo::TestQmoMatrix(void)
{
    auto mol = MoleculeFromArrays(ncenters_, atoms_);
    auto primary = BasisSetFromArrays(mol, ncenters_, primary_nshellspercenter_, primary_shells_, false);
    auto aux = BasisSetFromArrays(mol, ncenters_, aux_nshellspercenter_, aux_shells_, false);

    DFTensor dft(primary, aux);
    int naux, nso2;
    size_t matsize = dft.QsoDimensions(naux, nso2);

    int nso = primary->nbf();
    dft.SetCMatrix(cmo_matrix_->pointer(0)[0], nso, false);
    dft.GenQso(true);

    double * mat = new double[matsize];


    // Test getting it all at once
    // dft.GetBatch(mat, matsize);

    // test getting in batches
    double * buf = new double[3*nso2];
    int curq = 0;
    int n;
    while(n = dft.GetBatch_Qmo(buf, 3*nso2))
    {
        std::copy(buf, buf+n*nso2, mat+curq*nso2);
        curq += n;
    }
    delete [] buf;

    TestMatrix(mat, matsize, qmo_test_);
    delete [] mat;
}

/*
void TestInfo::TestERI(void)
{
    std::cout << "Testing ERI generation from Qso\n";

    auto mol = MoleculeFromArrays(ncenters_, atoms_);
    auto primary = BasisSetFromArrays(mol, ncenters_, primary_nshellspercenter_, primary_shells_, false);
    auto aux = BasisSetFromArrays(mol, ncenters_, aux_nshellspercenter_, aux_shells_, false);

    DFTensor dft(primary, aux);

    int nbf = primary->nbf();
    int nbf2 = nbf*nbf;
    int naux = aux->nbf();

    size_t matsize = naux * nbf2;
    double * qso = new double[matsize];
    dft.Qso(qso, matsize);


    std::shared_ptr<BasisSet> zero = BasisSet::zero_ao_basis_set();

    std::shared_ptr<TwoBodyAOInt> referi = GetERI(primary, primary, primary, primary);

    std::cout << std::setprecision(8);

    double *dfbuf = new double[10000]; // up to (ff|ff)

    for(int p = 0; p < primary->nshell(); p++)
    {
        for(int q = 0; q < primary->nshell(); q++)
        {
            for(int r = 0; r < primary->nshell(); r++)
            {
                for(int s = 0; s < primary->nshell(); s++)
                {
                    //Libint
                    int nreferi = referi->compute_shell(p,q,r,s);

                    //Density fitting
                    int ndferi = dft.CalculateERI(qso, matsize, p, q, r, s, dfbuf, 10000); 

                    std::cout << "Quartet " << p << " " << q << " " << r << " " << s << "\n";
                    for(int i = 0; i < ndferi; i++)
                        std::cout << "  Reference: " << referi->buffer()[i] << "\n";
                    for(int i = 0; i < ndferi; i++)
                        std::cout << "Density Fit: " << dfbuf[i] << "\n";
                    std::cout << "\n";

                }
            }
        }
    }

    delete [] dfbuf;
    delete [] qso;
}
*/

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
    failures +=  Test(out, "# of ao", test.nao);
    failures +=  Test(out, "# of bf", test.nbf);


    if(failures > 0)
    {
        out << "\n...skipping remainder of tests since dimensions didn't check out\n\n";
        return failures;
    }

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
    failures +=  Test(out, "# of elements", test.length);

    bool testelements = (failures > 0) ? false : true;
        
    failures +=  Test(out, "sum", test.sum);
    failures +=  Test(out, "checksum", test.checksum);

    if(!testelements)
    {
        out << "\n...skipping element testing since dimensions didn't check out\n\n";
        return failures;
    }

    for(int i = 0; i < 100; i++)
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
    mol_failures += Test(out, "# of atoms", molecule_test_.natom);
    mol_failures += Test(out, "# of atoms (+ dummies)", molecule_test_.nallatom);
    mol_failures += Test(out, "Schoenflies Symbol", molecule_test_.schoen);

    for(int i = 0; i < ncenters_; i++)
    {
        for(int c = 0; c < 3; c++)
        {
            stringstream ss;
            ss << "Center " << i << " (" << ((c == 0) ? 'x' : ((c == 1) ? 'y' : 'z')) << ")";
            mol_failures += Test(out, ss.str(), molecule_test_.xyz[c].at(i));
        }
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

    out << "\n\n\n";
    failures += PrintMatrixResults(out, "QMO", qmo_test_);

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


}
} //close namespace panache::testing

