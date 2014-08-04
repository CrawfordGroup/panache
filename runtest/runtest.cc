#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <array>
#include <memory>
#include <utility>
#include <cmath>

#include "panache/DFTensor.h"
#include "panache/DFTensor2.h"
#include "panache/SimpleMatrix.h"
#include "panache/Output.h"
#include "panache/c_convert.h" // int_t comes in through here

#define QSO_ELEMENT_THRESHOLD 1e-11
#define QSO_SUM_THRESHOLD 1e-8
#define QSO_CHECKSUM_THRESHOLD 1.0

#define QMO_ELEMENT_THRESHOLD 1e-10
#define QMO_SUM_THRESHOLD 1e-6
#define QMO_CHECKSUM_THRESHOLD 2.0

using namespace panache;
using namespace std;

void PrintUsage(void)
{
    cout << "\n"
         << "Libpanache Testing Utility\n"
         << "\n"
         << "Usage: runtest [opt] <dir>\n"
         << "\n"
         << "Options:\n"
         << "-v           Verbose printing\n"
         << "-d           Write Qso to disk (rather than in core)\n"
         << "-b           Get Qso/Qmo in batches\n"
         << "-t           Use transpose of C matrix\n"
         << "-h           Print help (you're looking at it\n"
         << "<dir>        Directory holding the test information\n"
         << "\n\n";
}


struct TestMatrixInfo
{
    int nrow, ncol;
    double sum;
    double checksum;

    // index is first element of pair, value is second
    array<pair<int, double>, 100> elements;
};




// Note - exceptions are turned on for the ifstream object
// so that any parsing errors just throw an exeption. Catch those,
// and throw an exception
shared_ptr<BasisSet> ReadBasisFile(shared_ptr<Molecule> mol, const string & filename)
{
    ifstream f(filename.c_str());

    if(!f.is_open())
        throw std::runtime_error("Cannot open basis file!");

    f.exceptions(std::ifstream::failbit |
                 std::ifstream::badbit  |
                 std::ifstream::eofbit);

    try
    {
        int_t ncenters, nshells, nbf, nao, nprim;
        f >> nshells >> nbf >> nao >> nprim >> ncenters;

        unique_ptr<int_t[]> nshellspercenter(new int_t[ncenters]);
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
        // Last argument - these are not normalized
        auto ret = BasisSetFromArrays(mol, ncenters, nshellspercenter.get(), shells.get(), false);

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




shared_ptr<Molecule> ReadMoleculeFile(const string & filename)
{
    ifstream f(filename.c_str());

    if(!f.is_open())
        throw runtime_error("Cannot open file!");

    f.exceptions(std::ifstream::failbit |
                 std::ifstream::badbit  |
                 std::ifstream::eofbit);
    try
    {
        int_t natoms, nallatoms;
        std::string schoen, fullpg; // fullpg is no longer used

        f >> natoms >> nallatoms;
        f.ignore(); // remove newline

        getline(f, schoen);
        getline(f, fullpg);

        unique_ptr<C_AtomCenter[]> atoms(new C_AtomCenter[nallatoms]);

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
        }

        auto ret = MoleculeFromArrays(nallatoms, atoms.get());

        for(int i = 0; i < nallatoms; i++)
            delete [] atoms[i].symbol;

        return ret;
    }
    catch(...)
    {
        throw std::runtime_error("Error parsing molecule file");
    }
}



/*! \brief Prints a row in the results table (a single test)
 *
 * \param [in] out Stream to print to
 * \param [in] name Descriptive name of the test
 * \param [in] ref  The reference (expected) value
 * \param [in] run  The value for this run
 * \param [in] diff The absolute difference
 * \param [in] thresh Threshold for failing the test
 * \param [in] result String of the result (usually "pass" or "FAIL")
 */
template<typename T>
static void PrintRow(const std::string name,
                     const T & ref,
                     const T & run,
                     const T & diff,
                     const T & thresh,
                     const std::string & result)
{
    cout << std::setw(35) << std::left << name;
    cout << std::setw(15) << std::left << ref;
    cout << std::setw(15) << std::left << run;
    cout << std::setw(15) << std::left << diff;
    cout << std::setw(15) << std::left << thresh;
    cout << std::setw(15) << std::left << result << "\n";
}




/*! \brief Overload for printing strings (usually header rows)
 *
 * Used to prevent ambiguities when compiling
 *
 * \param [in] name Descriptive name of the test
 * \param [in] ref  The reference (expected) value
 * \param [in] run  The value for this run
 * \param [in] diff The absolute difference
 * \param [in] thresh Threshold for failing the test
 * \param [in] result String of the result (usually "pass" or "FAIL")
 */
static void PrintRow(const char * name,
                     const char * ref,
                     const char * run,
                     const char * diff,
                     const char * thresh,
                     const char * result)
{
    cout << std::setw(35) << std::left << name;
    cout << std::setw(15) << std::left << ref;
    cout << std::setw(15) << std::left << run;
    cout << std::setw(15) << std::left << diff;
    cout << std::setw(15) << std::left << thresh;
    cout << std::setw(15) << std::left << result << "\n";
}


/*! \brief Prints a line of 104 dashes
 *
 * \param [in] out Stream to print to
 */
static void PrintSeparator(void)
{
    cout << std::string(104, '-') << std::endl;
}



TestMatrixInfo ReadMatrixInfo(const string & filename)
{
    ifstream f(filename.c_str());

    if(!f.is_open())
        throw runtime_error("Cannot open matrix file!");

    f.exceptions(std::ifstream::failbit |
                 std::ifstream::badbit  |
                 std::ifstream::eofbit);
    try
    {
        TestMatrixInfo tmi;

        f >> tmi.nrow >> tmi.ncol >> tmi.sum >> tmi.checksum;

        for(int i = 0; i < 100; i++)
            f >> tmi.elements[i].first >> tmi.elements[i].second;

        return tmi;
    }
    catch(...)
    {
        throw runtime_error("Error parsing matrix file");
    }
}


template<typename T> int TestAndPrint(const string & title, T val, T ref, T threshold, bool alwaysprint)
{
    bool ok = true;

    T diff = std::abs(val - ref);
    if(diff > threshold)
        ok = false;

    if(alwaysprint || !ok)
        PrintRow(title, ref, val, diff, threshold, (ok ? "pass" : "FAIL"));

    // return 1 if it fails
    return (ok ? 0 : 1);
}

int TestAndPrint(const string & title, const string & val, string & ref, bool alwaysprint)
{
    bool ok = true;

    if(val != ref)
        ok = false;

    if(alwaysprint || !ok)
        PrintRow(title, ref, val, string(ok ? "-" : "Y"), string("-"), (ok ? "pass" : "FAIL"));

    // return 1 if it fails
    return (ok ? 0 : 1);
}



int TestMatrix(const string & title, const string & reffile,
               double * mat, int matsize,
               double sum_threshold, double checksum_threshold, double element_threshold,
               bool verbose)
{
    TestMatrixInfo tmi = ReadMatrixInfo(reffile);
    int nfailures = 0;

    PrintRow("Name", "Reference", "This run", "Diff", "Threshold", "Pass/Fail");
    PrintSeparator();


    nfailures += TestAndPrint("# of elements", matsize, tmi.nrow*tmi.ncol, 0, verbose);

    // DON'T TEST ELEMENTS IF DIMENSIONS ARE WRONG
    if(nfailures > 0)
    {
        cout << "... not continuing with matrix test since dimensions are wrong!\n";
        return nfailures;
    }

    double sum = 0;
    double checksum = 0;

    for(int i = 0; i < matsize; i++)
    {
        sum += mat[i];
        checksum += mat[i]*static_cast<double>(i+1);
    }

    nfailures += TestAndPrint("sum", sum, tmi.sum, sum_threshold, verbose);
    nfailures += TestAndPrint("checksums", checksum, tmi.checksum, checksum_threshold, verbose);

    for(int i = 0; i < 100; i++)
    {
        stringstream ss;
        ss << "Element " << tmi.elements[i].first;

        nfailures += TestAndPrint(ss.str(), mat[tmi.elements[i].first], tmi.elements[i].second,
                                  element_threshold, verbose);
    }

    if(!verbose && nfailures == 0)
        cout << "           ... no failed tests!\n";
    cout << "\n";
    cout << "***********************************************************************\n";
    cout << "Matrix \"" << title << "\" result: " << (nfailures ? "FAIL" : "PASS");

    if(nfailures)
        cout << " (" << nfailures << " failures)";

    cout << "\n";
    cout << "***********************************************************************\n";

    return nfailures;
}


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





int main(int argc, char ** argv)
{
    int ret = 0;

    try
    {

        if(argc == 1)
        {
            PrintUsage();
            return 0;
        }

        string dir;

        bool verbose = false;
        bool inmem = true;
        bool transpose = false;
        int batchsize = 0;

        for(int i = 1; i < argc; i++)
        {
            string starg(argv[i]);
            if(starg == "-h")
            {
                PrintUsage();
                return 0;
            }
            else if(starg == "-b")
                batchsize = 3;
            else if(starg == "-d")
                inmem = false;
            else if(starg == "-v")
                verbose = true;
            else if(starg == "-t")
                transpose = true;
            else
            {
                // add trailing slash if needed
                if(starg[starg.length()-1] != '/')
                    starg.append("/");
                dir = starg;
            }
        }

        if(verbose)
            panache::output::SetOutput(&cout);

        // Get the test description
        string desc_path(dir);
        desc_path.append("desc");

        fstream desc(desc_path.c_str());
        if(!desc.is_open())
        {
            cout << "\nFatal Error: Cannot open " << desc_path << "\"\n\n";
            return -1;
        }

        string desc_str;
        getline(desc, desc_str);

        cout << "--------------------------------------\n";
        cout << "Results for test: " << desc_str << "\n";
        cout << "--------------------------------------\n";

        string primary_basis_filename(dir);
        string aux_basis_filename(dir);
        string molecule_filename(dir);
        string nocc_filename(dir);
        string qso_filename(dir);
        string qov_filename(dir);
        string qoo_filename(dir);
        string qvv_filename(dir);
        string qmo_filename(dir);
        string cmo_filename(dir);

        primary_basis_filename.append("basis.primary");
        aux_basis_filename.append("basis.aux");
        molecule_filename.append("geometry");
        nocc_filename.append("nocc");
        qso_filename.append("qso");
        qmo_filename.append("qmo");
        qov_filename.append("qov");
        qoo_filename.append("qoo");
        qvv_filename.append("qvv");
        cmo_filename.append("cmat");

        auto mol = ReadMoleculeFile(molecule_filename);
        auto primary = ReadBasisFile(mol, primary_basis_filename);
        auto aux = ReadBasisFile(mol, aux_basis_filename);
        auto cmat = ReadCMatrixFile(cmo_filename);

        if(transpose)
        {
            std::shared_ptr<SimpleMatrix> cmatt(new SimpleMatrix(cmat->ncol(), cmat->nrow()));
            for(size_t i = 0; i < cmat->nrow(); i++)
            for(size_t j = 0; j < cmat->ncol(); j++)
                (*cmatt)(j,i) = (*cmat)(i,j);

            std::swap(cmatt, cmat);
            // original cmat (now in cmatt), will be destructed here
        }

        int nso = primary->nbf();
        int naux = aux->nbf();
        int nso2 = nso*nso;
        int nsotri = nso*(nso+1)/2;
        int nocc = ReadNocc(nocc_filename);
        int nmo = nso;
        int nmo2 = nmo*nmo;

        DFTensor dft(primary, aux, "Test.mat", 0);
        dft.SetCMatrix(cmat->pointer(), nmo, transpose);
        dft.SetNOcc(nocc); 

        DFTensor2 dft2(primary, aux, "/tmp", 0);
        dft2.SetCMatrix(cmat->pointer(), nmo, transpose);
        dft2.SetNOcc(nocc);

        if(inmem)
            dft2.GenQTensors(15, DFTensor2::QStorage::INMEM);
        else
            dft2.GenQTensors(15, DFTensor2::QStorage::ONDISK);


        size_t matsize = dft.QsoDimensions(naux, nso2);
        size_t mattmpsize = nsotri * naux;

        size_t buffsize = matsize;

        if(batchsize)
            buffsize = batchsize * nso2;


        dft.GenQso(inmem);
        int curq = 0;
        int n;

        ///////////
        // Test Qso
        ///////////
        unique_ptr<double[]> mat(new double[matsize]);
        unique_ptr<double[]> mattmp(new double[mattmpsize]);
        unique_ptr<double[]> outbuf(new double[buffsize]);
        dft.SetOutputBuffer(outbuf.get(), buffsize);

        while((n = dft.GetBatch_Qso()))
        {
            std::copy(outbuf.get(), outbuf.get() + n*nso2, mat.get() + curq*nso2);
            curq += n;
        }

        ret += TestMatrix("QSO", qso_filename,
                          mat.get(), matsize,
                          QSO_SUM_THRESHOLD, QSO_CHECKSUM_THRESHOLD, QSO_ELEMENT_THRESHOLD,
                          verbose);

        dft.ResetBatches();


        std::fill(mat.get(), mat.get()+matsize, 99);
        dft2.GetQBatch_Qso(mattmp.get(), mattmpsize, 0);

        // expand
        for(int q = 0; q < naux; q++)
        for(int i = 0; i < nso; i++)
        for(int j = 0; j <= i; j++)
            mat[q*nso2+i*nso+j] = mat[q*nso2+j*nso+i] = mattmp[q*nsotri+i*(i+1)/2+j];

        ret += TestMatrix("QSO2", qso_filename,
                          mat.get(), matsize,
                          QSO_SUM_THRESHOLD, QSO_CHECKSUM_THRESHOLD, QSO_ELEMENT_THRESHOLD,
                          verbose);

        

        curq = 0;

        ///////////
        // Test Qmo
        ///////////

        if(batchsize)
            buffsize = batchsize * nmo2;

        // Note - nmo = nso in this case, so no need to resize mat
        //outbuf = unique_ptr<double[]>(new double[buffsize]);
        //mat = unique_ptr<double[]>(new double[naux*nmo2]);
        //dft.SetOutputBuffer(outbuf.get(), buffsize);
 

        while((n = dft.GetBatch_Qmo()))
        {
            std::copy(outbuf.get(), outbuf.get() + n*nmo2, mat.get() + curq*nmo2);
            curq += n;
        }

        ret += TestMatrix("QMO", qmo_filename,
                          mat.get(), matsize,
                          QMO_SUM_THRESHOLD, QMO_CHECKSUM_THRESHOLD, QMO_ELEMENT_THRESHOLD,
                          verbose);

        dft.ResetBatches();

        std::fill(mat.get(), mat.get()+matsize, 99);
        dft2.GetQBatch_Qmo(mat.get(), matsize, 0);
        ret += TestMatrix("QMO2", qmo_filename,
                          mat.get(), matsize,
                          QMO_SUM_THRESHOLD, QMO_CHECKSUM_THRESHOLD, QMO_ELEMENT_THRESHOLD,
                          verbose);

        curq = 0;


        ///////////
        // Test Qov
        ///////////
        int nvir = nmo - nocc;
        matsize = nocc*nvir*naux;

        if(batchsize)
            buffsize = batchsize * nocc*nvir;

        outbuf = unique_ptr<double[]>(new double[buffsize]);
        mat = unique_ptr<double[]>(new double[matsize]);
        dft.SetOutputBuffer(outbuf.get(), buffsize);


        while((n = dft.GetBatch_Qov()))
        {
            std::copy(outbuf.get(), outbuf.get() + n*nocc*nvir, mat.get() + curq*nocc*nvir);
            curq += n;
        }

        ret += TestMatrix("QOV", qov_filename,
                          mat.get(), matsize,
                          QMO_SUM_THRESHOLD, QMO_CHECKSUM_THRESHOLD, QMO_ELEMENT_THRESHOLD,
                          verbose);

        dft.ResetBatches();

        std::fill(mat.get(), mat.get()+matsize, 99);
        dft2.GetQBatch_Qov(mat.get(), matsize, 0);
        ret += TestMatrix("QOV2", qov_filename,
                          mat.get(), matsize,
                          QMO_SUM_THRESHOLD, QMO_CHECKSUM_THRESHOLD, QMO_ELEMENT_THRESHOLD,
                          verbose);

        curq = 0;




        ///////////
        // Test Qoo
        ///////////
        matsize = nocc*nocc*naux;

        if(batchsize)
            buffsize = batchsize * nocc*nocc;

        outbuf = unique_ptr<double[]>(new double[buffsize]);
        mat = unique_ptr<double[]>(new double[matsize]);
        dft.SetOutputBuffer(outbuf.get(), buffsize);

        while((n = dft.GetBatch_Qoo()))
        {
            std::copy(outbuf.get(), outbuf.get() + n*nocc*nocc, mat.get() + curq*nocc*nocc);
            curq += n;
        }

        ret += TestMatrix("QOO", qoo_filename,
                          mat.get(), matsize,
                          QMO_SUM_THRESHOLD, QMO_CHECKSUM_THRESHOLD, QMO_ELEMENT_THRESHOLD,
                          verbose);

        dft.ResetBatches();

        std::fill(mat.get(), mat.get()+matsize, 99);
        dft2.GetQBatch_Qoo(mat.get(), matsize, 0);
        ret += TestMatrix("QOO2", qoo_filename,
                          mat.get(), matsize,
                          QMO_SUM_THRESHOLD, QMO_CHECKSUM_THRESHOLD, QMO_ELEMENT_THRESHOLD,
                          verbose);

        curq = 0;



        ///////////
        // Test Qvv
        ///////////
        matsize = nvir*nvir*naux;

        if(batchsize)
            buffsize = batchsize * nvir*nvir;

        outbuf = unique_ptr<double[]>(new double[buffsize]);
        mat = unique_ptr<double[]>(new double[matsize]);
        dft.SetOutputBuffer(outbuf.get(), buffsize);

        while((n = dft.GetBatch_Qvv()))
        {
            std::copy(outbuf.get(), outbuf.get() + n*nvir*nvir, mat.get() + curq*nvir*nvir);
            curq += n;
        }

        ret += TestMatrix("QVV", qvv_filename,
                          mat.get(), matsize,
                          QMO_SUM_THRESHOLD, QMO_CHECKSUM_THRESHOLD, QMO_ELEMENT_THRESHOLD,
                          verbose);

        dft.ResetBatches();

        std::fill(mat.get(), mat.get()+matsize, 99);
        dft2.GetQBatch_Qvv(mat.get(), matsize, 0);
        ret += TestMatrix("QVV2", qvv_filename,
                          mat.get(), matsize,
                          QMO_SUM_THRESHOLD, QMO_CHECKSUM_THRESHOLD, QMO_ELEMENT_THRESHOLD,
                          verbose);

        curq = 0;







        cout << "\n\n";
        cout << "*************************************************\n"
             << "*************************************************\n"
             << "OVERALL RESULT: " << (ret ? "FAIL" : "PASS") << "\n";
        if(ret)
            cout << "    ( " << ret << " failures)\n";
        cout << "*************************************************\n"
             << "*************************************************\n";


        dft2.PrintTimings();


    }

    catch(const exception & ex)
    {
        cout << "\n*****************"
             << "\nCAUGHT EXCEPTION!"
             << "\n" << ex.what()
             << "\n*****************"
             << "\n\n";

        ret = -1;
    }

    return ret;
}





