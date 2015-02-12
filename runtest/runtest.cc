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
#include "panache/SimpleMatrix.h"
#include "panache/Output.h"
#include "panache/c_convert.h"
#include "panache/Flags.h"
#include "panache/Iterator.h"

#define CHOLESKY_DELTA 1e-3

#define QSO_ELEMENT_THRESHOLD 1e-11
#define QSO_SUM_THRESHOLD 1e-8
#define QSO_CHECKSUM_THRESHOLD 1.0

#define QMO_ELEMENT_THRESHOLD 1e-9
#define QMO_SUM_THRESHOLD 1e-6
#define QMO_CHECKSUM_THRESHOLD 2.0

using namespace panache;
using namespace std;

#include "panache/Parallel.h"

#ifdef PANACHE_MPI
#include <mpi.h>
#endif

ostream * out;

void PrintUsage(void)
{
    *out << "\n"
         << "Libpanache Testing Utility\n"
         << "\n"
         << "Usage: runtest [opt] <dir>\n"
         << "\n"
         << "Options:\n"
         << "-v           Verbose printing\n"
         << "-d           Write Q tensors to disk (rather than in core)\n"
         << "-k           Keep Q tensors on disk when done\n"
         << "-c           Use Cyclops Tensor Framework\n"
         << "-b           Number of batches to get at a time (default = all)\n"
         << "-t           Use transpose of C matrix\n"
         << "-g           Generate tests from basis/molecule info\n"
         << "-C           Disable cholesky runs\n"
         << "-S           Skip testing (useful for benchmarking)\n"
         << "-X           Skip getting batches + testing (useful for benchmarking)\n"
         << "-r           Read tensor from disk\n"
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
    *out << std::setw(35) << std::left << name;
    *out << std::setw(15) << std::left << ref;
    *out << std::setw(15) << std::left << run;
    *out << std::setw(15) << std::left << diff;
    *out << std::setw(15) << std::left << thresh;
    *out << std::setw(15) << std::left << result << "\n";
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
    *out << std::setw(35) << std::left << name;
    *out << std::setw(15) << std::left << ref;
    *out << std::setw(15) << std::left << run;
    *out << std::setw(15) << std::left << diff;
    *out << std::setw(15) << std::left << thresh;
    *out << std::setw(15) << std::left << result << "\n";
}


/*! \brief Prints a line of 104 dashes
 *
 * \param [in] out Stream to print to
 */
static void PrintSeparator(void)
{
    *out << std::string(104, '-') << std::endl;
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
    *out << "***********************************************************************\n";
    *out << "Matrix \"" << title << "\"" << "\n";
    *out << "***********************************************************************\n";

    PrintRow("Name", "Reference", "This run", "Diff", "Threshold", "Pass/Fail");
    PrintSeparator();


    nfailures += TestAndPrint("# of elements", matsize, tmi.nrow*tmi.ncol, 0, verbose);

    // DON'T TEST ELEMENTS IF DIMENSIONS ARE WRONG
    if(nfailures > 0)
    {
        *out << "... not continuing with matrix test since dimensions are wrong!\n\n\n\n\n\n";
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
    nfailures += TestAndPrint("checksum", checksum, tmi.checksum, checksum_threshold, verbose);

    for(int i = 0; i < 100; i++)
    {
        stringstream ss;
        ss << "Element " << tmi.elements[i].first;

        nfailures += TestAndPrint(ss.str(), mat[tmi.elements[i].first], tmi.elements[i].second,
                                  element_threshold, verbose);
    }

    if(!verbose && nfailures == 0)
        *out << "           ... no failed tests!\n";
    *out << "\n";
    *out << "***********************************************************************\n";
    *out << "Matrix \"" << title << "\" result: " << (nfailures ? "FAIL" : "PASS");

    if(nfailures)
        *out << " (" << nfailures << " failures)";

    *out << "\n";
    *out << "***********************************************************************\n";
    *out << "\n\n\n\n\n";
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




int RunTestMatrix(ThreeIndexTensor & dft, const string & title,
                  int batchsize, int tensorflag,
                  const string & reffile,
                  double sum_threshold, double checksum_threshold, double element_threshold,
                  bool skiptest, bool verbose)
{
    int ret = 0;
    int naux, ndim1, ndim2;
    dft.TensorDimensions(tensorflag, naux, ndim1, ndim2);

    // may be packed
    int ndim12 = dft.QBatchSize(tensorflag);

    int matsize = naux * ndim1 * ndim2; // always expanded
    int bufsize = ndim12 * naux;
    if(batchsize)
        bufsize = ndim12 * batchsize;
        
    unique_ptr<double[]> mat(new double[matsize]);
    unique_ptr<double[]> outbuf(new double[bufsize]);
    //std::fill(mat.get(), mat.get()+matsize, 0.0);
    //std::fill(outbuf.get(), outbuf.get()+bufsize, 0.0);

    // First, do by q
    ThreeIndexTensor::IteratedQTensorByQ iqtq = dft.IterateByQ(tensorflag, outbuf.get(), bufsize);
    while(iqtq)
    {
        int curq = iqtq.q();

        // tests are always done on unpacked matrices
        if(dft.IsPacked(tensorflag))
        {
            for(int i = 0; i < ndim1; i++)
            for(int j = 0; j <= i; j++)
                mat[curq*ndim1*ndim2+i*ndim2+j] 
              = mat[curq*ndim1*ndim2+j*ndim1+i] 
              = iqtq[i*(i+1)/2+j];
        }
        else
            std::copy(iqtq.Get(), iqtq.Get() + ndim12, mat.get() + curq*ndim12);

        ++iqtq;
    }

    string titleq(title);
    titleq.append(" (by Q)");
    if(skiptest)
        *out << "** Skipping testing " << titleq << "\n\n";
    else
        ret += TestMatrix(titleq, reffile,
                          mat.get(), matsize,
                          sum_threshold, checksum_threshold, element_threshold,
                          verbose);

    // Now do by ij
    bufsize = matsize;
    if(batchsize)
        bufsize = naux * batchsize;

    outbuf = unique_ptr<double[]>(new double[bufsize]);
    //std::fill(mat.get(), mat.get()+matsize, 0.0);
    //std::fill(outbuf.get(), outbuf.get()+bufsize, 0.0);


    // Note - The reference matrices are always stored "by q". So some
    // index math is appropriate
    ThreeIndexTensor::IteratedQTensorByIJ iqtij = dft.IterateByIJ(tensorflag, outbuf.get(), bufsize);
    while(iqtij)
    {
        int i = iqtij.i();
        int j = iqtij.j();

        // tests are always done on unpacked matrices
        if(dft.IsPacked(tensorflag))
        {
            // expand
            for(int q = 0; q < naux; q++)
                mat[q*ndim1*ndim2 + i * ndim2 + j]
              = mat[q*ndim1*ndim2 + j * ndim1 + i]
              = iqtij[q];
        }
        else
        {
            for(int q = 0; q < naux; q++)
                mat[q*ndim1*ndim2 + i * ndim2 + j] = iqtij[q];
        }

        ++iqtij;
    }

    string titleij(title);
    titleij.append(" (by IJ)");
    if(skiptest)
        *out << "** Skipping testing " << titleij << "\n\n";
    else
        ret += TestMatrix(titleij, reffile,
                          mat.get(), matsize,
                          sum_threshold, checksum_threshold, element_threshold,
                          verbose);

    return ret;

}


void GenTestMatrix(ThreeIndexTensor & dft, const string & title,
                  int tensorflag, int batchsize,
                  const string & reffile,
                  bool verbose)
{
    using namespace std;

    *out << "** Generating Matrix \"" << title << "\"" << "\n";

    // starts of the same as RunTestMatrix
    int naux, ndim1, ndim2;
    dft.TensorDimensions(tensorflag, naux, ndim1, ndim2);

    // may be packed
    int ndim12 = dft.QBatchSize(tensorflag);

    int bufsize = ndim12 * naux;
    if(batchsize)
        bufsize = ndim12 * batchsize;
        
    unique_ptr<double[]> outbuf(new double[bufsize]);
    unique_ptr<double[]> outbuf2(new double[ndim1*ndim2]); // expanded piece


    // only go by q

    double sum = 0;
    double checksum = 0;
    size_t index = 0;

    array<pair<size_t, double>, 50> largest, smallest;
    for(int i = 0; i < 50; i++)
    {
        largest[i].second = -numeric_limits<double>::max();
        smallest[i].second = numeric_limits<double>::max();
    }

    ThreeIndexTensor::IteratedQTensorByQ iqtq = dft.IterateByQ(tensorflag, outbuf.get(), bufsize);
    while(iqtq)
    {
        // tests are always done on unpacked matrices
        if(dft.IsPacked(tensorflag))
        {
            for(int i = 0; i < ndim1; i++)
            for(int j = 0; j <= i; j++)
                outbuf2[i*ndim1+j] 
              = outbuf2[j*ndim1+i] 
              = iqtq[i*(i+1)/2+j];
        }
        else
            copy(iqtq.Get(), iqtq.Get() + ndim12, outbuf2.get());

        // now we have the expanded matrices. Check for the largest & smallest elements
        // and update the checksums
        for(int64_t i = 0; i < ndim1*ndim2; i++)
        {
            double val = outbuf2[i];
            sum += val;
            checksum += val*static_cast<double>(index+1);

            if(val > largest[0].second)
            {
                largest[0].first = index;
                largest[0].second = val;
                sort(largest.begin(), largest.end(), [](pair<size_t, double> p1, pair<size_t, double> p2)
                {
                    return p1.second < p2.second;
                });
            }
            if(val < smallest[49].second)
            {
                smallest[49].first = index;
                smallest[49].second = val;
                sort(smallest.begin(), smallest.end(), [](pair<size_t, double> p1, pair<size_t, double> p2)
                {
                    return p1.second < p2.second;
                });
            }

            index++;
        }

        ++iqtq;
    }


    // write out the matrix
    ofstream outfile(reffile.c_str(), ofstream::trunc);
    outfile << naux << "\n" << ndim1*ndim2 << "\n";
    outfile << setprecision(20);
    outfile << sum << "\n" << checksum << "\n";
    for(auto & it : largest)
        outfile << it.first << " " << it.second << "\n";
    for(auto & it : smallest)
        outfile << it.first << " " << it.second << "\n";

    // file closed automatically
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
        bool byq = false;
        bool transpose = false;
        int batchsize = 0;
        bool cyclops = false;
        bool disk = false;
        bool docholesky = true;
        bool generate = false;
        bool skiptest = false;
        bool skipgetbatch = false;
        bool keepdisk = false;
        bool readdisk = false;

        int i = 1;
        while(i < argc)
        {
            string starg(GetNextArg(i, argc, argv));
            if(starg == "-v")
                verbose = true;
            else if(starg == "-t")
                transpose = true;
            else if(starg == "-q")
                byq = true;
            else if(starg == "-b")
                batchsize = GetIArg(i, argc, argv);
            else if(starg == "-C")
                docholesky = false;
            else if(starg == "-S")
                skiptest = true;
            else if(starg == "-d")
                disk = true;
            else if(starg == "-k")
                keepdisk = true;
            else if(starg == "-r")
                readdisk = true;
            else if(starg == "-c")
                cyclops = true;
            else if(starg == "-g")
                generate = true;
            else if(starg == "-X")
            {
                skipgetbatch = true;
                skiptest = true;
            }
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

        if(generate && !docholesky)
            throw std::runtime_error("Generate must include cholesky!");

        if(generate && skiptest)
            throw std::runtime_error("Generate and skiptest doesn't make any sense!");

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
        if(generate)
            *out << "Generating tests for : " << desc_str << "\n";
        else
            *out << "Results for test: " << desc_str << "\n";
        *out << "--------------------------------------\n";

        auto mol = ReadMoleculeFile(dir + "geometry");
        auto primary = ReadBasisFile(mol, dir + "basis.primary");
        auto cmat = ReadCMatrixFile(dir + "cmat");

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
        int nocc = ReadNocc(dir + "nocc");
        int nmo = nso;

        DFTensor dft(primary, dir + "basis.aux.gbs", "/tmp/df",
                     DFOPT_COULOMB | DFOPT_EIGINV,
                     BSORDER_PSI4, 0);

        // *** We are only testing Qso from CHTensor               *** //
        // *** But generating them all (to test for memory issues) *** //
        CHTensor cht(primary, CHOLESKY_DELTA, "/tmp/ch", BSORDER_PSI4, 0);

        dft.SetCMatrix(cmat->pointer(), nmo, transpose);
        cht.SetCMatrix(cmat->pointer(), nmo, transpose);
        dft.SetNOcc(nocc);
        cht.SetNOcc(nocc);

        int dfqflags = (QGEN_QSO | QGEN_QMO | QGEN_QOO | QGEN_QOV | QGEN_QVV);
        int chqflags = (QGEN_QSO | QGEN_QMO | QGEN_QOO | QGEN_QOV | QGEN_QVV);

        int qstore = 0;
        if(byq)
            qstore |= QSTORAGE_BYQ;

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

        if(docholesky)
            cht.GenQTensors(chqflags, qstore);


        if(generate)
        {
            ///////////
            // Gen Qso
            ///////////
            GenTestMatrix(dft, "QSO", QGEN_QSO, batchsize,
                          dir + "qso", verbose);
    
            ///////////
            // Gen Qmo
            ///////////
            GenTestMatrix(dft, "QMO", QGEN_QMO, batchsize,
                          dir + "qmo", verbose);
    
            ///////////
            // Gen Qoo
            ///////////
            GenTestMatrix(dft, "QOO", QGEN_QOO, batchsize,
                          dir + "qoo", verbose);
    
            ///////////
            // Gen Qov
            ///////////
            GenTestMatrix(dft, "QOV", QGEN_QOV, batchsize,
                          dir + "qov", verbose);
    
            ///////////
            // Gen Qvv
            ///////////
            GenTestMatrix(dft, "QVV", QGEN_QVV, batchsize,
                          dir + "qvv", verbose);

            ///////////////////////
            // Test Cholesky QSO
            ///////////////////////
            GenTestMatrix(cht, "CHQSO", QGEN_QSO, batchsize,
                           dir + "chqso", verbose);
        }
        else if(!skipgetbatch)
        {
            ///////////
            // Test Qso
            ///////////
            ret += RunTestMatrix(dft, "QSO",
                                 batchsize, QGEN_QSO,
                                 dir + "qso", 
                                 QSO_SUM_THRESHOLD, QSO_CHECKSUM_THRESHOLD, QSO_ELEMENT_THRESHOLD,
                                 skiptest, verbose);
    
            ///////////
            // Test Qmo
            ///////////
            ret += RunTestMatrix(dft, "QMO",
                                 batchsize, QGEN_QMO,
                                 dir + "qmo", 
                                 QMO_SUM_THRESHOLD, QMO_CHECKSUM_THRESHOLD, QMO_ELEMENT_THRESHOLD,
                                 skiptest, verbose);
    
            ///////////
            // Test Qoo
            ///////////
            ret += RunTestMatrix(dft, "QOO",
                                 batchsize, QGEN_QOO,
                                 dir + "qoo", 
                                 QMO_SUM_THRESHOLD, QMO_CHECKSUM_THRESHOLD, QMO_ELEMENT_THRESHOLD,
                                 skiptest, verbose);
    
            ///////////
            // Test Qov
            ///////////
            ret += RunTestMatrix(dft, "QOV",
                                 batchsize, QGEN_QOV,
                                 dir + "qov",
                                 QMO_SUM_THRESHOLD, QMO_CHECKSUM_THRESHOLD, QMO_ELEMENT_THRESHOLD,
                                 skiptest, verbose);
    
            ///////////
            // Test Qvv
            ///////////
            ret += RunTestMatrix(dft, "QVV",
                                 batchsize, QGEN_QVV,
                                 dir + "qvv",
                                 QMO_SUM_THRESHOLD, QMO_CHECKSUM_THRESHOLD, QMO_ELEMENT_THRESHOLD,
                                 skiptest, verbose);
    
            ///////////////////////
            // Test Cholesky QSO
            ///////////////////////
            if(docholesky)
            {
                ret += RunTestMatrix(cht, "CHQSO",
                                     batchsize, QGEN_QSO,
                                     dir + "chqso",
                                     QSO_SUM_THRESHOLD, QSO_CHECKSUM_THRESHOLD, QSO_ELEMENT_THRESHOLD,
                                     skiptest, verbose);
            }
        }


        if(!skiptest)
        {
            *out << "\n\n"
                 << "*************************************************\n"
                 << "*************************************************\n";
    
            if(generate)
            {
                *out << " REFERENCE FILES GENERATED\n";
            }
            else
            {
                *out << "OVERALL RESULT: " << (ret ? "FAIL" : "PASS") << "\n";
                if(ret)
                    *out << "    ( " << ret << " failures)\n";
            }
    
            *out << "*************************************************\n"
                 << "*************************************************\n";
        }

        // even if not verbose, print the timints
        if(!verbose)
            panache::output::SetOutput(&*out);

        dft.PrintTimings();

        if(docholesky)
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





