#ifndef PANACHE_TESTING_H
#define PANACHE_TESTING_H

#define TEST_EXPONENT_THRESHOLD 1e-11
#define TEST_COEFFICIENT_THRESHOLD 1e-11
#define TEST_MOLECULE_XYZ_THRESHOLD 1e-11
#define TEST_MOLECULE_Z_THRESHOLD 1e-11
#define TEST_MOLECULE_MASS_THRESHOLD 1e-11
#define TEST_QSO_ELEMENT_THRESHOLD 1e-8

#include <stdexcept>
#include <string>
#include <vector>
#include <cmath>
#include <array>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <memory>
#include "c_interface.h"

using std::vector;
using std::string;
using std::array;

namespace panache
{

class Matrix;
typedef std::shared_ptr<Matrix> SharedMatrix;

namespace testing
{


/*! \brief An exception class thrown if there
 *         is an error parsing a test file
 */
class TestingParserException : public std::exception
{
private:
    string str_;        //!< Some descriptive string, also containing the filename


public:
    /*! \brief Constructor
     *
     * \param [in] desc Some descriptive string
     * \param [in] filename The file in which the error occurred
     */
    TestingParserException(const std::string & desc, const std::string & filename)
    {
        std::stringstream ss;
        ss << desc << "\nFilename: \'" << filename << "\'";
        str_ = ss.str();
    }

    /*! \brief Prints the exception information
     */
    virtual const char * what() const noexcept
    {
        return str_.c_str();
    }
};





typedef std::runtime_error TestingSanityException; //!< Used for internal checks of the TestingInfo class





/*! \brief Holds the info for a single test
 *
 *  Each object will hold the information for testing a single number, the
 *  type of which depends on data_type
 *
 *  \tparam data_type The type of data point. Usually int or double.
 */
template<typename data_type>
struct TestResult
{
    data_type reference;     //!< The reference value to be compared against
    data_type thisrun;       //!< The data actually obtained
    data_type diff;          //!< The difference
    data_type threshold;     //!< The threshold above which the test fails


    /*! \brief Constructor
     *
     * Takes an optional threshold, which by default is zero.
     *
     * \param [in] ref The reference (expected) value
     * \param [in] thr The threshold above which the test would fail.
     */
    TestResult(data_type ref, data_type thr = static_cast<data_type>(0))
    {
        diff = thisrun = static_cast<data_type>(0);
        set(ref, thr);
    }


    /*! \brief Default constructor
     */
    TestResult()
    {
        reference = thisrun = diff = threshold = static_cast<data_type>(0);
    }


    /*! \brief Sets the reference and threshold information
     *
     *  Unless specified, the threshold is zero.
     */
    void set(data_type ref,
             data_type thr = static_cast<data_type>(0))
    {
        reference = ref;
        threshold = thr;
        diff = thisrun = static_cast<data_type>(0);
    }


    /*! \brief Check the test to see if it passes
     *
     * \return True if the test passes (difference is less
     *         than the threshold), false otherwise
     */
    bool check(void) const
    {
        if(diff > threshold)
            return false;
        else
            return true;
    }


    /*! \brief Set the results obtained for a run
     *
     *  Will also calculate the (absolute) difference.
     */
    void set_thisrun(data_type r)
    {
        thisrun = r;
        diff = (thisrun - reference);
        if(diff < 0)
            diff *= static_cast<data_type>(-1);
    }
};


/*! \brief Specialization of TestResult for string comparison
 */
template<>
struct TestResult<std::string>
{
    std::string reference;     //!< The reference string to be compared against
    std::string thisrun;       //!< The string actually obtained
    bool diff;                 //!< True if different, false if same


    /*! \brief Constructor
     *
     * \param [in] ref The reference (expected) string
     */
    TestResult(std::string ref)
    {
        diff = false;
        set(ref);
    }


    /*! \brief Default constructor
     */
    TestResult()
    {
        reference = thisrun = "";
        diff = false;
    }


    /*! \brief Sets the reference string
     */
    void set(std::string ref)
    {
        reference = ref;
        diff = false;
    }


    /*! \brief Check the test to see if it passes
     *
     * \return True if the test passes (strings are the same), false otherwise
     */
    bool check(void) const
    {
        return (!diff);
    }


    /*! \brief Set the results obtained for a run
     *
     *  Will also check for string differences
     */
    void set_thisrun(std::string r)
    {
        diff = (thisrun != reference);
    }
};


typedef TestResult<double> DTestResult;  //!< Test result for a double type
typedef TestResult<int> ITestResult;     //!< Test result for a (signed) integer type
typedef TestResult<std::string> STestResult;     //!< Test result for a string





/*! \brief A object containing a test for a single molecule/basis
 *         set.
 *
 *   This class tests:
 *         - Generation of a Molecule object
 *         - Generation of a BasisSet object
 *             - Primary basis set
 *             - Auxiliary basis set
 *         - Generation of 3-index integrals
 *         - Generation of the Q matrix
 */
class TestInfo
{
private:



    /*! \brief Holds information for a test of a single basis set conversion
     *
     * The purpose is to test the code that generates a BasisSet object
     * from C arrays or from other information.
     */
    struct BasisTest
    {
        ITestResult nshell;    //!<  The total number of shells
        ITestResult nbf;       //!<  The total number of basis functions (possibly spherical)
        ITestResult nao;       //!<  Number of cartesian AOs
        ITestResult nprim;     //!<  Number of primitives

        vector<ITestResult> center_nshell;  //!< The number of shells on each center
        vector<ITestResult> shell_nprim;    //!< The number of primitives in each shell
        vector<ITestResult> shell_am;       //!< Angular momentum of each shell
        vector<ITestResult> shell_ispure;   //!< Is pure or not
        vector<DTestResult> exp;         //!< The exponents of the primitives
        vector<DTestResult> coef;        //!< Contraction coefficients of the primitives
    };


    /*! \brief Holds information for a test of a molecule conversion
     *
     * The purpose is to test the code that genrates a Molecule object
     * from C arrays or from other information.
     */
    struct MoleculeTest
    {
        ITestResult natom;        //!< Number of centers in the molecule
        ITestResult nallatom;    //!< Number of centers in the molecule (including dummies)
        STestResult schoen;      //!< Schoenflies symbol
        STestResult fullpg;      //!< Full point group symbol
        vector<DTestResult> Z;   //!< Z-number of each center
        vector<DTestResult> mass; //!< Mass of each center

        /*! \brief XYZ coordinates
         *
         * Elements of the array class correspond to x, y, or z, with the
         * vector being (hopefully) being the number of centers. ie
         *
         * \code{.cpp}
         * xyz[0][0] // x coordinate of center 0
         * xyz[1][0] // y coordinate of center 0
         * xyz[2][4] // z coordinate of center 4
         * \endcode
         */
        array<vector<DTestResult>, 3> xyz;
    };



    /*! \brief Holds information for testing a single matrix element
     */
    struct MatrixElementTest
    {
        size_t index;   //!< Index of the element
        DTestResult test; //!< Value at that index
    };


    /*! \brief Holds information for testing a matrix
     */
    struct MatrixTest
    {
        ITestResult nrow; //!< Number of rows
        ITestResult ncol; //!< Number of columns
        DTestResult sum; //!< Sum of the Qso matrix
        DTestResult checksum; //!< Checksum of the Qso matrix
        array<MatrixElementTest, 100> elements; //!< 100 elements of the Qso matrix
    };


    string testname_;  //!< Name of the test (usually "molecule basis")

    BasisTest primary_test_;    //!< Tests the primary basis set
    BasisTest aux_test_;        //!< Tests the auxiliary basis set
    MoleculeTest molecule_test_; //!< Tests the molecule

    MatrixTest qso_test_; //!< Tests the Qso matrix


    // These are used internally to know what to test, etc
    // Usually read in from the test files or otherwise deduced
    int ncenters_;   //!< Number of centers in the moleulce

    int primary_nshells_;  //!< Number of shells in the primary basis
    int aux_nshells_; //!< Number of shells in the auxiliary basis set

    int primary_nprim_; //!< Number of primitives in the primary basis set
    int aux_nprim_;     //!< Number of primitives in the auxiliary basis set

    int * primary_nshellspercenter_;  //!< Number of shells on each center for the primary basis set
    int * aux_nshellspercenter_;      //!< Number of shells on each center for the auxiliary basis set

    C_ShellInfo * primary_shells_;    //!< Shell information for the primary basis set
    C_ShellInfo * aux_shells_;        //!< Shell information for the auxiliary basis set

    C_AtomCenter * atoms_; //!< Specification of the molecule


    // Parsers
    /*! \brief Reads a basis set in from a file
     *
     *  The return value is the number of centers and the length of the nshellspercenter
     *  array. The length of the shells array is the sum of the elements of the nshellspercenter
     *  array.
     *
     *  \param [in] filename The full path to the file
     *  \param [out] nshellspercenter Pointer that will be set to newly-allocated memory
     *                                with each element corresponding to the number of shells
     *                                on that center. Must be deleted afterwards.
     *
     *  \param [out] shells Pointer that will be set to newly-allocated memory
     *                      with each element containing the information for a shell.
     *                      Must be deleted afterwards.
     *
     *  \param [in,out] test Used to store the information that was read in from the
     *                       file for later comparison with a BasisSet object.
     *
     *  \return The number of centers.
     */
    static int ReadBasisFile(const string & filename,
                             int * &nshellspercenter,
                             C_ShellInfo * &shells,
                             BasisTest & test);



    /*! \brief Reads molecule information from a file
     *
     *  The return value is the number of centers and the length of the atoms array.
     *
     *  \param [in] filename The full path to the file
     *  \param [out] atoms Pointer that will be set to newly-allocated memory
     *                     with each element corresponding to an atom or center.
     *
     *  \param [out] shells Pointer that will be set to newly-allocated memory
     *                      with each element containing the information for a shell.
     *                      Must be deleted afterwards.
     *
     *  \param [in,out] test Used to store the information that was read in from the
     *                       file for later comparison with a Molecule object.
     *
     *  \return The number of centers.
     */
    static int ReadMolecule(const string & filename,
                            C_AtomCenter * &atoms,
                            MoleculeTest & test);



    /*! \brief Reads matrix information a file
     *
     *  Reads the sum, checksum, and list of elements to test (with corresponding values).
     *
     *  \param [in] filename The full path to the file
     *  \param [out] test A MatrixTest object to put the test information in
     *  \prarm [in] threshold The desired threshold for checking matrix elements
     */
    static void ReadMatrixInfo(const string & filename,
                               MatrixTest & test, double threshold);


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
    static void PrintRow(std::ostream & out,
                         const std::string name,
                         const T & ref,
                         const T & run,
                         const T & diff,
                         const T & thresh,
                         const std::string & result)
    {
        out << std::setw(35) << std::left << name;
        out << std::setw(15) << std::left << ref;
        out << std::setw(15) << std::left << run;
        out << std::setw(15) << std::left << diff;
        out << std::setw(15) << std::left << thresh;
        out << std::setw(15) << std::left << result << "\n";
    }


    /*! \brief Overload for printing strings (usually header rows)
     *
     * Used to prevent ambiguities when compiling
     *
     * \param [in] out Stream to print to
     * \param [in] name Descriptive name of the test
     * \param [in] ref  The reference (expected) value
     * \param [in] run  The value for this run
     * \param [in] diff The absolute difference
     * \param [in] thresh Threshold for failing the test
     * \param [in] result String of the result (usually "pass" or "FAIL")
     */
    static void PrintRow(std::ostream & out,
                         const char * name,
                         const char * ref,
                         const char * run,
                         const char * diff,
                         const char * thresh,
                         const char * result)
    {
        out << std::setw(35) << std::left << name;
        out << std::setw(15) << std::left << ref;
        out << std::setw(15) << std::left << run;
        out << std::setw(15) << std::left << diff;
        out << std::setw(15) << std::left << thresh;
        out << std::setw(15) << std::left << result << "\n";
    }


    /*! \brief Prints a line of 100 dashes
     *
     * \param [in] out Stream to print to
     */
    static void PrintSeparator(std::ostream & out)
    {
        out << std::string(100, '-') << std::endl;
    }


    /*! \brief Checks the results of a test and prints a row with the results
     *
     * \param [in] out Stream to print to
     * \param [in] name A descriptive name of the test
     * \param [in] test The test to check and print
     * \return 0 if the test passes, otherwise 1
     */
    template<typename T>
    static int Test(std::ostream & out,
                    const std::string & name,
                    const TestResult<T> & test)
    {
        const char * pass = test.check() ? "pass" : "FAIL";
        PrintRow(out, name, test.reference, test.thisrun, test.diff, test.threshold, pass);

        if(test.check())
            return 0;
        else
            return 1;
    }


    /*! \brief Override for Test for testing strings
     *
     * \param [in] out Stream to print to
     * \param [in] name A descriptive name of the test
     * \param [in] test The test to check and print
     * \return 0 if the test passes, otherwise 1
     */
    static int Test(std::ostream & out,
                    const std::string & name,
                    const STestResult & test)
    {
        const char * pass = test.check() ? "pass" : "FAIL";
        const std::string diff = test.diff ? "T" : "F";
        PrintRow(out, name, test.reference, test.thisrun, diff, std::string("-"), pass);

        if(test.check())
            return 0;
        else
            return 1;
    }


    /*! \brief Tests the conversion of a single basis set
     *
     * Split out into its own function since we test two basis sets. Information
     * such as the number of centers is read from the class data members.
     *
     * \param [in] nshells The number of shells in the basis set
     * \param [in] nshellspercenter The number of shells on each center
     * \param [in] shells The information for each shell
     * \param [in] test A BasisTest object to place the results of this run
     */
    void TestBasisConversion(int nshells, int * nshellspercenter, C_ShellInfo * shells, BasisTest & test);



    /*! \brief Prints the results of a conversion of a single basis set
     *
     * Split out into its own function since we test two basis sets. Information
     * such as the number of centers is read from the class data members.
     *
     * \param [in] out The stream to print to
     * \param [in] type The type of basis ("primary" or "auxiliary")
     * \param [in] nshells The number of shells in the basis set
     * \param [in] nprim The number of primitives in the basis set
     * \param [in] test A BasisTest object to test
     * \param [in] verbose Verbose results
     */
    int PrintBasisResults(std::ostream & out, const string & type,
                          int nshells, int nprim,
                          const BasisTest & test,
                          bool verbose = false);



    /*! \brief Prints the results of a matrix test
     *
     * Split out into its own function since we may test multiple matrices.
     *
     * \param [in] out The stream to print to
     * \param [in] type The type of basis ("qso", etc)
     * \param [in] test A MatrixTest object to test
     * \param [in] verbose Verbose results
     */
    int PrintMatrixResults(std::ostream & out, const string & type,
                           const MatrixTest & test, bool verbose = false);


    /*! \brief Tests a matrix
     *
     * The test is carried out on the matrix as if it were an array, since the matrix is stored
     * in contiguous memory.
     *
     * \param [in] mat The matrix to test
     * \param [in] test The MatrixTest object to put the results in
     */
    static void TestMatrix(SharedMatrix mat,
                           MatrixTest & test);



    // disable copying and assignment
    TestInfo & operator=(const TestInfo & rhs) = delete;
    TestInfo(const TestInfo & t) = delete;


public:
    /*! \brief Constructor
     *
     *  \param [in] testname A descriptive name for the test (usually "molecule basis")
     *  \param [in] dir A directory containing the test files
     */
    TestInfo(const std::string & testname, const std::string & dir);

    /*! \brief Test conversion of C arrays into BasisSet objects
     *
     *  This tests both the primary and auxiliary basis sets.
     */
    void TestBasisConversion(void);


    /*! \brief Test conversion of C arrays into Molecule objects
     */
    void TestMoleculeConversion(void);


    /*! \brief Tests the Qso Matrix
     */
    void TestQsoMatrix(void);


    /*! \brief Print the results of all the tests.
     *
     *  Be sure to run them first!
     *
     *  \param [in] out Stream to print the results to
     *  \param [in] verbose Be verbose about everything
     *  \return The number of failed tests
     */
    int PrintResults(std::ostream & out, bool verbose = false);



#ifdef PANACHE_DEVELOPER_GENERATE
    /*! \brief Generate some test files
     *
     * Generates files related to matrix elements and ERI values
     */
    void Generate(void);
#endif



    ~TestInfo();






    ///////////////////////////////////////////////////////////////////
    // used only by me to bootstrap testing using information from psi
    ///////////////////////////////////////////////////////////////////
    C_ShellInfo * primary_shells(void)
    {
        return primary_shells_;
    }
    C_ShellInfo * aux_shells(void)
    {
        return aux_shells_;
    }

    int * primary_nshellspercenter(void)
    {
        return primary_nshellspercenter_;
    }
    int * aux_nshellspercenter(void)
    {
        return aux_nshellspercenter_;
    }

    C_AtomCenter * atoms(void)
    {
        return atoms_;
    }
    int ncenters(void)
    {
        return ncenters_;
    }
    ///////////////////////////////////////////////////////////////////

};

}
} //close namespace panache::testing

#endif


