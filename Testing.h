#ifndef PANACHE_TESTING_H
#define PANACHE_TESTING_H

#define TEST_EXPONENT_THRESHOLD 1e-11
#define TEST_COEFFICIENT_THRESHOLD 1e-11
#define TEST_MOLECULE_XYZ_THRESHOLD 1e-11
#define TEST_MOLECULE_Z_THRESHOLD 1e-11
#define TEST_MOLECULE_MASS_THRESHOLD 1e-11


#include <string>
#include <vector>
#include <cmath>
#include <array>
#include <iostream>
#include <iomanip>
#include "c_interface.h"

using std::vector;
using std::string;
using std::array;

namespace panache
{
namespace testing
{

class TestInfo
{
private:

    template<typename data_type>
    struct TestResult
    {
        data_type reference;
        data_type result;
        data_type diff;
        data_type threshold;


        TestResult(data_type ref, data_type thr)
        {
            diff = result = static_cast<data_type>(0);
            set(ref, thr);
        }


        TestResult(data_type ref)
            : TestResult(ref, static_cast<data_type>(0))
        {
        }

        TestResult()
        {
            reference = result = diff = threshold = static_cast<data_type>(0);
        }

        void set(data_type ref,
                 data_type thr = static_cast<data_type>(0))
        {
            reference = ref;
            threshold = thr;
            diff = result = static_cast<data_type>(0);
        }


        bool test(void) const
        {
            if(diff > threshold)
                return false;
            else
                return true;
        }

        void set_result(data_type r)
        {
            result = r;
            diff = abs(result - reference);
        }
    };


    struct BasisTest
    {
        TestResult<int> total_nshell;    // total # of shells
        TestResult<int> total_nprim;     // total # of primitives
        TestResult<int> total_ncenters;      // # of centers

        vector<TestResult<int>> center_nshell; // # of shells per center
        vector<TestResult<int>> shell_nprim;  // # of primitives per shell
        vector<TestResult<double>> exp; // exponents of the primitives
        vector<TestResult<double>> coef; // coefficients of the primitives
    };


    struct MoleculeTest
    {
        TestResult<int> total_ncenters;      // # of centers
        array<vector<TestResult<double>>, 3> xyz;
        vector<TestResult<double>> Z;
        vector<TestResult<double>> mass;
    };


    // Name of the test
    string testname_;

    BasisTest primary_test_, aux_test_;
    MoleculeTest molecule_test_;

    int ncenters_, primary_nshells_, aux_nshells_;

    // basis set information
    int * primary_nshellspercenter_; // of length ncenter
    int * aux_nshellspercenter_; // of length ncenter
    C_ShellInfo * primary_shells_;
    C_ShellInfo * aux_shells_;

    // molecule
    C_AtomCenter * atoms_;


    // Parsers
    static int ReadBasisFile(const string & filename,
                             int * &nshellspercenter,
                             C_ShellInfo * &shells,
                             BasisTest & test);

    static int ReadMolecule(const string & filename,
                            C_AtomCenter * &atoms,
                            MoleculeTest & test);


    // Printing of results
    template<typename T>
    static void PrintColumn(std::ostream & out, const T & val, int width = 15)
    {
        out << std::setw(width) << std::left << val;
    }

    template<typename T>
    static void PrintRow(std::ostream & out,
                                const std::string name,
                                const T & ref,
                                const T & run,
                                const T & diff,
                                const T & thresh,
                                const std::string & result)
    {
        PrintColumn(out, name, 25);
        PrintColumn(out, ref);
        PrintColumn(out, run);
        PrintColumn(out, diff);
        PrintColumn(out, thresh);
        PrintColumn(out, result);
        out << std::endl;
    }

    // for printing headers
    static void PrintRow(std::ostream & out,
                                const char * name,
                                const char * ref,
                                const char * run,
                                const char * diff,
                                const char * thresh,
                                const char * result)
    {
        PrintColumn(out, name, 25);
        PrintColumn(out, ref);
        PrintColumn(out, run);
        PrintColumn(out, diff);
        PrintColumn(out, thresh);
        PrintColumn(out, result);
        out << std::endl;
    }

    static void PrintSeparator(std::ostream & out)
    {
        out << std::string(100, '-') << std::endl;
    }

    template<typename T>
    static bool Test(std::ostream & out,
                            const std::string & name,
                            TestResult<T> & test)
    {
        const char * pass = test.test() ? "pass" : "FAIL";
        PrintRow(out, name,
                 test.reference,
                 test.result,
                 test.diff,
                 test.threshold,
                 pass);

        return test.test();

    }


    // disable copying and assignment
    TestInfo & operator=(const TestInfo & rhs);
    TestInfo(const TestInfo & t);

public:
    TestInfo(const std::string & testname, const std::string & dir);

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

    void TestBasisConversion(void);
    void TestMoleculeConversion(void);
    bool PrintResults(std::ostream & out, bool verbose = false); 

    ~TestInfo();
};




}
} //close namespace panache::testing

#endif

