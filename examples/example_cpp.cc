#include <iostream>
#include <fstream>
#include <sstream>

#include "panache/DFTensor.h"
#include "panache/CHTensor.h"
#include "panache/Output.h"
#include "panache/c_convert.h"
#include "panache/SimpleMatrix.h"
#include "examples/common_cpp.h"

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
         << "Usage: example_cpp [opt] <dir>\n"
         << "\n"
         << "Options:\n"
         << "-v           Verbose printing\n"
         << "-d           Write Q tensors to disk (rather than in core)\n"
         << "-c           Cholesky (default is density fitting)\n" 
         << "-r           Read tensor from disk\n"
         << "-h           Print help (you're looking at it\n"
         << "<dir>        Directory holding the test information\n"
         << "\n\n";
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
        bool disk = false;
        bool keepdisk = false;
        bool readdisk = false;
        bool cholesky = false;

        int i = 1;
        while(i < argc)
        {
            string starg(GetNextArg(i, argc, argv));
            if(starg == "-v")
                verbose = true;
            else if(starg == "-c")
                cholesky = true;
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

        std::unique_ptr<ThreeIndexTensor> ten;
        if(cholesky)
        {
            ten = std::unique_ptr<ThreeIndexTensor>(
                         new CHTensor(primary, CHOLESKY_DELTA, "/tmp/ch",
                         BSORDER_PSI4, 0)
                         );
        }
        else
        {
            ten = std::unique_ptr<ThreeIndexTensor>(
                         new DFTensor(primary, dir + "basis.aux.gbs", "/tmp/df",
                         DFOPT_COULOMB | DFOPT_EIGINV,
                         BSORDER_PSI4, 0)
                         );
        }

        ten->SetCMatrix(cmat->pointer(), nmo, false);
        ten->SetNOcc(nocc);

        int qflags = QGEN_QOV;

        int qstore = 0;

        if(disk)
            qstore |= QSTORAGE_ONDISK;
        if(keepdisk)
            qstore |= QSTORAGE_KEEPDISK;
        if(readdisk)
            qstore |= QSTORAGE_READDISK;
        else
            qstore |= QSTORAGE_INMEM;

        ten->GenQTensors(qflags, qstore);


        int naux, ndim1, ndim2;

        // ignoring return value
        ten->TensorDimensions(QGEN_QOV, naux, ndim1, ndim2);


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
            int ia = ten->CalcIndex(QGEN_QOV, i, a);
            ten->GetBatch(QGEN_QOV, buf_iap, naux, ia);

            for(int j = 0; j < nocc; j++)
            for(int b = 0; b < nvir; b++)
            {
                int jb = ten->CalcIndex(QGEN_QOV, j, b);
                ten->GetBatch(QGEN_QOV, buf_jbp, naux, jb);
    
                int ib = ten->CalcIndex(QGEN_QOV, i, b);
                int ja = ten->CalcIndex(QGEN_QOV, j, a);
    
                ten->GetBatch(QGEN_QOV, buf_ibp, naux, ib);
                ten->GetBatch(QGEN_QOV, buf_jap, naux, ja);
    
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

        ten->PrintTimings();

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





