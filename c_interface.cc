#include "c_interface.h"
#include "BasisSet.h"
#include "Molecule.h"
#include "DFTensor.h"

extern "C" {

    double * C_QAO(int ncenters, int nao, int nmo,
                 C_AtomCenter * atoms,
                 int * primary_nshellspercenter, struct C_ShellInfo * primary_shells,
                 int * aux_nshellspercenter, struct C_ShellInfo * aux_shells,
                 int * nrow_out, int * ncol_out)
    {
        // Molecule
        std::shared_ptr<panache::Molecule> molecule(new panache::Molecule());

        for(auto i = 0; i < ncenters; i++)
        {
            molecule->add_atom(
                atoms[i].Z,
                atoms[i].center[0],
                atoms[i].center[1],
                atoms[i].center[2],
                atoms[i].symbol,
                atoms[i].mass,
                0);
        }

        molecule->reset_point_group("c1");
        molecule->set_full_point_group();

        //std::cout << "New Molecule: Symmetry " << molecule->schoenflies_symbol() << "(full " << molecule->full_point_group(    ) << ")\n";
        //for(int i = 0; i < molecule->nallatom(); i++)
        //{
        //    panache::Vector3 c = molecule->xyz(i);
        //    std::cout << "Center " << i << " : " << c[0] << " , " << c[1] << " , " << c[2] << "\n";
        //}



        // Construct the basis set info
        std::vector<std::vector<panache::ShellInfo>> pshellmap, ashellmap;

        int primary_counter = 0;
        int aux_counter = 0;

        for(int i = 0; i < ncenters; ++i)
        {
            std::vector<panache::ShellInfo> ashellvec,pshellvec;
            panache::Vector3 cen;

            for(int j = 0; j < 3; j++)
                cen[j] = atoms[i].center[j];

            for(int j = 0; j < primary_nshellspercenter[i]; j++)
            {
                std::vector<double> c,e;
                for(int k = 0; k < primary_shells[primary_counter].nprim; k++)
                {
                    c.push_back(primary_shells[primary_counter].coef[k]);
                    e.push_back(primary_shells[primary_counter].exp[k]);
                }
                pshellvec.push_back(panache::ShellInfo(primary_shells[primary_counter].am, c, e,
                                                       primary_shells[primary_counter].ispure ? panache::ShellInfo::GaussianType::Pure : panache::ShellInfo::GaussianType::Cartesian,
                                                       i, cen, 0, panache::ShellInfo::PrimitiveType::Unnormalized));
                primary_counter++;
            }
            pshellmap.push_back(pshellvec);


            for(int j = 0; j < aux_nshellspercenter[i]; j++)
            {
                std::vector<double> c,e;
                for(int k = 0; k < aux_shells[aux_counter].nprim; k++)
                {
                    c.push_back(aux_shells[aux_counter].coef[k]);
                    e.push_back(aux_shells[aux_counter].exp[k]);
                }
                ashellvec.push_back(panache::ShellInfo(aux_shells[aux_counter].am, c, e,
                                                       aux_shells[aux_counter].ispure ? panache::ShellInfo::GaussianType::Pure : panache::ShellInfo::GaussianType::Cartesian,
                                                       i, cen, 0, panache::ShellInfo::PrimitiveType::Unnormalized));
                aux_counter++;
            }
            ashellmap.push_back(ashellvec);
        }

        std::shared_ptr<panache::BasisSet> primaryBasis(new panache::BasisSet(molecule, pshellmap));
        std::shared_ptr<panache::BasisSet> auxiliaryBasis(new panache::BasisSet(molecule, ashellmap));

        panache::DFTensor dft(primaryBasis, auxiliaryBasis, nao, nmo);
        auto mat = dft.Qso();

        // convert the resulting matrix
        *nrow_out = mat->rowspi()[0];
        *ncol_out = mat->colspi()[0]; 
        //return mat->give_up();

        double * ret = (double *)malloc(mat->rowspi()[0]*mat->colspi()[0] * sizeof(double));
        memcpy(ret, &(mat->pointer(0)[0][0]), mat->rowspi()[0]*mat->colspi()[0]*sizeof(double));
        return ret;

    }


    void free_matrix(double * mat)
    {
        free(mat);
    }

}

