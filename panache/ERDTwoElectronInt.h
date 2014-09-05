/*! \file
 * \brief Base class for libERD-based ERI (header)
 * \author Benjamin Pritchard (ben@bennyp.org)
 */

#ifndef PANACHE_ERDTWOELECTRONINT_H
#define PANACHE_ERDTWOELECTRONINT_H

#include "panache/TwoBodyAOInt.h"

typedef int F_INT;
typedef int F_BOOL;

namespace panache {


/*!
 * \brief Implements the libERD backend for calculation
 *        of electron repulsion integrals.
 */
class ERDTwoElectronInt : public TwoBodyAOInt
{


protected:
    
    double *new_cc_1_;  //!< The list of renormalized contraction coefficients for center 1
    double *new_cc_2_;  //!< The list of renormalized contraction coefficients for center 2
    double *new_cc_3_;  //!< The list of renormalized contraction coefficients for center 3
    double *new_cc_4_;  //!< The list of renormalized contraction coefficients for center 4
    
    double *alpha_1_;  //!< All basis set 1 exponents, stored as a flat array
    double *alpha_2_;  //!< All basis set 2 exponents, stored as a flat array
    double *alpha_3_;  //!< All basis set 3 exponents, stored as a flat array
    double *alpha_4_;  //!< All basis set 4 exponents, stored as a flat array
    
    double *xyz_1_;  //!< The x,y,z coordinates for each shell in basis set 1
    double *xyz_2_;  //!< The x,y,z coordinates for each shell in basis set 2
    double *xyz_3_;  //!< The x,y,z coordinates for each shell in basis set 3
    double *xyz_4_;  //!< The x,y,z coordinates for each shell in basis set 4
    
    double *cc_;  //!< The list of contraction coefficients
    
    double *alpha_;  //!< The list of exponents

    
    size_t d_buffer_size_;  //!< The current size of the integral buffer
    size_t i_buffer_size_;  //!< The current size of the integer scratch space
    F_INT *iscratch_;  //!< The integer scratch space
    double *dscratch_;  //!< The double scratch space, which has junk at the start, and integrals at the end

    int *pgto_offsets_1_;  //!< The address of the first contraction coefficient for each shell on center 1
    int *pgto_offsets_2_;  //!< The address of the first contraction coefficient for each shell on center 2
    int *pgto_offsets_3_;  //!< The address of the first contraction coefficient for each shell on center 3
    int *pgto_offsets_4_;  //!< The address of the first contraction coefficient for each shell on center 4
    
    int *npgto_1_;  //!< The number of primitive GTOs per shell in basis set 1
    int *npgto_2_;  //!< The number of primitive GTOs per shell in basis set 2
    int *npgto_3_;  //!< The number of primitive GTOs per shell in basis set 3
    int *npgto_4_;  //!< The number of primitive GTOs per shell in basis set 4
    
    int *am_1_;  //!< The angular momentum of each shell in basis set 1
    int *am_2_;  //!< The angular momentum of each shell in basis set 2
    int *am_3_;  //!< The angular momentum of each shell in basis set 3
    int *am_4_;  //!< The angular momentum of each shell in basis set 4
    
    
    F_INT buffer_offset_;  //!< The start address in the target integral buffer
    F_INT ccbeg_[4];  //!< The first primitive in each contracted function for shells P, Q, R, and S
    F_INT ccend_[4];  //!< The last primitive in each contracted function for shells P, Q, R, and S
    F_BOOL screen_;  //!< Whether to apply screening or not, within ERD
    F_BOOL spheric_;  //!< Whether ERD should use spherical harmonic basis functions
    bool same_bs_;  //!< Not relating to the monotony of integral computations, but whether the basis sets are all the same


    /*!
     * \brief Calculates the maximum amount of scratch space needed
     */ 
    void compute_scratch_size();


public:
    /*!
     * \brief Constructor 
     *
     * \param [in] bs1 The basis set on the 1st center
     * \param [in] bs2 The basis set on the 2nd center
     * \param [in] bs3 The basis set on the 3rd center
     * \param [in] bs4 The basis set on the 4th center
     */
    ERDTwoElectronInt(const SharedBasisSet bs1,
               const SharedBasisSet bs2,
               const SharedBasisSet bs3,
               const SharedBasisSet bs4);

    virtual ~ERDTwoElectronInt();


    // See TwoBodyAOInt::compute_shell
    virtual size_t compute_shell(int sh1, int sh2, int sh3, int sh4);
};

} // close namespace panache

#endif // PANACHE_ERDTWOELECTRONINT_H

