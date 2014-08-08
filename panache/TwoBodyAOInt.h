/*! \file
 * \brief Base class for two-body AO integrals (header)
 * \author Benjamin Pritchard (ben@bennyp.org)
 */

#ifndef PANACHE_TWOBODYAOINT_H
#define PANACHE_TWOBODYAOINT_H

#include <memory>

namespace panache {

class IntegralFactory;
class AOShellCombinationsIterator;
class BasisSet;
typedef std::shared_ptr<BasisSet> SharedBasisSet;
class GaussianShell;


/*!
 * \brief Base class for two-body AO integrals
 */
class TwoBodyAOInt
{
protected:

    SharedBasisSet bs1_;  //!< Basis set on center 1
    SharedBasisSet bs2_;  //!< Basis set on center 2
    SharedBasisSet bs3_;  //!< Basis set on center 3
    SharedBasisSet bs4_;  //!< Basis set on center 4

    const SharedBasisSet original_bs1_;  //!< Original basis set on center 1 (before any reordering)
    const SharedBasisSet original_bs2_;  //!< Original basis set on center 2 (before any reordering)
    const SharedBasisSet original_bs3_;  //!< Original basis set on center 3 (before any reordering)
    const SharedBasisSet original_bs4_;  //!< Original basis set on center 4 (before any reordering)

    
    double *target_;            //!< Buffer to hold the final integrals.
    int curr_buff_size_;        //!< Number of integrals in the current buffer
    double *tformbuf_;          //!< Buffer to hold the transformation intermediates.
    double *source_;            //!< Buffer to hold the initially computed integrals.
    int max_unique_quartets_;   //!< Maximum number of unique quartets needed to compute a set of SO's
    int natom_;                 //!< Number of atoms.
    bool force_cartesian_;      //!< Whether to force integrals to be generated in the Cartesian (AO) basis;
    unsigned char buffer_offsets_[4];  //!< The order of the derivative integral buffers, after permuting shells



    /*!
     * \brief Copies integrals, permuting if necessary
     *
     * \p s should contain the integrals over shells (sh1 sh2 | sh3 sh4)
     *
     * \param [in] s The source to be copied
     * \param [in] t The destination
     * \param [in] sh1 Shell 1 whose integrals are in \p s
     * \param [in] sh2 Shell 2 whose integrals are in \p s
     * \param [in] sh3 Shell 3 whose integrals are in \p s
     * \param [in] sh4 Shell 4 whose integrals are in \p s
     * \param [in] p12 If shells 1 and 2 should be permuted
     * \param [in] p34 If shells 3 and 4 should be permuted
     * \param [in] p13p24 If shells 1 & 3 as well as 2 & 4 should be permuted
     */
    void permute_target(double *s, double *t, int sh1, int sh2, int sh3, int sh4, bool p12, bool p34, bool p13p24);



    /*!
     * \brief Copies integrals, changing the order from 1234 to 1243
     *
     * The size of \p s should be \p nbf1 * \p nbf2 * \p nbf3 * \p nbf4
     *
     * \param [in] s The source to be copied
     * \param [in] t The destination
     * \param [in] nbf1 The number of basis functions in the first shell
     * \param [in] nbf2 The number of basis functions in the second shell
     * \param [in] nbf3 The number of basis functions in the third shell
     * \param [in] nbf4 The number of basis functions in the fourth shell
     */
    void permute_1234_to_1243(double *s, double *t, int nbf1, int nbf2, int nbf3, int nbf4);



    /*!
     * \brief Copies integrals, changing the order from 1234 to 2134
     *
     * The size of \p s should be \p nbf1 * \p nbf2 * \p nbf3 * \p nbf4
     *
     * \param [in] s The source to be copied
     * \param [in] t The destination
     * \param [in] nbf1 The number of basis functions in the first shell
     * \param [in] nbf2 The number of basis functions in the second shell
     * \param [in] nbf3 The number of basis functions in the third shell
     * \param [in] nbf4 The number of basis functions in the fourth shell
     */
    void permute_1234_to_2134(double *s, double *t, int nbf1, int nbf2, int nbf3, int nbf4);



    /*!
     * \brief Copies integrals, changing the order from 1234 to 2143
     *
     * The size of \p s should be \p nbf1 * \p nbf2 * \p nbf3 * \p nbf4
     *
     * \param [in] s The source to be copied
     * \param [in] t The destination
     * \param [in] nbf1 The number of basis functions in the first shell
     * \param [in] nbf2 The number of basis functions in the second shell
     * \param [in] nbf3 The number of basis functions in the third shell
     * \param [in] nbf4 The number of basis functions in the fourth shell
     */
    void permute_1234_to_2143(double *s, double *t, int nbf1, int nbf2, int nbf3, int nbf4);



    /*!
     * \brief Copies integrals, changing the order from 1234 to 3412
     *
     * The size of \p s should be \p nbf1 * \p nbf2 * \p nbf3 * \p nbf4
     *
     * \param [in] s The source to be copied
     * \param [in] t The destination
     * \param [in] nbf1 The number of basis functions in the first shell
     * \param [in] nbf2 The number of basis functions in the second shell
     * \param [in] nbf3 The number of basis functions in the third shell
     * \param [in] nbf4 The number of basis functions in the fourth shell
     */
    void permute_1234_to_3412(double *s, double *t, int nbf1, int nbf2, int nbf3, int nbf4);



    /*!
     * \brief Copies integrals, changing the order from 1234 to 4312
     *
     * The size of \p s should be \p nbf1 * \p nbf2 * \p nbf3 * \p nbf4
     *
     * \param [in] s The source to be copied
     * \param [in] t The destination
     * \param [in] nbf1 The number of basis functions in the first shell
     * \param [in] nbf2 The number of basis functions in the second shell
     * \param [in] nbf3 The number of basis functions in the third shell
     * \param [in] nbf4 The number of basis functions in the fourth shell
     */
    void permute_1234_to_4312(double *s, double *t, int nbf1, int nbf2, int nbf3, int nbf4);



    /*!
     * \brief Copies integrals, changing the order from 1234 to 3421
     *
     * The size of \p s should be \p nbf1 * \p nbf2 * \p nbf3 * \p nbf4
     *
     * \param [in] s The source to be copied
     * \param [in] t The destination
     * \param [in] nbf1 The number of basis functions in the first shell
     * \param [in] nbf2 The number of basis functions in the second shell
     * \param [in] nbf3 The number of basis functions in the third shell
     * \param [in] nbf4 The number of basis functions in the fourth shell
     */
    void permute_1234_to_3421(double *s, double *t, int nbf1, int nbf2, int nbf3, int nbf4);



    /*!
     * \brief Copies integrals, changing the order from 1234 to 4321
     *
     * The size of \p s should be \p nbf1 * \p nbf2 * \p nbf3 * \p nbf4
     *
     * \param [in] s The source to be copied
     * \param [in] t The destination
     * \param [in] nbf1 The number of basis functions in the first shell
     * \param [in] nbf2 The number of basis functions in the second shell
     * \param [in] nbf3 The number of basis functions in the third shell
     * \param [in] nbf4 The number of basis functions in the fourth shell
     */
    void permute_1234_to_4321(double *s, double *t, int nbf1, int nbf2, int nbf3, int nbf4);

    /*!
     * \brief Constructor (to be called only by derived classes) 
     *
     * \param [in] original_bs1 The basis set on the 1st center
     * \param [in] original_bs2 The basis set on the 2nd center
     * \param [in] original_bs3 The basis set on the 3rd center
     * \param [in] original_bs4 The basis set on the 4th center
     */
    TwoBodyAOInt(const SharedBasisSet original_bs1, 
                 const SharedBasisSet original_bs2,
                 const SharedBasisSet original_bs3,
                 const SharedBasisSet original_bs4);

public:
    TwoBodyAOInt(const TwoBodyAOInt & t) = delete;
    TwoBodyAOInt(const TwoBodyAOInt && t) = delete;
    TwoBodyAOInt & operator=(const TwoBodyAOInt && t) = delete;



    virtual ~TwoBodyAOInt();

    /*!
     * \brief Basis set on center one
     */
    SharedBasisSet basis() { return original_bs1_; }



    /*!
     * \brief Basis set on center one
     */
    SharedBasisSet basis1() { return original_bs1_; }



    /*!
     * \brief Basis set on center two
     */
    SharedBasisSet basis2() { return original_bs2_; }



    /*!
     * \brief Basis set on center three
     */
    SharedBasisSet basis3() { return original_bs3_; }



    /*!
     * \brief Basis set on center four
     */
    SharedBasisSet basis4() { return original_bs4_; }




    /*!
     * \brief Sets whether we're forcing this object to always generate Cartesian integrals
     * 
     * \param [in] t_f Set to true to force cartesian, false otherwise
     */
    void set_force_cartesian(bool t_f) { force_cartesian_ = t_f; }




    /*!
     * \brief Get the buffer where the integrals are placed
     */
    const double *buffer() const { return target_; }




    /*!
     * \brief Compute integrals between 4 shells. Result are obtained via buffer()
     *
     * \param [in] shellIter Iterator representing the shells to compute
     * \return Number of computed integrals
     */
    virtual size_t compute_shell(const AOShellCombinationsIterator & shellIter) = 0;




    /*!
     * \brief Compute the integrals between 4 shells.
     *
     * \param [in] sh1 First shell
     * \param [in] sh2 Second shell
     * \param [in] sh3 Third shell
     * \param [in] sh4 Fourth shell
     * \return Number of computed integrals
     */
    virtual size_t compute_shell(int sh1, int sh2, int sh3, int sh4) = 0;


    /*!
     * \brief Transform cartesian to pure spherical harmonics
     *
     * Results go back to buffer_.
     *
     * \param [in] sh1 First shell
     * \param [in] sh2 Second shell
     * \param [in] sh3 Third shell
     * \param [in] sh4 Fourth shell
     * \param [in] nchunk Chunk size (???)
     */
    void pure_transform(int sh1, int sh2, int sh3, int sh4, int nchunk);
};


/*!
 * \brief A shared TwoBodyAOInt object
 */
typedef std::shared_ptr<TwoBodyAOInt> SharedTwoBodyAOInt;

}

#endif //PANACHE_TWOBODYAOINT_H
