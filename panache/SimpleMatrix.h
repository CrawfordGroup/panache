/*!
 * \file
 * \brief Simple 2D Matrix
 * \author Benjamin Pritchard (ben@bennyp.org)
 */

#ifndef PANACHE_SIMPLEMATRIX_H
#define PANACHE_SIMPLEMATRIX_H

#include <memory>

namespace panache
{

/*! \brief A simple 2D matrix with continuous storage
 *
 * It doesn't get much more simple than this. Notice that a
 * lot of features are missing - I've only implemented
 * what I really need.
 *
 */
class SimpleMatrix
{
private:

    double * data_;         //!< The actual matrix data
    unsigned int nrow_;     //!< Number of rows in the matrix
    unsigned int ncol_;     //!< Number of columns in the matrix

    /*!
     * \brief Deletes the data in the matrix, checking to see
     *        if data is there first.
     */
    void Delete_(void)
    {
        if(data_ != nullptr)
        {
            delete [] data_;
            data_ = nullptr;
            nrow_ = ncol_ = 0;
        }
    }

public:

    // Remove these functions
    SimpleMatrix(const SimpleMatrix & rhs) = delete;
    SimpleMatrix & operator=(const SimpleMatrix & rhs) = delete;


    /*!
     * \brief Create a matrix
     *
     * Note that this does not zero the matrix contents
     *
     * \param [in] nrow Number of rows in the matrix
     * \param [in] ncol Number of columns in the matrix
     *
     */
    SimpleMatrix(unsigned int nrow, unsigned int ncol)
        : data_(nullptr),nrow_(nrow),ncol_(ncol)
    {
        allocate(nrow,ncol);
    }


    /*!
     * \brief Create an empty matrix
     */
    SimpleMatrix(void)
        : data_(nullptr),nrow_(0),ncol_(0)
    {
    }


    /*!
     * \brief Destructs the matrix and frees memory
     */
    ~SimpleMatrix()
    {
        Delete_();
    }


    /*!
     * \brief Returns the number of rows in the matrix
     */
    unsigned int nrow(void) const
    {
        return nrow_;
    }


    /*!
     * \brief Returns the number of columns in the matrix
     */
    unsigned int ncol(void) const
    {
        return ncol_;
    }


    /*!
     * \brief Access data element [i,j] (by reference)
     *
     * \param [in] i Row index
     * \param [in] j Column index
     */
    double & operator() (unsigned int i, unsigned int j)
    {
        return data_[i*ncol_ + j];
    }


    /*!
     * \brief Access data element [i,j] (by reference, const version)
     *
     * \param [in] i Row index
     * \param [in] j Column index
     */
    const double & operator() (unsigned int i, unsigned int j) const
    {
        return data_[i*ncol_ + j];
    }



    /*!
     * \brief (Re)Allocate memory
     *
     * Currently stored data is deleted first. If the overall
     * size is the same, reallocation does not take place.
     *
     * \param [in] nrow Number of rows in the matrix
     * \param [in] ncol Number of columns in the matrix
     *
     */
    void allocate(unsigned int nrow, unsigned int ncol)
    {

        // allocate only if the size changes
        if((nrow*ncol) != (nrow_ * ncol_))
        {
            Delete_();
            data_ = new double[nrow*ncol];
        }

        nrow_ = nrow;
        ncol_ = ncol;

    }


    /*!
     * \brief Sets all matrix elements to zero
     */
    void zero(void)
    {
        for(size_t i = 0; i < nrow_*ncol_; i++)
            data_[i] = 0.0;
    }


    /*!
     * \brief Returns a raw pointer to the matrix
     *
     * \return A raw pointer to the matrix data.
     */
    double * pointer(void)
    {
        return data_;
    }

};

} // close namespace Panache

#endif

