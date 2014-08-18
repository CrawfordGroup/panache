/*! \file
 * \brief Iterator over combined orbital indices
 * \author Benjamin Pritchard (ben@bennyp.org)
 */

#ifndef PANACHE_IJITERATOR_H
#define PANACHE_IJITERATOR_H

namespace panache {

/*!
 * \brief Iterator over combined orbital indices
 *
 * Can handle packed or unpacked
 */
class IJIterator
{
private:
    int i_;       //!< Current i value
    int j_;       //!< Current j value
    int ij_;      //!< Current combined index
    int ni_;      //!< Dimension along i
    int nj_;      //!< Dimension along j
    int nij_;     //!< Total number of combined indices
    bool valid_;  //!< Is currently valid or not
    bool packed_; //!< Use packed indices
    int perm_;    //!< Permutational symmetry factor

    /*!
     * \brief Fill in valid_ member
     */
    void validate_(void)
    {
        valid_ = (ij_ < nij_ && ij_ >= 0);
    }

public:

    /*!
     * \brief Constructor
     *
     * \param [in] ni Dimension along first index
     * \param [in] nj Dimension along second index
     * \param [in] packed Use packed indices
     */
    IJIterator(int ni, int nj, bool packed = false) 
             : ni_(ni),nj_(nj),nij_(packed ? (ni_*(ni_+1))/2 : ni_*nj_),packed_(packed)
    { 
        i_ = j_ = ij_  = 0; 
        valid_ = true;
    }

    /*!
     * \brief Copy Constructor
     */ 
    IJIterator(IJIterator & iji)
    {
        i_ = iji.i_;
        j_ = iji.j_;
        ij_ = iji.ij_;
        ni_ = iji.ni_;
        nj_ = iji.nj_;
        valid_ = iji.valid_;
        packed_ = packed_;
    }

    /*!
     * \brief Get the current value of i
     */ 
    int i(void) const { return i_; }

    /*!
     * \brief Get the current value of j
     */ 
    int j(void) const { return j_; }

    /*!
     * \brief Get the current combined index
     */ 
    int ij(void) const { return ij_; }

    /*!
     * \brief Get whether or not these are packed indices
     */ 
    bool packed(void) const { return packed_; }

    /*!
     * \brief Get the permutational symmetry factor
     */ 
    bool perm(void) const { return perm_; }

    /*!
     * \brief Test if the iterator is valid or not
     */
    operator bool() const { return valid_; }

    /*!
     * \brief Prefix increment operator
     *
     * Move to the next combined index
     */
    IJIterator & operator++()
    {
        return (*this += 1);
    }

    /*!
     * \brief Prefix decrement operator
     *
     * Move to the previous combined index
     */
    IJIterator & operator--()
    {
        return (*this -= 1);
    }


    /*!
     * \brief Postfix increment operator
     *
     * Move to the next combined index
     */
    IJIterator operator++(int)
    {
        IJIterator n(*this);
        ++(*this);
        return n; 
    }

    /*!
     * \brief Postfix decrement operator
     *
     * Move to the previous combined index
     */
    IJIterator operator--(int)
    {
        IJIterator n(*this);
        --(*this);
        return n; 
    }


    /*!
     * \brief Move forward
     *
     * Jump ahead some number of combined indices
     */
    IJIterator & operator+=(int ij)
    {
        if(!valid_)
            return *this;

        if(packed_)
        {
            while(ij > 0)
            {
                j_++;
                if(j_ > i_)
                {
                    i_++;
                    j_ = 0;
                }
                ij--;
                ij_++;
            }
        }
        else
        {
            i_ += (j_ + ij)/nj_;
            j_ = (j_ + ij) % nj_;
            ij_ += ij;
        }

        validate_();

        // Set the permutational symmetry factor
        if(valid_ && packed_ && i_ != j_)
            perm_ = 2;
        else
            perm_ = 1;  // will be set to 1 if !valid_, but shouldn't care

        return *this;
    }

    /*!
     * \brief Move backward
     *
     * Jump back some number of combined indices
     */
    IJIterator & operator-=(int ij)
    {
        if(!valid_)
            return *this;

        if(packed_)
        {
            while(ij > 0)
            {
                ij_--;
                j_--;
                if(j_ < 0)
                {
                    i_--;
                    j_ = i_;

                    if(i_ < 0)
                        break;
                }
                ij--;
            }
        }
        else
        {
            int idec = ij / nj_;
            int jdec = ij % nj_;

            if(jdec > j_)
            {
                idec++;
                j_ = nj_;
            }

            i_ -= idec;
            j_ -= jdec;
            ij_ -= ij;
        }

        validate_();

        // Set the permutational symmetry factor
        if(valid_ && packed_ && i_ != j_)
            perm_ = 2;
        else
            perm_ = 1;  // will be set to 1 if !valid_, but shouldn't care

        return *this;
    }


    /*!
     * \brief Move forward
     *
     * Jump ahead some number of combined indices
     */
    IJIterator operator+(int ij)
    {
        IJIterator n(*this);
        n += ij;
        return n;
    }

    /*!
     * \brief Move backward
     *
     * Jump backward some number of combined indices
     */
    IJIterator operator-(int ij)
    {
        IJIterator n(*this);
        n -= ij;
        return n;
    }
};

} // close namespace panache

#endif


