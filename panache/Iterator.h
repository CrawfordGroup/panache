/*! \file
 * \brief Iterator over combined orbital indices
 * \author Benjamin Pritchard (ben@bennyp.org)
 */

#ifndef PANACHE_ITERATOR_H
#define PANACHE_ITERATOR_H

namespace panache {

template<typename ITTYPE>
class IndexIterator
{
private:
    ITTYPE data_;


protected:
    const ITTYPE & Iterator(void) const { return data_; }


public:
    IndexIterator(const ITTYPE & it) : data_(it) { }



    
    /*!
     * \brief Get the "master" index
     */ 
    int Index(void) const { return data_.Index(); }

    /*!
     * \brief Get whether or not these are packed indices
     */ 
    bool Packed(void) const { return data_.Packed(); }

    /*!
     * \brief Get the permutational symmetry factor
     */ 
    bool Perm(void) const { return data_.Perm(); }

    bool Valid(void) const { return data_.Valid(); }

    /*!
     * \brief Test if the iterator is valid or not
     */
    operator bool() const { return data_.Valid(); }

    /*!
     * \brief Prefix increment operator
     *
     * Move to the next combined index
     */
    IndexIterator & operator++()
    {
        return (*this += 1);
    }

    /*!
     * \brief Prefix decrement operator
     *
     * Move to the previous combined index
     */
    IndexIterator & operator--()
    {
        return (*this -= 1);
    }


    /*!
     * \brief Postfix increment operator
     *
     * Move to the next combined index
     */
    IndexIterator operator++(int)
    {
        IndexIterator n(*this);
        ++(*this);
        return n; 
    }

    /*!
     * \brief Postfix decrement operator
     *
     * Move to the previous combined index
     */
    IndexIterator operator--(int)
    {
        IndexIterator n(*this);
        --(*this);
        return n; 
    }


    /*!
     * \brief Move forward
     *
     * Jump ahead some number of combined indices
     */
    IndexIterator operator+(int ij)
    {
        IndexIterator n(*this);
        n += ij;
        return n;
    }

    /*!
     * \brief Move backward
     *
     * Jump backward some number of combined indices
     */
    IndexIterator operator-(int ij)
    {
        IndexIterator n(*this);
        n -= ij;
        return n;
    }

    IndexIterator & operator+=(int a)
    {
        if(a < 0)
            data_ -= (-1*a);
        else
            data_ += a;

        return *this;
    }   

    IndexIterator & operator-=(int a)
    {
        if(a < 0)
            data_ += (-1*a);
        else
            data_ -= a;
    }   

};

struct IJIteratorType
{
    int i;       //!< Current i value
    int j;       //!< Current j value
    int ij;      //!< Current combined index
    int ni;      //!< Dimension along i
    int nj;      //!< Dimension along j
    int nij;     //!< Total number of combined indices
    bool packed;

    IJIteratorType(int ni, int nj, bool packed)
        : ni(ni),nj(nj),packed(packed)
    {
        nij = (packed ? (ni*(ni+1))/2 : ni*nj);
        i = j = ij = 0;
    }

    IJIteratorType(const IJIteratorType & rhs) = default;


    int Index(void) const
    {
        return ij; 
    }

    bool Valid(void) const
    {
        return (ij < nij && ij >= 0);
    }

    short Perm(void) const
    {
        if(packed)
            return (i == j ? 1 : 2);
        else
            return 1; 
    }

    bool Packed(void) const
    {
        return packed;
    }


    /*!
     * \brief Move forward
     *
     * Jump ahead some number of combined indices
     */
    IJIteratorType & operator+=(int a)
    {
        if(!Valid())
            return *this;

        if(packed)
        {
            while(a > 0)
            {
                j++;
                if(j > i)
                {
                    i++;
                    j = 0;
                }
                a--;
                ij++;
            }
        }
        else
        {
            i += (j + a)/nj;
            j = (j + a) % nj;
            ij += a;
        }

        return *this;
    }


    /*!
     * \brief Move backward
     *
     * Jump back some number of combined indices
     */
    IJIteratorType & operator-=(int a)
    {
        if(!Valid())
            return *this;

        if(packed)
        {
            while(a > 0)
            {
                ij--;
                j--;
                if(j < 0)
                {
                    i--;
                    j = i;

                    if(i < 0)
                        break;
                }
                a--;
            }
        }
        else
        {
            int idec = a / nj;
            int jdec = a % nj;

            if(jdec > j)
            {
                idec++;
                j = nj;
            }

            i -= idec;
            j -= jdec;
            ij -= a;
        }

        return *this;
    }
};


struct QIteratorType
{
    int q;       //!< Current q value
    int nq;      //!< Number of aux functions
    bool packed;

    QIteratorType(int nq, bool packed)
        : nq(nq),packed(packed)
    {
        q = 0;
    }

    QIteratorType(const QIteratorType & rhs) = default;
    
    int Index(void) const
    {
        return q; 
    }

    bool Valid(void) const
    {
        return (q >= 0 && q < nq);
    }

    short Perm(void) const
    {
        return 1;
    }

    bool Packed(void) const
    {
        return packed;
    }

    /*!
     * \brief Move forward
     *
     * Jump ahead some number of combined indices
     */
    QIteratorType & operator+=(int a)
    {
        q += a;
        return *this;
    }


    /*!
     * \brief Move backward
     *
     * Jump back some number of combined indices
     */
    QIteratorType & operator-=(int a)
    {
        q -= a;
        return *this;
    }
};

class IJIterator : public IndexIterator<IJIteratorType>
{
public:
    typedef std::function<int(const IJIterator &)> GetBatchFunc;

    IJIterator(const IJIteratorType ijt) : IndexIterator(ijt) { }
    IJIterator(int ni, int nj, bool packed) : IndexIterator(IJIteratorType(ni, nj, packed)) { }

    int i(void) const { return Iterator().i; } 
    int j(void) const { return Iterator().j; } 
    int ij(void) const { return Iterator().ij; } 
};


class QIterator : public IndexIterator<QIteratorType>
{
public:
    typedef std::function<int(const QIterator &)> GetBatchFunc;

    QIterator(const QIteratorType & qi) : IndexIterator(qi) { }
    QIterator(int nq, bool packed) : IndexIterator(QIteratorType(nq, packed)) { }

    int q(void) const { return Iterator().q; } 
};


} // close namespace panache

#endif


