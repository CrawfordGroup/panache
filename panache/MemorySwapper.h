#ifndef PANACHE_MEMORYSWAPPER_H
#define PANACHE_MEMORYSWAPPER_H

#include <algorithm>

namespace panache {
namespace reorder {


class MemorySwapper
{
public:
    virtual void swap(double *a, double *b) = 0;
};


class TotalMemorySwapper : public MemorySwapper
{
public:
    TotalMemorySwapper(size_t size) 
           : _size(size), _buff(new double[_size])
    {}

    ~TotalMemorySwapper()
    {
        delete [] _buff;
    }

    void swap(double *a, double *b)
    {
        using std::copy;

        copy(a, a+_size, _buff);
        copy(b, b+_size, a);
        copy(_buff, _buff+_size, b);
    }

private:
    size_t _size;
    double * _buff;
};



class LimitedMemorySwapper : public MemorySwapper
{
public:
    LimitedMemorySwapper(size_t size, size_t buffsize) 
                      : _size(size),
                        _buffsize(size <= buffsize ? size : buffsize),
                        _buff(new double[_buffsize])
    { }

    ~LimitedMemorySwapper()
    {
        delete [] _buff;
    }

    void swap(double *a, double *b)
    {
        using std::copy;

        size_t left = _size;
        size_t tocopy = 0;

        do
        {
            tocopy = (left <= _buffsize) ? left : _buffsize;
            copy(a, a+tocopy, _buff);
            copy(b, b+tocopy, a);
            copy(_buff, _buff+tocopy, b);
            left -= tocopy;
            a += tocopy;
            b += tocopy;

        } while(left > 0);
    }

private:
    size_t _size, _buffsize;
    double * _buff;
};

} } // close namespace panache::reorder

#endif
