/*! \file
 * \brief Some basis math functions (source)
 * \author Benjamin Pritchard (ben@bennyp.org)
 */

#include <vector>
#include "panache/Math.h"

// lookup tables in an anonymous namespace
namespace {
  std::vector<double> factorial_table_;
  std::vector<double> double_factorial_table_;
  std::vector<int64_t> packed_row_table_;
  bool init_ = false;
}

namespace panache {
namespace math {

void initmath(void)
{
    if(!init_)
    {
        factorial_table_.resize(FACTORIAL_TABLE_SIZE);
        double_factorial_table_.resize(DBL_FACTORIAL_TABLE_SIZE);
        packed_row_table_.resize(PACKED_ROW_TABLE_SIZE);
        
        factorial_table_[0] = 1.0;
        factorial_table_[1] = 1.0;
        for(int i = 2; i < FACTORIAL_TABLE_SIZE; i++)
            factorial_table_[i] = factorial_table_[i-1]*static_cast<double>(i);
    
        double_factorial_table_[0] = 1.0;
        double_factorial_table_[1] = 1.0;
        for(int i = 2; i < DBL_FACTORIAL_TABLE_SIZE; i++)
            double_factorial_table_[i] = double_factorial_table_[i-2]*static_cast<double>(i);

        // these are offset by one    
        packed_row_table_[0] = 1;
        for(int i = 1; i < PACKED_ROW_TABLE_SIZE; i++)
            packed_row_table_[i] = packed_row_table_[i-1] + (i+1);
    
        init_ = true;
    }
}


double factorial(int n)
{
    if(!init_)
        initmath();

    double r = 1.0;

    while(n >= FACTORIAL_TABLE_SIZE)
    {
        r *= static_cast<double>(n);
        n--;
    }
    return r*factorial_table_[n];
}

double double_factorial(int n)
{
    if(!init_)
        initmath();

    double r = 1.0;
    while(n >= DBL_FACTORIAL_TABLE_SIZE)
    {
        r *= static_cast<double>(n);
        n -= 2;
    }
    return r*double_factorial_table_[n];
}

double double_factorial_nminus1(int n)
{
    if(n == 0)
        return 1.0;
    else
        return double_factorial(n-1);
}

double combinations(int n, int k)
{
    if(!init_)
        initmath();

   double comb;

   if (n == k) return (1.0) ;
   else if (k > n) return(0.0) ;
   else if (k == 0) return(1.0) ;
   comb = factorial(n) / (factorial(k) * factorial(n-k)) ;

   return(comb) ;
}

std::pair<int64_t, int64_t> decomposeij_packed(int64_t ij)
{
    if(!init_)
        initmath();

    if(ij < packed_row_table_.back())
    {
        int64_t i = 0;
        while(ij >= packed_row_table_[i])
            i++;

        return std::pair<int64_t, int64_t>(i, ij-packed_row_table_[i]+i+1);
    }
    else
    {
        // start at end of table
        int64_t i = PACKED_ROW_TABLE_SIZE;
        int64_t ival = packed_row_table_.back() + (i+1);

        // keep going
        while(ij >= ival)
            ival += (++i+1); // hmmmmm

        return std::pair<int64_t, int64_t>(i, ij-ival+i+1);
    }
}

} } //close namespace panache::Math
