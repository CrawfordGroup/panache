
namespace panache {
namespace math {


//! \todo lookup table
double factorial(int n)
{
    double r = 1.0;
    while(n > 1)
    {
        r *= n;
        n--;
    }
    return r;
}

double double_factorial(int n)
{
    double r = 1.0;
    while(n > 1)
    {
        r *= n;
        n -= 2;
    }
    return r;
}

double double_factorial_nminus1(int n)
{
    double r = 1.0;
    n--;
    while(n > 1)
    {
        r *= n;
        n -= 2;
    }
    return r;
}

double combinations(int n, int k)
{
   double comb;

   if (n == k) return (1.0) ;
   else if (k > n) return(0.0) ;
   else if (k == 0) return(1.0) ;
   comb = factorial(n) / (factorial(k) * factorial(n-k)) ;

   return(comb) ;
}

}
} //close namespace panache::Math
