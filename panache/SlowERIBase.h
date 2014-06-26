/*! \file
 * \brief Class for the calculation of slow but accurate ERI (header)
 * \author Benjamin Pritchard (ben@bennyp.org)
 */

#ifndef PANACHE_SLOWERIBASE_H
#define PANACHE_SLOWERIBASE_H

#define MAXFAC 100
#define EPS 1.0E-17


class SlowERIBase
{
private:
    double * df;   //!< Holds a lookup of double factorials
    double * fac;  //!< Holds a lookup of factorials
    double ** bc;  //!< Holds a lookup of binomial combinations

    /*!
     * \brief Calculates the fundamental
     *
     * This function computes infamous integral Fn(t). For its definition 
     * see Obara and Saika paper, or Shavitt's chapter in the 
     * Methods in Computational Physics book (see reference below). 
     * This piece of code is from Dr. Justin Fermann's program CINTS
     *
     * \todo Document \p n and \p t
     *
     * \param [in] F Fundamentals are placed in this buffer
     * \param [in] n ???
     * \param [in] t ???
     */
    void calc_f(double *F, int n, double t);

    /*!
     * \brief Claculates the normalization constant of a primitive
     * \f$x^l y^m z^n \exp^(-\alpha r^2) \f$
     *
     * \param [in] l1 Exponent on x
     * \param [in] m1 Exponent on y
     * \param [in] n1 Exponent on z
     * \param [in] alpha1 Exponent of the gaussian function
     * \param [in] A Array of 3 doubles representing the cartesian coordinates
     * \return Normalization constant for the given primitive
     */
    double norm_const(unsigned int l1, unsigned int m1, unsigned int n1,
                      double alpha1, const double* A);

    /*!
     * \brief Creates an array and initializes all elements to zero
     *
     * \param [in] size Length of the array
     * \return Pointer to a newly-allocated array
     */
    double* init_array(unsigned long int size);


    /*!
     * \brief Frees an array created by init_array
     *
     * \param [in] array Pointer to the array to free
     */ 
    void free_array(double* array);


public:
    /*!
     * \brief Constructor
     *
     * Initializes some private members
     */ 
    SlowERIBase();

    
    /*!
     * \brief Destructor
     *
     * Frees associated memory
     */ 
    ~SlowERIBase();


    /*!
     * \brief Calculate an electron repulsion integral
     *
     * Stable, but slow. Each (cartesian) basis function primitive
     * is given as \f$x^l y^m z^n \exp^(-\alpha r^2) \f$
     *
     * \param [in] l1 First center, exponent on x 
     * \param [in] m1 First center, exponent on y
     * \param [in] n1 First center, exponent on z
     * \param [in] alpha1 First center, exponent of gaussian
     * \param [in] A First center, array of three doubles representing the cartesian coordinates
     *
     * \param [in] l2 Second center, exponent on x 
     * \param [in] m2 Second center, exponent on y
     * \param [in] n2 Second center, exponent on z
     * \param [in] alpha2 Second center, exponent of gaussian
     * \param [in] B Second center, array of three doubles representing the cartesian coordinates
     *
     * \param [in] l3 Third center, exponent on x 
     * \param [in] m3 Third center, exponent on y
     * \param [in] n3 Third center, exponent on z
     * \param [in] alpha3 Third center, exponent of gaussian
     * \param [in] C Third center, array of three doubles representing the cartesian coordinates
     *
     * \param [in] l4 Fourth center, exponent on x 
     * \param [in] m4 Fourth center, exponent on y
     * \param [in] n4 Fourth center, exponent on z
     * \param [in] alpha4 Fourth center, exponent of gaussian
     * \param [in] D Fourth center, array of three doubles representing the cartesian coordinates
     *
     * \param [in] norm_flag If nonzero, normalization will be performed
     *
     * \return The ERI for the given four centers
     */
    double eri(int l1, int m1, int n1, double alpha1, const double* A,
               int l2, int m2, int n2, double alpha2, const double* B,
               int l3, int m3, int n3, double alpha3, const double* C,
               int l4, int m4, int n4, double alpha4, const double* D,
               int norm_flag);

};

#endif

