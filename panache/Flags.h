/*! \file
 *  \brief Flags for some parameters passed to DFTensor
 *  \author Benjamin Pritchard (ben@bennyp.org)
 */

// mostly done this way (rather than enum) because of the c interface

#ifndef PANACHE_FLAGS_H
#define PANACHE_FLAGS_H

    /*! \name Flags specifing which tensors to calculate */
    ///@{
    #define QGEN_QMO 1 //!< Generate Qmo
    #define QGEN_QOO 2 //!< Generate Qoo
    #define QGEN_QOV 4 //!< Generate Qov
    #define QGEN_QVV 8 //!< Generate Qvv
    ///@}

    /*! \name Flags specifing how tensors should be stored */
    ///@{
    #define QSTORAGE_PACKED 1  //!< Internal use only
    #define QSTORAGE_INMEM  2  //!< Store in memory (core)
    #define QSTORAGE_ONDISK 4  //!< Store on disk
    #define QSTORAGE_ONFLY  8  //!< Generate on-the-fly (not yet implemented)
    #define QSTORAGE_BYQ    16 //!< Store with Q as the first (slowest) index
    ///@}

    /*! \name Flags specifing basis function ordering */
    ///@{
    #define BSORDER_PSI4   1 //!< Order as Psi4 does
    #define BSORDER_GAMESS 2 //!< Order as GAMESS does
    ///@}

#endif
