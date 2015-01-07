/*! \file
 *  \brief Flags for some parameters passed to ThreeIndexTensor
 *  \author Benjamin Pritchard (ben@bennyp.org)
 */

// mostly done this way (rather than enum) because of the c interface

#ifndef PANACHE_FLAGS_H
#define PANACHE_FLAGS_H

    /*! \name Flags specifing which tensors to calculate */
    ///@{
    #define QGEN_DFQSO 1 //!< Select DF-based so tensor
    #define QGEN_CHQSO 2 //!< Select Cholesky-based so tensor
    #define QGEN_QMO   4 //!< Select Qmo
    #define QGEN_QOO   8 //!< Select Qoo
    #define QGEN_QOV   16 //!< Select Qov
    #define QGEN_QVV   32 //!< Select Qvv
    ///@}

    /*! \name Flags specifing how tensors should be stored */
    ///@{
    #define QSTORAGE_PACKED  1  //!< Internal use only
    #define QSTORAGE_ONFLY   2  //!< Generate on-the-fly (not yet implemented)
    #define QSTORAGE_BYQ     4 //!< Store with Q as the first (slowest) index
    #define QSTORAGE_INMEM   8  //!< Store in memory (core)
    #define QSTORAGE_ONDISK  16  //!< Store on disk

    #ifdef PANACHE_CYCLOPS
    #define QSTORAGE_CYCLOPS 32  //!< Use Cyclops library
    #endif
    ///@}

    /*! \name Flags specifing basis function ordering */
    ///@{
    #define BSORDER_PSI4   0 //!< Order as Psi4 does
    #define BSORDER_GAMESS 1 //!< Order as GAMESS does
    #define BSORDER_DALTON 2 //!< Order as DALTON does
    ///@}

#endif
