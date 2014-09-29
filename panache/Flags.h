/*! \file
 *  \brief Flags for some parameters passed to ThreeIndexTensor
 *  \author Benjamin Pritchard (ben@bennyp.org)
 */

// mostly done this way (rather than enum) because of the c interface

#ifndef PANACHE_FLAGS_H
#define PANACHE_FLAGS_H

    /*! \name Flags specifing which tensors to calculate */
    ///@{
    #define QGEN_QSO 1 //!< Generate Qso
    #define QGEN_QMO 2 //!< Generate Qmo
    #define QGEN_QOO 4 //!< Generate Qoo
    #define QGEN_QOV 8 //!< Generate Qov
    #define QGEN_QVV 16 //!< Generate Qvv
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
    #define BSORDER_PSI4   1 //!< Order as Psi4 does
    #define BSORDER_GAMESS 2 //!< Order as GAMESS does
    ///@}

#endif
