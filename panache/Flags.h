/*! \file
 *  \brief Flags for some parameters passed to ThreeIndexTensor
 *  \author Benjamin Pritchard (ben@bennyp.org)
 */

// mostly done this way (rather than enum) because of the c interface

#ifndef PANACHE_FLAGS_H
#define PANACHE_FLAGS_H

    /*! \name Flags specifing which type of approximation this is */
    ///@{
    #define QTYPE_DFQSO 1  //!< Select DF-based so tensor
    #define QTYPE_CHQSO 2  //!< Select Cholesky-based so tensor
    ///@}

    /*! \name Flags specifing which tensors to calculate */
    ///@{
    #define QGEN_QSO   1  //!< Select Qso
    #define QGEN_QMO   2  //!< Select Qmo
    #define QGEN_QOO   4  //!< Select Qoo
    #define QGEN_QOV   8  //!< Select Qov
    #define QGEN_QVV   16 //!< Select Qvv
    ///@}

    /*! \name Flags specifing how tensors should be stored */
    ///@{
    #define QSTORAGE_PACKED  1     //!< Internal use only
    #define QSTORAGE_ONFLY   2     //!< Generate on-the-fly (not yet implemented)
    #define QSTORAGE_BYQ     4     //!< Store with Q as the first (slowest) index
    #define QSTORAGE_INMEM   8     //!< Store in memory (core)
    #define QSTORAGE_ONDISK  16    //!< Store on disk
    #define QSTORAGE_KEEPDISK  32  //!< Don't erase the file on disk when done. If QSTORAGE_INMEM, writes to disk when done.
    #define QSTORAGE_READDISK  64  //!< Read the file previously saved with QSTORAGE_KEEPDISK

    #ifdef PANACHE_CYCLOPS
    #define QSTORAGE_CYCLOPS 2048  //!< Use Cyclops library
    #endif
    ///@}

    /*! \name Flags specifing basis function ordering */
    ///@{
    #define BSORDER_PSI4   0 //!< Order as Psi4 does
    #define BSORDER_GAMESS 1 //!< Order as GAMESS does
    #define BSORDER_DALTON 2 //!< Order as DALTON does
    ///@}


    /*! \name Flags specifing density fitting metrics */
    ///@{
    #define DFOPT_COULOMB 1 //!< Coulomb metric

    #define DFOPT_EIGINV  2048  //!< Take the inverse sqrt of the metric
    #define DFOPT_CHOINV  4096 //!< Cholesky inverse of the metric
    ///@}

#endif
