/*! \file
 *  \brief Flags for some parameters passed to DFTensor
 *  \author Benjamin Pritchard (ben@bennyp.org)
 */

// mostly done this way (rather than enum) because of the c interface

#ifndef PANACHE_FLAGS_H
#define PANACHE_FLAGS_H


    #define QFLAGS_QMO 1
    #define QFLAGS_QOO 2
    #define QFLAGS_QOV 4
    #define QFLAGS_QVV 8

    #define BSORDER_PSI4   1
    #define BSORDER_GAMESS 2

    #define QSTORAGE_INMEM  1
    #define QSTORAGE_ONDISK 2
    #define QSTORAGE_ONFLY  4
    #define QSTORAGE_BYQ    8
    #define QSTORAGE_PACKED 16 //!< Internal use only


#endif
