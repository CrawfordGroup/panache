/*! \file
 *  \brief Parser for basis set files (header)
 *  \author Benjamin Pritchard (ben@bennyp.org)
 */

#ifndef PANACHE_BASISSETPARSER_H
#define PANACHE_BASISSETPARSER_H

#include <vector>
#include <string>
#include "panache/Exception.h"
#include "panache/GaussianShell.h"

namespace panache {

class BasisSet;

/*!
 * \brief Abstract class for parsing basis sets from a text file.
 *
 */
class BasisSetParser
{

protected:
    std::string filename_; //!< Full path to file to read


public:

    bool force_puream_or_cartesian_;   //!< If the parser needs to force spherical or cartesian
    bool forced_is_puream_;            //!< Is the forced value to use puream?  (Otherwise force Cartesian if \p force_puream_or_cartesian_ is true).

    /*!
     * \brief Constructs a BasisSetParser object
     */ 
    BasisSetParser();


    /*!
     * \brief Constructs a BasisSetParser object
     *
     * \param [in] forced_puream Force spherical basis functions 
     */
    BasisSetParser(bool forced_puream);



    /*!
     * \brief Destructor
     */ 
    virtual ~BasisSetParser();


    /*!
     *  \brief Reads a file into an array of strings.
     *  \param [in] filename Full path of the file to load
     *  \return Contents of the file as a vector of strings.
     */
    std::vector<std::string> load_file(const std::string& filename);


    /*! 
     * \brief Take a multiline string and convert it to a vector of strings
     *
     * \param [in] data A multi-line string to split
     * \return Vector of strings containing all the lines from \p data
     *
     */
    std::vector<std::string> string_to_vector(const std::string& data);



    /*!
     * \brief Given a file (as a vector of strings), parse for the needed for atom.
     *
     * \param [in] symbol Atomic/Center symbol to look up
     * \param [in] dataset Basis set information to look through 
     * \return Shell information for this symbol
     */
    virtual std::vector<ShellInfo> parse(const std::string& symbol, const std::vector<std::string>& dataset) = 0;

};


/*! 
 *   \brief Class for reading in basis sets formatted for Gaussian.
 */
class Gaussian94BasisSetParser : public BasisSetParser
{
public:

    /*!
     * \brief Default constructor 
     */ 
    Gaussian94BasisSetParser(): BasisSetParser() {}

    /*!
     * \brief Constructor with option to force spherical basis functions
     *
     * \param [in] forced_puream Set to true to force spherical basis functions
     */ 
    Gaussian94BasisSetParser(bool forced_puream): BasisSetParser(forced_puream) {}

    /*!
     * \brief Given a file (as a vector of strings), parse for the needed for atom.
     *
     * \param [in] symbol Atomic/Center symbol to look up
     * \param [in] dataset Basis set information to look through 
     * \return Shell information for this symbol
     */
    virtual std::vector<ShellInfo> parse(const std::string& symbol, const std::vector<std::string>& dataset);
};


/*!
 * \brief Shared pointer to a basis set parser
 */
typedef std::shared_ptr<BasisSetParser> SharedBasisSetParser;


} // end namespace panache 

#endif
