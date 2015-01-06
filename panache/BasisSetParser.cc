/*! \file
 *  \brief Parser for basis set files (source)
 *  \author Benjamin Pritchard (ben@bennyp.org)
 */


#include <fstream>
#include <algorithm>
#include <sstream>

#include "panache/ShellInfo.h"
#include "panache/GaussianShell.h"
#include "panache/BasisSet.h"
#include "panache/BasisSetParser.h"

// the third parameter of from_string() should be
// one of std::hex, std::dec or std::oct
template <class T>
static bool from_string(T& t,
                        const std::string& s,
                        std::ios_base& (*f)(std::ios_base&))
{
    std::istringstream iss(s);
    return !(iss >> f >> t).fail();
}

static std::string to_lower_copy(std::string str)
{
    std::transform(str.begin(), str.end(), str.begin(), ::tolower);
    return str;
}

// shamelessly stolen from stack overflow
// trim from start
static inline std::string ltrim(std::string s)
{
    s.erase(s.begin(), std::find_if(s.begin(), s.end(), std::not1(std::ptr_fun<int, int>(std::isspace))));
    return s;
}

// trim from end
static inline std::string rtrim(std::string s)
{
    s.erase(std::find_if(s.rbegin(), s.rend(), std::not1(std::ptr_fun<int, int>(std::isspace))).base(), s.end());
    return s;
}

// trim from both ends
static inline std::string trim(std::string s)
{
    return ltrim(rtrim(s));
}



namespace panache
{

BasisSetParser::BasisSetParser()
{
    force_puream_or_cartesian_ = false;
    forced_is_puream_ = false;
}

BasisSetParser::BasisSetParser(bool forced_puream)
{
    force_puream_or_cartesian_ = true;
    forced_is_puream_ = forced_puream;
}

BasisSetParser::~BasisSetParser()
{
}

std::vector<std::string> BasisSetParser::load_file(const std::string& filename)
{
    filename_ = filename;

    // Loads an entire file.
    std::vector<std::string> lines;

    // temp variable
    std::string text;

    // Stream to use
    std::ifstream infile(filename.c_str());

    if (!infile)
        throw RuntimeError("BasisSetParser::parse: Unable to open basis set file: " + filename);

    while (getline(infile, text).good())
        lines.push_back(text);

    return lines;
}

std::vector<std::string> BasisSetParser::string_to_vector(const std::string &data)
{
    std::vector<std::string> ret;
    std::stringstream stream(data);
    std::string line;

    while(stream.good())
    {
        stream >> line;
        ret.push_back(line);
    }

    return ret;
}

std::vector<ShellInfo>
Gaussian94BasisSetParser::parse(const std::string& symbol, const std::vector<std::string> &lines)
{
    using namespace std;

    // s, p and s, p, d can be grouped together in Pople-style basis sets
    const string sp("SP"), spd("SPD");

    //                     a  b  c  d  e  f  g  h  i  j  k  l  m  n  o  p  q  r  s  t  u  v  w  x  y  z
    char shell_to_am[] = {-1,-1,-1, 2,-1, 3, 4, 5, 6,-1, 7, 8, 9,10,11, 1,12,13, 0,14,15,16,17,18,19,20};

    // Basis type.
    ShellInfo::GaussianType gaussian_type = ShellInfo::GaussianType::Pure;

    if (force_puream_or_cartesian_ && !forced_is_puream_)
        gaussian_type = ShellInfo::GaussianType::Cartesian;


    // Need a dummy center for the shell.
    Vector3 center;

    std::vector<ShellInfo> shell_list;

    int lineno = 0;
    bool found = false;
    bool newentry = true;

    while (lineno < static_cast<long int>(lines.size()))
    {
        string line = trim(lines[lineno++]);
        std::vector<std::string> splitline = string_to_vector(line);

        // Ignore blank lines
        if (line.empty())
            continue;

        // Look for ShellInfo::GaussianType::Cartesian or Spherical
        if (!force_puream_or_cartesian_)
        {
            if (to_lower_copy(line) == "cartesian")
            {
                gaussian_type = ShellInfo::GaussianType::Cartesian;
                continue;
            }
            else if (to_lower_copy(line) == "spherical")
            {
                gaussian_type = ShellInfo::GaussianType::Pure;
                continue;
            }
        } // end case where puream setting wasn't forced by caller

        // Do some matches
        if (line[0] == '!')
            continue;

        if (line == "****")
        {
            newentry = true;
            continue;
        }

        if (newentry)
        {
            newentry = false;

            // Check the captures and see if this basis set is for the atom we need.
            found = false;

            if (symbol == splitline[0])
            {
                found = true;

                // Read in the next line
                line = trim(lines[lineno++]);
                splitline = string_to_vector(line);

                // Need to do the following until we match a "****" which is the end of the basis set
                while (line != "****")
                {
                    // Match shell information
                    if (std::isalpha(splitline[0][0]))
                    {
                        std::string shell_type = splitline[0];
                        int nprimitive;
                        double scale;

                        std::transform(shell_type.begin(), shell_type.end(), shell_type.begin(), ::toupper);

                        double exponent, contraction;

                        if (!from_string<int>(nprimitive, splitline[1], std::dec))
                            throw RuntimeError("Gaussian94BasisSetParser::parse: Unable to convert number of primitives:\n" + line);
                        if (!from_string<double>(scale, splitline[2], std::dec))
                            throw RuntimeError("Gaussian94BasisSetParser::parse: Unable to convert scale factor:\n" + line);

                        if (shell_type.size() == 1)
                        {
                            int am = (int)shell_to_am[shell_type[0] - 'A'];

                            std::vector<double> exponents(nprimitive);
                            std::vector<double> contractions(nprimitive);

                            for (int p=0; p<nprimitive; ++p)
                            {
                                line = trim(lines[lineno++]);
                                splitline = string_to_vector(line);

                                int idx;
                                while((idx=line.find_first_of('D')) >= 0 )
                                {
                                    line.replace( idx, 1, "e" );
                                }
                                while((idx=line.find_first_of('d')) >= 0 )
                                {
                                    line.replace( idx, 1, "e" );
                                }

                                if (!from_string<double>(exponent, splitline[0], std::dec))
                                    throw RuntimeError("Gaussian94BasisSetParser::parse: Unable to convert exponent:\n" + line);
                                if (!from_string<double>(contraction, splitline[1], std::dec))
                                    throw RuntimeError("Gaussian94BasisSetParser::parse: Unable to convert contraction:\n" + line);

                                // Scale the contraction
                                contraction *= scale;

                                // Save the information.
                                exponents[p] = exponent;
                                contractions[p] = contraction;
                            }

                            //                                printf("Adding new shell. nprimitive = %d\n", nprimitive);
                            // We have a full shell, push it to the basis set
                            shell_list.push_back(ShellInfo(am, contractions, exponents, gaussian_type, 0, center));
                        }
                        else if (shell_type.size() == 2)
                        {
                            // This is to handle instances of SP, PD, DF, FG, ...
                            int am1 = (int)shell_to_am[shell_type[0] - 'A'];
                            int am2 = (int)shell_to_am[shell_type[1] - 'A'];

                            std::vector<double> exponents(nprimitive);
                            std::vector<double> contractions1(nprimitive);
                            std::vector<double> contractions2(nprimitive);

                            for (int p=0; p<nprimitive; ++p)
                            {
                                line = trim(lines[lineno++]);
                                splitline = string_to_vector(line);

                                int idx;
                                while((idx=line.find_first_of('D')) >= 0 )
                                {
                                    line.replace( idx, 1, "e" );
                                }
                                while((idx=line.find_first_of('d')) >= 0 )
                                {
                                    line.replace( idx, 1, "e" );
                                }

                                if (!from_string<double>(exponent, splitline[0], std::dec))
                                    throw RuntimeError("Gaussian94BasisSetParser::parse: Unable to convert exponent:\n" + line);
                                if (!from_string<double>(contraction, splitline[1], std::dec))
                                    throw RuntimeError("Gaussian94BasisSetParser::parse: Unable to convert first contraction:\n" + line);

                                // Scale the contraction
                                contraction *= scale;

                                // Save the information
                                exponents[p] = exponent;
                                contractions1[p] = contraction;

                                // Do the other contraction
                                if (!from_string<double>(contraction, splitline[2], std::dec))
                                    throw RuntimeError("Gaussian94BasisSetParser::parse: Unable to convert second contraction:\n" + line);

                                // Scale the contraction
                                contraction *= scale;

                                // Save the information
                                contractions2[p] = contraction;
                            }

                            //                                printf("Adding 2 new shells. nprimitive = %d\n", nprimitive);
                            shell_list.push_back(ShellInfo(am1, contractions1, exponents, gaussian_type, 0, center));
                            shell_list.push_back(ShellInfo(am2, contractions2, exponents, gaussian_type, 0, center));
                        }
                        else
                        {
                            throw RuntimeError("Gaussian94BasisSetParser::parse: Unable to parse basis sets with spd, or higher grouping\n");
                        }
                    }
                    else
                    {
                        throw RuntimeError("Gaussian94BasisSetParser::parse: Expected shell information, but got:\n" + line);
                    }
                    line = trim(lines[lineno++]);
                    splitline = string_to_vector(line);
                }
                break;
            }
        }
    }
    if (found == false)
        throw RuntimeError("Gaussian94BasisSetParser::parser: Unable to find the basis set for " + symbol + " in " + filename_);

    // The constructor, or the caller, should refresh the basis set.
    return shell_list;
}

} // end namespace panache


