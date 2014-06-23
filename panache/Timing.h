/*!
 * \file
 * \brief Simple timer class using C++11 clocks
 * \author Benjamin Pritchard (ben@bennyp.org)
 */


#ifndef PANACHE_TIMING_H
#define PANACHE_TIMING_H

#include <chrono>


namespace panache
{


/*!
 * \brief Functions for timing
 */
namespace timing
{

using panacheclock = std::chrono::high_resolution_clock;


/*!
 * \brief Simple timer class
 *
 * This class uses the C++11 high-resolution clock facility
 * to implement a simple timer
 */
class Timer
{
private:
    panacheclock::time_point begin;   //!< The last time the clock was started
    panacheclock::time_point end;     //!< The (last) time the clocked was stopped
    panacheclock::duration dur;       //!< Accumulated duration (all the time between starts and stops)
    unsigned long int nmeasurements;  //!< Number of times the timer has been stopped

public:
    /*!
     * \brief Creates a timer
     */
    Timer()
    {
        Reset();
    }


    /*!
     * \brief Start the timer
     *
     * Sets the start time to now
     */
    void Start(void)
    {
        begin = panacheclock::now();
    }


    /*!
     * \brief Stop the timer
     *
     * Adds the time between now and the last start to the
     * accumulated duration, and increments nmeasurements
     */
    void Stop(void)
    {
        end = panacheclock::now();
        dur += (end - begin);
        nmeasurements++;
    }


    /*!
     * \brief Adds the elapsed time to accumulated duration
     *
     * Also increases the number of measurements
     */
    void Lap(void)
    {
        dur += panacheclock::now() - begin;
        nmeasurements++;
    }

    /*!
     * \brief Zeroes the duration and number of measurements
     */
    void Reset(void)
    {
        dur = panacheclock::duration::zero();
        nmeasurements = 0;
    }


    /*!
     * \brief Returns the accumulated duration in seconds
     */
    unsigned long int Seconds(void) const
    {
        return std::chrono::duration_cast<std::chrono::seconds>(dur).count();
    }


    /*!
     * \brief Returns the accumulated duration in milliseconds
     */
    unsigned long int Milliseconds(void) const
    {
        return std::chrono::duration_cast<std::chrono::milliseconds>(dur).count();
    }


    /*!
     * \brief Returns the accumulated duration in microseconds
     */
    unsigned long int Microseconds(void) const
    {
        return std::chrono::duration_cast<std::chrono::microseconds>(dur).count();
    }


    /*!
     * \brief Returns the number of durations added
     */
    unsigned long int TimesCalled(void) const
    {
        return nmeasurements;
    }

};


}
} // close namespace panache::timing

#endif // PANACHE_TIMING_H

