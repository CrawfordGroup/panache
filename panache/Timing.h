#ifndef PANACHE_TIMING_H
#define PANACHE_TIMING_H

#include <chrono>


namespace panache
{
namespace timing
{

using panacheclock = std::chrono::high_resolution_clock;

class Timer
{
private:
    panacheclock::time_point begin, end;
    panacheclock::duration dur;
    unsigned long int nmeasurements;

public:
    Timer()
    {
        Reset();
    }

    void Start(void)
    {
        begin = panacheclock::now();
    }

    void Stop(void)
    {
        end = panacheclock::now();
        dur += (end - begin);
        nmeasurements++;
    }

    void Lap(void)
    {
        dur += panacheclock::now() - begin;
        nmeasurements++;
    }

    void Reset(void)
    {
        dur = panacheclock::duration::zero();
        nmeasurements = 0;
    }

    unsigned long int Seconds(void) const
    {
        return std::chrono::duration_cast<std::chrono::seconds>(dur).count();
    }

    unsigned long int Milliseconds(void) const
    {
        return std::chrono::duration_cast<std::chrono::milliseconds>(dur).count();
    }

    unsigned long int Microseconds(void) const
    {
        return std::chrono::duration_cast<std::chrono::microseconds>(dur).count();
    }

    unsigned long int TimesCalled(void) const
    {
        return nmeasurements;
    }

};


}} // close namespace panache::timing

#endif // PANACHE_TIMING_H
