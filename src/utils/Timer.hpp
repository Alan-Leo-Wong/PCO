//
// Created by lei on 23-9-15.
//
#pragma once

#include "Config.hpp"

#include <iostream>

#if defined(WIN32) || defined(_WIN32)
#define WINDOWS_LEAN_AND_MEAN

#include <windows.h>

#undef min
#undef max
#elif defined(__linux__)

#include <sys/time.h>
#include <ctime>

#endif // WIN32

NAMESPACE_BEGIN(PCO)
    namespace utils::timer {
        class TimerInterface {
        public:
            TimerInterface() = default;

            virtual ~TimerInterface() {}

        public:
            // Start time measurement
            virtual void start() = 0;

            // Stop time measurement
            virtual void stop() = 0;

            // Reset time counters to zero
            virtual void reset() = 0;

            // Get diff time when stop timer
            virtual double getDiffTime() = 0;

            // Time in msec. after start. If the stop watch is still running (i.e. there
            // was no call to stop()) then the elapsed time is returned, otherwise the
            // time between the last start() and stop call is returned
            virtual double getAllTime() = 0;

            // Mean time to date based on the number of times the stopwatch has been
            // _stopped_ (ie finished sessions) and the current total time
            virtual double getAverageTime() = 0;
        };

#if defined(WIN32) || defined(_WIN32)

        class WinTimer : public TimerInterface {
        public:
            WinTimer() : start_time(), end_time(),
                         diff_time(0.0), total_time(0.0),
                         isRunning(false), sessions(0),
                         isSetFreq(false) {
                if (!isSetFreq) {
                    isSetFreq = true;
                    // get the tick frequency from the OS
                    QueryPerformanceFrequency(&freq);
                }
            }

            ~WinTimer() {}

        public:
            inline void start();

            inline void stop();

            inline void reset();

            inline double getDiffTime();

            inline double getAllTime();

            inline double getAverageTime();

        private:
            LARGE_INTEGER start_time;

            LARGE_INTEGER end_time;

            // tick frequency
            LARGE_INTEGER freq;

            double diff_time;

            double total_time;

            // flag if the stop watch is running
            bool isRunning;

            // Number of times clock has been started and stopped to allow averaging
            int sessions;


            // flag if the frequency has been set
            bool isSetFreq;
        };

        inline void WinTimer::start() {
            QueryPerformanceCounter(reinterpret_cast<LARGE_INTEGER *>(&start_time));
            isRunning = true;
        }

        inline void WinTimer::stop() {
            QueryPerformanceCounter(reinterpret_cast<LARGE_INTEGER *>(&end_time));
            diff_time = (double) (((double) (end_time.QuadPart) -
                                   (double) (start_time.QuadPart)) * 1000.0 /
                                  freq.QuadPart);
            total_time += diff_time;
            sessions++;
            isRunning = false;
        }

        inline void WinTimer::reset() {
            diff_time = 0.;
            total_time = 0.;
            sessions = 0;
            /* Does not change the timer running state but does recapture this point
               in time as the current start time if it is running. */
            if (isRunning)
                QueryPerformanceCounter(reinterpret_cast<LARGE_INTEGER *>(&start_time));
        }

        ////////////////////////////////////////////////////////////////////////////////
        //! Time in msec. after start. If the stop watch is still running (i.e. there
        //! was no call to stop()) then returns 0.0,
        //! otherwise the current elapsed time difference alone is returned.
        ////////////////////////////////////////////////////////////////////////////////
        inline double WinTimer::getDiffTime() {
            if (!isRunning) return diff_time;
            else {
                printf("Warning: You must STOP timer first!\n");
                return .0;
            }
        }

        ////////////////////////////////////////////////////////////////////////////////
        //! Time in msec. after start. If the stop watch is still running (i.e. there
        //! was no call to stop()) then the elapsed time is returned added to the
        //! current diff_time sum, otherwise the current summed time difference alone
        //! is returned.
        ////////////////////////////////////////////////////////////////////////////////
        inline double WinTimer::getAllTime() {
            double retval = total_time;
            if (isRunning) {
                LARGE_INTEGER temp;
                QueryPerformanceCounter(reinterpret_cast<LARGE_INTEGER *>(&temp));
                retval += (double) (((double) (temp.QuadPart) -
                                     (double) (start_time.QuadPart)) * 1000 /
                                    freq.QuadPart);
            }
            return retval;
        }

        inline double WinTimer::getAverageTime() {
            /* flag 'isRunning' is false when call the function */
            return (sessions > 0) ? (total_time / sessions) : 0.0;
        }

#elif defined(__linux__)

        class LinuxTimer : public TimerInterface {
        public:
            LinuxTimer() : start_time(),
                           diff_time(0.0), total_time(0.0),
                           isRunning(false), sessions(0) {}

            ~LinuxTimer() {}

        public:
            inline void start();

            inline void stop();

            inline void reset();

            inline double getAllTime();

            inline double getAverageTime();

            inline double getDiffTime();

        private:
            struct timeval start_time;

            double diff_time;

            double total_time;

            // flag if the stop watch is running
            bool isRunning;

            // Number of times clock has been started and stopped to allow averaging
            int sessions;
        };

        inline void LinuxTimer::start() {
            gettimeofday(&start_time, 0);
            isRunning = true;
        }

        inline void LinuxTimer::stop() {
            diff_time = getDiffTime();
            total_time += diff_time;
            sessions++;
            isRunning = false;
        }

        inline void LinuxTimer::reset() {
            diff_time = 0.;
            total_time = 0.;
            sessions = 0;
            /* Does not change the timer running state but does recapture this point
               in time as the current start time if it is running. */
            if (isRunning)
                gettimeofday(&start_time, 0);
        }

        inline double LinuxTimer::getAllTime() {
            double retval = total_time;
            if (isRunning)
                retval += getDiffTime();
            return retval;
        }

        inline double LinuxTimer::getAverageTime() {
            /* flag 'isRunning' is false when call the function */
            return (sessions > 0) ? (total_time / sessions) : 0.0;
        }

        inline double LinuxTimer::getDiffTime() {
            struct timeval t_time;
            gettimeofday(&t_time, 0);

            // time difference in milli-seconds
            return (double) (1e3 * (t_time.tv_sec - start_time.tv_sec) +
                             (1e-3 * (t_time.tv_usec - start_time.tv_usec)));
        }

#endif // WIN32

        inline bool createTimer(TimerInterface **timer_interface) {
#if defined(WIN32) || defined(_WIN32) || defined(WIN64) || defined(_WIN64)
            *timer_interface = reinterpret_cast<TimerInterface *>(new WinTimer());
#else
            *timer_interface = reinterpret_cast<TimerInterface *>(new LinuxTimer());
#endif
            return (*timer_interface) ? true : false;
        }

        inline bool deleteTimer(TimerInterface **timer_interface) {
            if (*timer_interface) {
                delete *timer_interface;
                *timer_interface = nullptr;
                return true;
            } else {
                printf("Can not DELETE a null timer!\n");
                return false;
            }
        }

        inline bool startTimer(TimerInterface **timer_interface) {
            if (*timer_interface) {
                (*timer_interface)->start();
                return true;
            } else {
                printf("Can not START a null timer!\n");
                return false;
            }
        }

        inline bool stopTimer(TimerInterface **timer_interface) {
            if (*timer_interface) {
                (*timer_interface)->stop();
                return true;
            } else {
                printf("Can not STOP a null timer!\n");
                return false;
            }
        }

        inline bool resetTimer(TimerInterface **timer_interface) {
            if (*timer_interface) {
                (*timer_interface)->reset();
                return true;
            } else {
                printf("Can not RESET a null timer!\n");
                return false;
            }
        }

        inline double getElapsedTime(TimerInterface **timer_interface) {
            if (*timer_interface) {
                return (*timer_interface)->getDiffTime();
            } else {
                printf("Can not GET TIME from a null timer!\n");
                return 0.0;
            }
        }

        ////////////////////////////////////////////////////////////////////////////////
        //! Total execution time for the timer over all runs since the last reset
        //! or timer creation.
        //! @param name  name of the timer to obtain the value of.
        ////////////////////////////////////////////////////////////////////////////////
        inline double getAllTimeValue(TimerInterface **timer_interface) {
            if (*timer_interface) {
                return (*timer_interface)->getAllTime();
            } else {
                printf("Can not GET TIME from a null timer!\n");
                return 0.0;
            }
        }

        inline double getAverageTimerValue(TimerInterface **timer_interface) {
            if (*timer_interface) {
                return (*timer_interface)->getAverageTime();
            } else {
                printf("Can not GET AVERAGE TIME from a null timer!\n");
                return 0.0;
            }
        }

    } // namespace utils::timer
NAMESPACE_END(PCO)