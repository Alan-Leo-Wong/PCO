#pragma once

#include "Config.hpp"
#include <cstdio>
#include <cstdlib>
#include <string>
#include <sstream>
#include <stdexcept>

NAMESPACE_BEGIN(PCO)

#define PCO_FMT_CSTR(description, ...) description
#define PCO_FMT_STR(description, ...) std::string(description)
#define PCO_FMT_PRINT(description, ...) /*std::printf("%s\n", description)*/
#define PCO_FMT_ARG(arg)

#define PCO_ENABLE_ENSURE_HANDLER

#if defined(PCO_DISABLE_ENSURES)

#define PCO_ENSURE(expr, ...) ((void)0)

#elif defined(PCO_ENABLE_ENSURE_HANDLER)

    void ensureFailed(char const *function, char const *file, int line,
                      char const *description);

#define PCO_ENSURE(expr, description, ...)                                  \
  ((expr)                                                                      \
       ? ((void)0)                                                             \
       : ::PCO::ensureFailed(PCO_FUNCTION, __FILE__, __LINE__,           \
                                PCO_FMT_CSTR(description, ##__VA_ARGS__)))

#else

#define PCO_DEDAULT_ENSURE_FAILURE_IMPL(function, file, line, description, \
                                           ...)                                \
  do {                                                                         \
    std::printf("Assert failed in function '%s', "                             \
                "file '%s', line %d.\n",                                       \
                function, file, line);                                         \
    PCO_FMT_PRINT(description, ##__VA_ARGS__);                              \
    std::abort();                                                              \
  } while (0)

#ifdef __CUDACC__
#define PCO_ENSURE(expr, description, ...)                                  \
  do {                                                                         \
    if (!(expr)) {                                                             \
      std::printf("Assert failed in function '%s', file '%s', line %d.\n",     \
                  PCO_FUNCTION, __FILE__, __LINE__);                        \
      std::printf("%s", description);                                          \
      /* there is no std::abort in cuda kernels, hence we just print the error \
       * message here*/                                                        \
    }                                                                          \
  } while (0)
#else
#define PCO_ENSURE(expr, ...)                                               \
  do {                                                                         \
    if (!(expr)) {                                                             \
      PCO_DEDAULT_ENSURE_FAILURE_IMPL(PCO_FUNCTION, __FILE__, __LINE__,  \
                                         ##__VA_ARGS__);                       \
    }                                                                          \
  } while (0)
#endif

#endif

    template<typename T>
    void throw_LA_ASSERT(const T msg, const char *file, int line) {
        std::ostringstream oss;
        oss << "Robustness error: " << file << " line: " << line << " " << msg;
        throw std::runtime_error(oss.str());
    }

#define LA_GET_3RD_ARG_HELPER(arg1, arg2, arg3, ...) arg3

#define LA_Assert_1(cond)                                  \
    if (!(cond)) {                                         \
        simplicial_arrangement::throw_LA_ASSERT("", __FILE__, __LINE__); \
    }

#define LA_Assert_2(cond, msg)                              \
    if (!(cond)) {                                          \
        simplicial_arrangement::throw_LA_ASSERT(msg, __FILE__, __LINE__); \
    }

#if defined(_MSC_VER) // special VS variadics

#define LA_MSVS_EXPAND(x) x

#define LA_ASSERT(...)                                                           \
    LA_MSVS_EXPAND(LA_GET_3RD_ARG_HELPER(__VA_ARGS__, LA_Assert_2, LA_Assert_1)) \
    LA_MSVS_EXPAND((__VA_ARGS__))

#else
#define LA_ASSERT(...) LA_GET_3RD_ARG_HELPER(__VA_ARGS__, LA_Assert_2, LA_Assert_1)(__VA_ARGS__)
#endif

#ifdef PCO_NON_ROBUST
#define PCO_ASSERT(x) do { LA_ASSERT(x); } while(0)
#else
#define PCO_ASSERT(x) do { assert(x); } while(0)
#endif

NAMESPACE_END(PCO)