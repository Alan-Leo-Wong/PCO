#pragma once

#if defined(__clang__) || defined(__GNUC__)
#  define PCO_INLINE __attribute__((always_inline)) inline
#elif defined(_MSC_VER)
#  define PCO_INLINE __forceinline
#endif

/* namespace macro */
#define PCO offset

#if !defined(NAMESPACE_BEGIN)
#  define NAMESPACE_BEGIN(name) namespace name {
#endif

#if !defined(NAMESPACE_END)
#  define NAMESPACE_END(name) }
#endif

#ifdef __GNUC__
#define PCO_FUNCTION __PRETTY_FUNCTION__
#elif defined(__clang__) || (_MSC_VER >= 1310)
#define PCO_FUNCTION __FUNCTION__
#else
#define PCO_FUNCTION "unknown"
#endif