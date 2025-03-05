#pragma once

#include <string>
#include <vector>
#include <iomanip>
#include <iostream>
#include <algorithm>

namespace donut {

    template<class T, class F, T... inds>
    constexpr void loop(std::integer_sequence<T, inds...>, F &&f) {
        (f(std::integral_constant<T, inds>{}), ...);
    }

    template<class T, T count, class F>
    constexpr void Loop(F &&f) {
        loop(std::make_integer_sequence<T, count>{}, std::forward<F>(f));
    }

    template<typename T>
    T max(const T &value) {
        return value;
    }

    template<typename T, typename...Args>
    T max(const T &value, Args &&...args) {
        return std::max(value, max(std::forward<Args>(args)...));
    }

    template<typename T>
    T min(const T &value) {
        return value;
    }

    template<typename T, typename...Args>
    T min(const T &value, Args &&...args) {
        return std::min(value, min(std::forward<Args>(args)...));
    }

    template<typename T, typename...Args>
    bool isSandwitchVal(const T &value, Args &&...args) {
        T minArg = min(std::forward<Args>(args)...);
        T maxArg = max(std::forward<Args>(args)...);
        //return (minArg < value && value < maxArg);
        return (minArg <= value && value <= maxArg);
    }

    template<typename T>
    void removeDupicates(std::vector<T> &vec) {
        auto endIter{vec.end()};
        for (auto iter = vec.begin(); iter != endIter; ++iter)
            endIter = std::remove(iter + 1, endIter, *iter);
        vec.erase(endIter, vec.end());
    }

} // namespace donut