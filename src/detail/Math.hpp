#pragma once

#include "Config.hpp"
#include "Handler.hpp"
#include "BasicDataType.hpp"
#include <type_traits>
#include <iostream>

NAMESPACE_BEGIN(offset)
    // Use standard mathematical constants' M_PI if available
#ifdef M_PI
    constexpr double PI = M_PI;
#else
    constexpr double PI = 3.1415926535897932384626433832795;
#endif

    template<class Scalar_>
    class MaxMetric {
    public:
        static Scalar_ impl(Scalar_ s0, Scalar_ s1) {
            using std::abs;
            return abs(s0 - s1);
        }
    };

    template<class Scalar_, int M, int N>
    class MaxMetric<Matrix<Scalar_, M, N>> {
    public:
        static Scalar_ impl(Matrix<Scalar_, M, N> const &p0,
                            Matrix<Scalar_, M, N> const &p1) {
            return (p0 - p1).template lpNorm<Eigen::Infinity>();
        }
    };

    template<class Scalar_>
    class SetToZero {
    public:
        static void impl(Scalar_ &s) { s = Scalar_(0); }
    };

    template<class Scalar_, int M, int N>
    class SetToZero<Matrix<Scalar_, M, N>> {
    public:
        static void impl(Matrix<Scalar_, M, N> &v) { v.setZero(); }
    };

    template<class T1, class Scalar_>
    class SetElementAt;

    template<class Scalar_>
    class SetElementAt<Scalar_, Scalar_> {
    public:
        static void impl(Scalar_ &s, Scalar_ value, int at) {
            PCO_ENSURE(at == 0, "is {}", PCO_FMT_ARG(at));
            s = value;
        }
    };

    template<class Scalar_, int N>
    class SetElementAt<Vector<Scalar_, N>, Scalar_> {
    public:
        static void impl(Vector<Scalar_, N> &v, Scalar_ value, int at) {
            //OFFSET_ENSURE(at >= 0 && at < N, "is {}", OFFSET_FMT_ARG(at));
            v[at] = value;
        }
    };

    template<class Scalar_>
    class SquaredNorm {
    public:
        static Scalar_ impl(Scalar_ const &s) { return s * s; }
    };

    template<class Scalar_, int N>
    class SquaredNorm<Matrix<Scalar_, N, 1>> {
    public:
        static Scalar_ impl(Matrix<Scalar_, N, 1> const &s) { return s.squaredNorm(); }
    };

    template<class Scalar_>
    class Transpose {
    public:
        static Scalar_ impl(Scalar_ const &s) { return s; }
    };

    template<class Scalar_, int M, int N>
    class Transpose<Matrix<Scalar_, M, N>> {
    public:
        static Matrix<Scalar_, M, N> impl(Matrix<Scalar_, M, N> const &s) {
            return s.transpose();
        }
    };

    /// Returns maximum metric between two points ``p0`` and ``p1``, with ``p0, p1``
    /// being matrices or a _Scalars.
    ///
    template<class T>
    auto maxMetric(T const &p0, T const &p1)
    -> decltype(MaxMetric<T>::impl(p0, p1)) {
        return MaxMetric<T>::impl(p0, p1);
    }

    /// Sets point ``p`` to zero, with ``p`` being a matrix or a Scalar_.
    ///
    template<class T>
    void setToZero(T &p) {
        return SetToZero<T>::impl(p);
    }

    /// Sets ``i``th component of ``p`` to ``value``, with ``p`` being a
    /// matrix or a Scalar_. If ``p`` is a Scalar_, ``i`` must be ``0``.
    ///
    template<class T, class Scalar_>
    void setElementAt(T &p, Scalar_ value, int i) {
        return SetElementAt<T, Scalar_>::impl(p, value, i);
    }

    /// Returns the squared 2-norm of ``p``, with ``p`` being a vector or a Scalar_.
    ///
    template<class T>
    auto squaredNorm(T const &p) -> decltype(SquaredNorm<T>::impl(p)) {
        return SquaredNorm<T>::impl(p);
    }

    /// Returns ``p.transpose()`` if ``p`` is a matrix, and simply ``p`` if m is a
    /// Scalar_.
    ///
    template<class T>
    auto transpose(T const &p) -> decltype(Transpose<T>::impl(T())) {
        return Transpose<T>::impl(p);
    }

    template<class Scalar_>
    struct IsFloatingPoint {
        static bool const value = std::is_floating_point<Scalar_>::value;
    };

    template<class Scalar_, int M, int N>
    struct IsFloatingPoint<Matrix<Scalar_, M, N>> {
        static bool const value = std::is_floating_point<Scalar_>::value;
    };

    template<class _Scalar_>
    struct Get_Scalar {
        using Scalar_ = _Scalar_;
    };

    template<class _Scalar_, int M, int N>
    struct Get_Scalar<Matrix<_Scalar_, M, N>> {
        using Scalar_ = _Scalar_;
    };

    /// If the Vector type is of fixed size, then IsFixedSizeVector::value will be
    /// true.
    template<typename Vector, int NumDimensions,
            typename = typename std::enable_if<
                    Vector::RowsAtCompileTime == NumDimensions &&
                    Vector::ColsAtCompileTime == 1>::type>
    struct IsFixedSizeVector : std::true_type {
    };

    namespace {
        // 'GEO_EPSILON' is a static definition in anonymous namespace; static is redundant here
        /*static */constexpr Scalar GEO_EPSILON = 1E-6;
    }

    template<typename T>
    int dcmp(const T &x, const T &epsilon) {
        if (std::abs(x) < epsilon) return 0;
        else return x < 0 ? -1 : 1;
    }

    PCO_INLINE double safe_acos(double value) {
        if (value <= -1.0) {
            return PI;
        } else if (value >= 1.0) {
            return 0;
        } else {
            return std::acos(value);
        }
    }

    template<typename T>
    void traverseSparseMat(const _SpMatrix<T> &spMat) {
        // _SpMatrix是列优先存储，所以优先遍历每列的非空元素，即使第一层for循环使用了rows()
        for (int k = 0; k < /*spMat.rows()*/spMat.outerSize(); ++k) {
            for (typename _SpMatrix<T>::InnerIterator it(spMat, k); it; ++it) {
                // it.value();
                // it.row();   // row index
                // it.col();   // col index(here it is equal to k)
                // it.index(); // inner index(here it is equal to it.row())
                std::cout << "row = " << it.row() << ", ";
                std::cout << "col = " << it.col() << ", ";
                std::cout << "value = " << it.value() << std::endl;
            }
        }
    }

    template<typename T>
    struct IndexedValue {
        int index;
        T value;

        IndexedValue() = default;

        IndexedValue(int idx, const T &val) : index(idx), value(val) {}

        // 按照值进行升序排序(大顶堆用这个, 从小到大排序(包括set)用这个)
        /*
         * 如 std::priority_queue<IndexedValue<int>, std::vector<IndexedValue<int>>,
         *    decltype(&IndexedValue<int>::asc_cmp)> maxHeap(&IndexedValue<int>::asc_cmp);
        */
        static bool asc_cmp(const IndexedValue<T> &a, const IndexedValue<T> &b) {
            return /*a.index == b.index ? */(a.value < b.value);
        }

        // 按照值进行降序排序(小顶堆用这个, 从大到小排序(包括set)用这个)
        static bool dsc_cmp(const IndexedValue<T> &a, const IndexedValue<T> &b) {
            return a.value > b.value;
        }

        bool operator<(const IndexedValue<T> &other) const {
            return (index == other.index) ? (value < other.value) : (index < other.index);
            //                return asc_cmp(*this, other);
        }

        template<typename T1>
        friend std::ostream &operator<<(std::ostream &os, const IndexedValue<T1> &iv);
    };

    template<typename T1>
    std::ostream &operator<<(std::ostream &os, const IndexedValue<T1> &iv) {
        os << "index = " << iv.index << ", value = " << iv.value;
        return os;
    }

    PCO_INLINE Vector3
    getBarycentric(const Vector3 &p, const Vector3 &tri_vert_0, const Vector3 &tri_vert_1, const Vector3 &tri_vert_2) {
        Vector3 baryCoord;
        // compute barycentric coordinates(fast but a bit inaccurate)
        // http://gamedev.stackexchange.com/a/23745
        {
            Vector3 v0 = tri_vert_1 - tri_vert_0;
            Vector3 v1 = tri_vert_2 - tri_vert_0;
            Vector3 v2 = p - tri_vert_0;
            Scalar d00 = v0.dot(v0);
            Scalar d01 = v0.dot(v1);
            Scalar d11 = v1.dot(v1);
            Scalar d20 = v2.dot(v0);
            Scalar d21 = v2.dot(v1);
            Scalar denom = d00 * d11 - d01 * d01;

            baryCoord(1) = (d11 * d20 - d01 * d21) / denom;
            baryCoord(2) = (d00 * d21 - d01 * d20) / denom;
            baryCoord(0) = 1.0 - baryCoord(1) - baryCoord(2);
        }
        return baryCoord;
    }

    PCO_INLINE Vector3 getBarycentric(const Vector3 &p, const MatrixX &tri_verts) {
        return getBarycentric(p, tri_verts.row(0), tri_verts.row(1), tri_verts.row(2));
    }

NAMESPACE_END(offset)