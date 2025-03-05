#pragma once

#include "Config.hpp"
#include <numeric>
#include <array>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <unordered_set>

NAMESPACE_BEGIN(PCO)
typedef unsigned int uint;

#ifdef SINGLE
typedef float Scalar;
const Scalar DINF = std::numeric_limits<float>::max();
const Scalar DMIN = std::numeric_limits<float>::min(); // the minimum positive value
const Scalar DNMIN = std::numeric_limits<float>::lowest();
#else
typedef double Scalar;
const Scalar DINF = std::numeric_limits<double>::max();
const Scalar DMIN = std::numeric_limits<double>::min(); // the minimum positive value
const Scalar DNMIN = std::numeric_limits<double>::lowest();
#endif

template<typename REAL, int Dim>
using Vector = typename Eigen::Matrix<REAL, Dim, 1>;

using VectorX = Eigen::Matrix<Scalar, Eigen::Dynamic, 1>;

using VectorXi = Eigen::VectorXi;

using Vector2 = Eigen::Matrix<Scalar, 2, 1>;

using Vector3 = Eigen::Matrix<Scalar, 3, 1>;

using Vector3f = Eigen::Matrix<float, 3, 1>;

using Vector4 = Eigen::Matrix<Scalar, 4, 1>;

using Vector2i = Eigen::Vector2i;

using Vector3i = Eigen::Vector3i;

using Array3 = Eigen::Array<Scalar, 3, 1>;

using Array4 = Eigen::Array<Scalar, 4, 1>;

using Array4i = Eigen::Array<Scalar, 4, 1>;

template<typename REAL, int Rows, int Cols>
using Matrix = typename Eigen::Matrix<REAL, Rows, Cols>;

using MatrixX = Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>;

using Matrix3 = Eigen::Matrix<Scalar, 3, 3>;

using MatrixXi = Eigen::MatrixXi;

template<typename REAL, int Num>
using ArrayX = typename Eigen::Array<REAL, Num, 1>;

template<typename REAL, int Rows, int Cols>
using ArrayXX = typename Eigen::Array<REAL, Rows, Cols>;

template<typename REAL>
using _SpMatrix = Eigen::SparseMatrix<REAL>;
using SpMatrix = _SpMatrix<Scalar>;
using SpMatrixi = _SpMatrix<int>;

template<typename REAL>
using _SpTriplet = Eigen::Triplet<REAL>;
using SpTriplet = _SpTriplet<Scalar>;
using SpITriplet = _SpTriplet<int>;

NAMESPACE_END(PCO)

namespace std {
	template<typename REAL, size_t N>
	struct hash<array<REAL, N> > {
		typedef array<REAL, N> argument_type;
		typedef size_t result_type;

		result_type operator()(const argument_type& a) const {
			hash<REAL> hasher;
			result_type h = 0;
			for (result_type i = 0; i < N; ++i) {
				h = h * 31 + hasher(a[i]);
			}
			return h;
		}
	};
}

template<>
struct std::less<PCO::ArrayX<PCO::Scalar, 4>> {
public:
	using Array4 = PCO::ArrayX<PCO::Scalar, 4>;

	bool operator()(const Array4& a, const Array4& b) const {
		for (int i = 0; i < a.size(); ++i) {
			if (fabs(a[i] - b[i]) < PCO::DMIN) continue;

			if (a[i] < b[i]) return true;
			else if (a[i] > b[i]) return false;
		}
		return false;
	}
};

template<>
struct std::less<PCO::VectorX> {
public:
	using VectorX = PCO::VectorX;

	bool operator()(const VectorX& a, const VectorX& b) const {
		for (int i = 0; i < a.size(); ++i) {
			if (fabs(a[i] - b[i]) < PCO::DMIN) continue;

			if (a[i] < b[i]) return true;
			else if (a[i] > b[i]) return false;
		}
		return false;
	}
};

template<>
struct std::less<PCO::Vector3> {
public:
	using Vector3 = PCO::Vector3;

	bool operator()(const Vector3& a, const Vector3& b) const {
		for (int i = 0; i < 3; ++i) {
			if (fabs(a[i] - b[i]) < PCO::DMIN) continue;

			if (a[i] < b[i]) return true;
			else if (a[i] > b[i]) return false;
		}
		return false;
	}
};

template<>
struct std::less<PCO::Vector3i> {
public:
    using Vector3i = PCO::Vector3i;

    bool operator()(const Vector3i& a, const Vector3i& b) const {
        for (int i = 0; i < 3; ++i) {
            if (a[i] == b[i]) continue;

            if (a[i] < b[i]) return true;
            else if (a[i] > b[i]) return false;
        }
        return false;
    }
};

template<>
struct std::less<PCO::Vector4> {
public:
	using Vector4 = PCO::Vector4;

	bool operator()(const Vector4& a, const Vector4& b) const {
		for (int i = 0; i < 4; ++i) {
			if (fabs(a[i] - b[i]) < PCO::DMIN) continue;

			if (a[i] < b[i]) return true;
			else if (a[i] > b[i]) return false;
		}
		return false;
	}
};
