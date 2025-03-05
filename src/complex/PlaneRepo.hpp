#pragma once

#include "detail/BasicDataType.hpp"

#include <vector>

NAMESPACE_BEGIN(PCO)
namespace complex {

	/**
	 * A plane is defined by the barycentric plane equation:
	 *     f0 * b0 + f1 * b1 + f2 * b2 + f3 * b3 = 0 // For 3D
	 *     f0 * b0 + f1 * b1 + f2 * b2 = 0           // For 2D
	 * where the b's are the barycentric variables, and f's are the
	 * plane equation coefficients.  We store the f's for each plane.
	 */
	using Plane = Eigen::Matrix<Scalar, DIM + 1, 1>;

	class PlaneRepo {
	public:
		PlaneRepo() noexcept {
			initialize_simplicial_planes();
		}

		explicit PlaneRepo(const std::vector<Plane>& planes)
			: m_planes(planes) {
			initialize_simplicial_planes();
		}

		void add_new_plane(const Plane& plane) {
			m_planes.emplace_back(plane);
		}

		const Plane& get_plane(size_t i) const {
			if (i <= DIM) {
				return m_simplex_planes[i];
			}
			else {
				return m_planes[i - DIM - 1];
			}
		}

		size_t get_num_planes() const { return m_planes.size(); }

	private:
		void initialize_simplicial_planes() {
			m_simplex_planes[0] = { 1, 0, 0, 0 };
			m_simplex_planes[1] = { 0, 1, 0, 0 };
			m_simplex_planes[2] = { 0, 0, 1, 0 };
			m_simplex_planes[3] = { 0, 0, 0, 1 };
		}

	private:
		std::array<Plane, DIM + 1> m_simplex_planes;
		std::vector<Plane> m_planes;
	};

} // namespace complex
NAMESPACE_END(PCO)
