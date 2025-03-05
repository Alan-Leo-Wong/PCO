//
// Created by Lei on 3/31/2024.
//

#ifndef PCO_SIMPLEXCUT_HPP
#define PCO_SIMPLEXCUT_HPP

#include "Config.hpp"
#include "Complex.hpp"
#include "PlaneRepo.hpp"

NAMESPACE_BEGIN(PCO)
namespace complex {

	class ComplexCut {
	public:
		ComplexCut() noexcept = default;

        ComplexCut(const Complex& _complex) : complex(_complex) {}

		ComplexCut(const PlaneRepo& _planes, const Complex& _complex)
			: planes(_planes), complex(_complex) {}

	private:
		int8_t cut_0_face(size_t vid,
			size_t plane_index);

		std::array<size_t, 3> cut_1_face(size_t eid,
			size_t plane_index,
			const std::vector<int8_t>& orientations);

		std::array<size_t, 3> cut_2_face(size_t fid,
			size_t plane_index,
			const std::vector<int8_t>& orientations,
			const std::vector<std::array<size_t, 3>>& subedges);

		std::array<size_t, 3> cut_3_face(size_t cid,
			size_t plane_index,
			const std::vector<std::array<size_t, 3>>& subfaces);

		void consolidate();

	public:
		/**
		 * Insert a plane into the existing complex.
		 *
		 * @param[in] plane_index     The index of the plane to be inserted.
		 * @param[in] is_positive     The index of the plane to be inserted.
		 *
		 * @return The index of an existing plane that is coplanar with the inserted
		 *         plane if exists. Otherwise, return `INVALID`.
		 */
		size_t add_plane(size_t plane_index, bool is_positive);

		size_t add_plane(const Plane& plane, bool is_positive);

		const Complex& get_complex() const { return complex; }

		Complex& get_complex() { return complex; }

		Complex&& export_complex()&& { return std::move(complex); }

	private:
		Complex complex;

		PlaneRepo planes;
	};

} // namespace complex
NAMESPACE_END(PCO)

#endif //PCO_SIMPLEXCUT_HPP
