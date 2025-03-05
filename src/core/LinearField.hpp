//
// Created by Lei on 3/30/2024.
//

#ifndef PCO_LINEARFIELD_HPP
#define PCO_LINEARFIELD_HPP

#include "Config.hpp"
#include "detail/Geometry.hpp"
#include "detail/BasicDataType.hpp"

NAMESPACE_BEGIN(PCO)
namespace core {
	using TetDis = Vector4; // We use the distance of a tet's four vertexes to represent a linear field
	using LFGradient = Vector3;
	using namespace geometry;

	struct LinearField {
	public:
		LinearField() = default;

		LinearField(const TetDis& _field, const Tet& tet) : field(_field) {
			initialize(tet);
		}

	private:
		void initialize(const Tet& tet) {
			MatrixX A(4, 4);

			A.row(0) = Vector4(tet(0, 0), tet(0, 1), tet(0, 2), 1);
			A.row(1) = Vector4(tet(1, 0), tet(1, 1), tet(1, 2), 1);
			A.row(2) = Vector4(tet(2, 0), tet(2, 1), tet(2, 2), 1);
			A.row(3) = Vector4(tet(3, 0), tet(3, 1), tet(3, 2), 1);

			gradient = (A.inverse() * field).head<3>().normalized();
		}

	public:
		TetDis field;
		LFGradient gradient;
	};

	bool isIntersect(const LinearField& lf_1, const LinearField& lf_2, double offsetDis);

	bool isIntersect(const LinearField& lf_1, const LinearField& lf_2, const Tet& tetQ, double offsetDis, int tetType = 0);

	bool isHaveValidUnion(const LinearField& lf_1, const LinearField& lf_2, double offsetDis);

	bool isHaveValidUnion(const LinearField& lf_1, const LinearField& lf_2, const Tet& tetQ);

	Scalar gradientAngle(const LinearField& lf_1, const LinearField& lf_2);

} // namespace core
NAMESPACE_END(PCO)

#endif //PCO_LINEARFIELD_HPP
