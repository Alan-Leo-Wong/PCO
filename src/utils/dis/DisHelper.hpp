#pragma once

#include "detail/Geometry.hpp"
#include "detail/BasicDataType.hpp"

#include <vector>

#include <omp.h>
//#include <PQP/Build.h>
#include <igl/signed_distance.h>
#include <igl/point_mesh_squared_distance.h>

NAMESPACE_BEGIN(PCO)
namespace utils::dishelper {

	PCO_INLINE
		Scalar simplexSqrDis(const Vector3& p, const MatrixX& V) {
		MatrixXi F(1, 3);
		F << 0, 1, 2;
		MatrixXi::Index i = 0;
		Scalar sqrD;
		MatrixX C;
		igl::point_simplex_squared_distance<3>(p, V, F, i, sqrD, C);
		return sqrD;
	}

	PCO_INLINE
		VectorX simplexMultiSqrDis(const MatrixX& ps, const MatrixX& V) {
		int num_p = ps.rows();
		VectorX dis(num_p);
#pragma omp parallel for
		for (int i = 0; i < num_p; ++i)
			dis(i) = simplexSqrDis(ps.row(i), V);
		return dis;
	}

	PCO_INLINE
		VectorX simplexMultiSqrDis(const MatrixX& ps, const MatrixX& V, MatrixX& C) {
		MatrixXi F(1, 3);
		F << 0, 1, 2;
		MatrixXi::Index index = 0;

		int num_p = ps.rows();
		VectorX dis(num_p);
#pragma omp parallel for
		for (int i = 0; i < num_p; ++i) {
			MatrixX c;
			igl::point_simplex_squared_distance<3>(ps.row(i), V, F, index, dis(i), c);
			C.row(i) = c.row(0);
		}
		return dis;
	}

	PCO_INLINE
		VectorX simplexLpDis(const MatrixX& ps, const MatrixX& V, int lp = 2) {
		MatrixXi F(1, 3);
		F << 0, 1, 2;
		MatrixXi::Index index = 0;

		int num_p = ps.rows();
		VectorX dis(num_p);
#pragma omp parallel for
		for (int i = 0; i < num_p; ++i) {
			Vector3 p = ps.row(i);
			MatrixX c;
			igl::point_simplex_squared_distance<3>(ps.row(i), V, F, index, dis(i), c);
			//if (lp == 2) continue;

			double x0, x1, x2;
			Vector3 _c = c.row(0);
			switch (lp) {
			case 1:
				dis(i) = std::abs(p(0) - _c(0)) + std::abs(p(1) - _c(1)) + std::abs(p(2) - _c(2));
				break;
			case 2:
				dis(i) = std::sqrt(dis(i));
				break;
			case 3:
				x0 = std::pow(p(0) - _c(0), 3);
				x1 = std::pow(p(1) - _c(1), 3);
				x2 = std::pow(p(2) - _c(2), 3);
				dis(i) = std::sqrt(std::pow((x0 + x1 + x2) * (x0 + x1 + x2), 1.0 / 3));
				break;
			case 4:
				x0 = std::pow(p(0) - _c(0), 4);
				x1 = std::pow(p(1) - _c(1), 4);
				x2 = std::pow(p(2) - _c(2), 4);
				dis(i) = std::pow(x0 + x1 + x2, 1.0 / 4);
				break;
			case 6:
				x0 = std::pow(p(0) - _c(0), 6);
				x1 = std::pow(p(1) - _c(1), 6);
				x2 = std::pow(p(2) - _c(2), 6);
				dis(i) = std::pow(x0 + x1 + x2, 1.0 / 6);
				break;
			case 8:
				x0 = std::pow(p(0) - _c(0), 8);
				x1 = std::pow(p(1) - _c(1), 8);
				x2 = std::pow(p(2) - _c(2), 8);
				dis(i) = std::pow(x0 + x1 + x2, 1.0 / 8);
				break;
			case 10:
				x0 = std::pow(p(0) - _c(0), 10);
				x1 = std::pow(p(1) - _c(1), 10);
				x2 = std::pow(p(2) - _c(2), 10);
				dis(i) = std::pow(x0 + x1 + x2, 1.0 / 10);
				break;
			case 16:
				x0 = std::pow(p(0) - _c(0), 16);
				x1 = std::pow(p(1) - _c(1), 16);
				x2 = std::pow(p(2) - _c(2), 16);
				dis(i) = std::pow(x0 + x1 + x2, 1.0 / 16);
				break;
			}
		}
		return dis;
	}
}

//namespace pqp_helper {
//	using namespace detail;
//	using namespace geometry;
//
//	inline void initPQPModel(const std::vector<Triangle<Vector3>>& modelTris, PQP_Model*& pqpModel)
//	{
//		pqpModel = new PQP_Model();
//		pqpModel->BeginModel();
//
//		for (int f_id = 0; f_id < modelTris.size(); f_id++)
//		{
//			PQP_REAL tri[3][3];
//			tri[0][0] = modelTris[f_id].p1(0), tri[0][1] = modelTris[f_id].p1(1), tri[0][2] = modelTris[f_id].p1(2);
//			tri[1][0] = modelTris[f_id].p2(0), tri[1][1] = modelTris[f_id].p2(1), tri[1][2] = modelTris[f_id].p2(2);
//			tri[2][0] = modelTris[f_id].p3(0), tri[2][1] = modelTris[f_id].p3(1), tri[2][2] = modelTris[f_id].p3(2);
//			pqpModel->AddTri(tri[0], tri[1], tri[2], f_id);
//		}
//
//		pqpModel->EndModel();
//	}
//
//	inline Scalar distanceFromPointToModel(const Vector3& queryPoint, PQP_Model* pqpModel)
//	{
//		PQP_DistanceResult res;
//		Scalar point[3] = { queryPoint(0), queryPoint(1), queryPoint(2) };
//		PQP_Distance(&res, pqpModel, point, 0.0, 0.0);
//		return res.distance;
//	}
//}

NAMESPACE_END(offset)