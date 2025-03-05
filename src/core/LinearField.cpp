//
// Created by Lei on 3/30/2024.
//

#include "LinearField.hpp"
#include "utils/Log.hpp"
#include <utils/Donut.hpp>

NAMESPACE_BEGIN(PCO)
namespace core {

	static constexpr int CROSS_EDGE_IDX[7][4] = { {0, 1, 2, -1},
													 {0, 3, 4, -1},
													 {1, 3, 5, -1},
													 {2, 4, 5, -1},
													 {1, 2, 3, 4},
													 {0, 2, 3, 5},
													 {0, 1, 4, 5} };

	static constexpr int TET_EDGE_PAIR[6][2] = { {0, 1},
													{0, 2},
													{0, 3},
													{1, 2},
													{1, 3},
													{2, 3} };

	inline int
		getFieldTypeInTet2(const TetDis& dis, double offsetDis) {

		int type = -1;
		if (!donut::isSandwitchVal(offsetDis, dis(0), dis(1), dis(2), dis(3))) return type;
		if (!donut::isSandwitchVal(offsetDis, dis(0), dis(1), dis(2))) {
			type = 3;
			if (dis(3) < offsetDis) type += 7;
		}
		else if (!donut::isSandwitchVal(offsetDis, dis(0), dis(1), dis(3))) {
			type = 2;
			if (dis(2) < offsetDis) type += 7;
		}
		else if (!donut::isSandwitchVal(offsetDis, dis(0), dis(2), dis(3))) {
			type = 1;
			if (dis(1) < offsetDis) type += 7;
		}
		else if (!donut::isSandwitchVal(offsetDis, dis(1), dis(2), dis(3))) {
			type = 0;
			if (dis(0) < offsetDis) type += 7;
		}

		if (type == -1) {
			type = 4;
			for (int i = 1; i <= 3; ++i) {
				if (!donut::isSandwitchVal(offsetDis, dis(0), dis(i))) {
					break;
				}
				++type;
			}
			if (dis(0) < offsetDis) type += 7;
		}
		return type;
	}


	inline
		std::vector<Vector3> getEdgePoints(const TetDis& dis, const Tet& tetQ, double offsetDis) {
		std::vector<Vector3> edgePoints;
		int disType = getFieldTypeInTet2(dis, offsetDis) % 7;
		for (int cross_edge : CROSS_EDGE_IDX[disType]) {
			if (cross_edge == -1) break;
			const auto& [vert_0, vert_1] = TET_EDGE_PAIR[cross_edge];
			double t = (offsetDis - dis(vert_0)) / (dis(vert_1) - dis(vert_0));
			Vector3 edgePoint = tetQ.row(vert_0) + t * (tetQ.row(vert_1) - tetQ.row(vert_0));
			edgePoints.emplace_back(edgePoint);
		}
		return edgePoints;
	}

	bool isIntersect(const LinearField& lf_1, const LinearField& lf_2, double offsetDis) {
		MatrixX A(3, 4);
		A.row(0) = lf_1.field;
		A.row(1) = lf_2.field;
		A.row(2) = TetDis::Ones();

		Eigen::FullPivLU<MatrixX> lu_decomp_A(A);
		LOG::DEBUG("A.rank() = {}.", lu_decomp_A.rank());

		MatrixX A_(3, 5);
		A_.block<3, 4>(0, 0) = A;
		A_.col(3) = Vector4(offsetDis, offsetDis, offsetDis, 1);
		Eigen::FullPivLU<MatrixX> lu_decomp_A_(A_);
		LOG::DEBUG("A_.rank() = {}.", lu_decomp_A_.rank());

		return (lu_decomp_A.rank() == lu_decomp_A_.rank());
	}

	bool isIntersect(const LinearField& lf_1, const LinearField& lf_2,
		const Tet& tet, double offsetDis, int tetType) {
		std::vector<Vector3> edgePoints = getEdgePoints(lf_1.field, tet, offsetDis);
		const VectorX other_param = lf_2.gradient;
		Vector3 edgePoint = edgePoints[0];
		Vector4 homo_coord = Vector4(edgePoint.x(), edgePoint.y(), edgePoint.z(), 1.0);

		if (tetType == 0) {
			bool isUpper = (other_param.dot(homo_coord)) > 0;
			bool isOn = (other_param.dot(homo_coord)) == 0;

			int numEdgePoints = edgePoints.size();
			for (int i = 1; i < numEdgePoints; ++i) {
				edgePoint = edgePoints[i];
				homo_coord = Vector4(edgePoint.x(), edgePoint.y(), edgePoint.z(), 1.0);
				// 等于也算相交(除非当前距离场完全与other这个距离场重合)
				if (!isOn &&
					(isUpper ^ ((other_param.dot(homo_coord)) > 0) || other_param.dot(homo_coord) == 0))
					return true;
				else if (isOn && (other_param.dot(homo_coord) != 0)) return true;
			}
			return false;
		}
		else
		{
			bool isLower = (other_param.dot(homo_coord)) < 0;
			bool isOn = (other_param.dot(homo_coord)) == 0;

			int numEdgePoints = edgePoints.size();
			for (int i = 1; i < numEdgePoints; ++i) {
				edgePoint = edgePoints[i];
				homo_coord = Vector4(edgePoint.x(), edgePoint.y(), edgePoint.z(), 1.0);
				if (!isOn &&
					(isLower ^ ((other_param.dot(homo_coord)) > 0) || other_param.dot(homo_coord) == 0))
					return true;
				else if (isOn && (other_param.dot(homo_coord) != 0)) return true;
			}
			return false;
		}
	}

	bool isHaveValidUnion(const LinearField& lf_1, const LinearField& lf_2, double offsetDis) {
		bool flag1 = false;
		for (int i = 0; i < 4; ++i) {
			if ((lf_1.field(i) > offsetDis) ==
				(lf_2.field(i) > offsetDis)) {
				flag1 = true;
				break;
			}
		}
		if (!flag1) return false;

		bool flag2 = false;
		for (int i = 0; i < 4; ++i) {
			if ((lf_1.field(i) < offsetDis) ==
				(lf_2.field(i) < offsetDis)) {
				flag2 = true;
				break;
			}
		}

		return flag2;
	}

	bool isHaveValidUnion(const LinearField& lf_1, const LinearField& lf_2, const Tet& tetQ) {
		MatrixX _tetQ(4, 4);
		if (tetQ.cols() == 3) {
			_tetQ.block(0, 0, 4, 3) = tetQ;
			//                    _tetQ.col(3) = Vector4::Ones();
			_tetQ.col(3).setOnes();
		}
		else _tetQ = tetQ;
		const VectorX param = lf_1.gradient;
		const VectorX other_param = lf_2.gradient;
		for (int i = 0; i < tetQ.rows(); ++i) {
			Vector4 homo_coord = _tetQ.row(i);
			if (homo_coord.dot(param) > 0 && homo_coord.dot(other_param) > 0) return true;
		}
		return false;
	}

	Scalar gradientAngle(const LinearField& lf_1, const LinearField& lf_2) {
		/*std::cout << "gradient#1: " << lf_1.gradient.transpose() << std::endl;
		std::cout << "gradient#2: " << lf_2.gradient.transpose() << std::endl;
		std::cout << "===========\n" << std::endl;*/
		return -lf_1.gradient.dot(lf_2.gradient);
	}

} // namespace core
NAMESPACE_END(PCO)