#include "PCOModel.hpp"
#include "detail/MortonLUT.hpp"
#include "utils/Log.hpp"
#include "utils/Timer.hpp"
#include "utils/Donut.hpp"
#include "utils/dis/DisHelper.hpp"
#include "complex/OffsetSurfaceBuilder.hpp"
#include "mesh/SurfaceMesh.hpp"

#ifdef ENABLE_VIEWER

#include "viewer/MeshViewer.hpp"

#endif

#include <omp.h>
#include <queue>
#include <numeric>
#include <utility>

#ifdef ERROR
#   undef ERROR
#endif

NAMESPACE_BEGIN(PCO)
namespace core {

	////////////////////////
	//    Constructors    //
	////////////////////////
	OffsetModel::OffsetModel(const std::string& filename,
		double offsetFactor_,
		int maxDepth_) noexcept
		: TriMesh(filename), offsetFactor(offsetFactor_), maxOctreeDepth(maxDepth_) {
		isOuterOffset = (offsetFactor > 0);
		offsetDis = (this->diagonalLengthOfBBox) * offsetFactor_ / 100.0;
		LOG::INFO("[PCO] Offset distance = {}.", offsetDis);
		abs_offsetDis = std::abs(offsetDis);

		setOffsetBoundingBox();

		octree = std::make_shared<Octree>(maxOctreeDepth);
		LOG::INFO("[PCO] Maximum octree depth = {}.", maxOctreeDepth);

		mergeFunc = [](TetDis& dis_1, const TetDis& dis_2) {
			dis_1 = dis_1.array().min(dis_2.array());
			};
		constMergeFunc = [](const TetDis& dis_1, const TetDis& dis_2) -> TetDis {
			return dis_1.array().min(dis_2.array());
			};
	}

	////////////////////////
	//   Initialization   //
	////////////////////////
	void
		OffsetModel::setOffsetBoundingBox() {
		Vector3 modelBBOrigin = (this->modelBoundingBox).boxOrigin;
		Vector3 modelBBWidth = (this->modelBoundingBox).boxWidth;
		Vector3 modelBBEnd = (this->modelBoundingBox).boxEnd;

		{
			Vector3 outerOffsetWidth =
				modelBBWidth + 3 * Vector3(abs_offsetDis, abs_offsetDis, abs_offsetDis);
			Vector3 outerOffsetBBOrigin =
				modelBBOrigin - 1.5 * Vector3(abs_offsetDis, abs_offsetDis, abs_offsetDis);
			Vector3 outerOffsetBBEnd =
				modelBBEnd + 1.5 * Vector3(abs_offsetDis, abs_offsetDis, abs_offsetDis);
			offsetBoundingBox = AABox<Vector3>(outerOffsetBBOrigin, outerOffsetBBEnd);
		}

		Vector3 minV = offsetBoundingBox.boxOrigin;
		Vector3 maxV = offsetBoundingBox.boxEnd;
		Vector3 lengths = maxV - minV; // check length of given bbox in every direction
		const double max_length = fmaxf(lengths.x(), fmaxf(lengths.y(), lengths.z())); // find max length
		for (unsigned int i = 0; i < 3; i++) { // for every direction (X,Y,Z)
			if (max_length == lengths[i]) {
				continue;
			}
			else {
				const double delta = max_length -
					lengths[i]; // compute difference between largest length and current (X,Y or Z) length
				offsetBoundingBox.boxOrigin[i] =
					minV[i] - (delta / 2.0f); // pad with half the difference before current min
				offsetBoundingBox.boxEnd[i] =
					maxV[i] + (delta / 2.0f); // pad with half the difference behind current max
			}
		}

		// Next snippet adresses the problem reported here: https://github.com/Forceflow/cuda_voxelizer/issues/7
		// Suspected cause: If a triangle is axis-aligned and lies perfectly on a voxel edge, it sometimes gets counted / not counted
		// Probably due to a numerical instability (division by zero?)
		// Ugly fix: we pad the bounding box on all sides by 1/10001th of its total length, bringing all triangles ever so slightly off-grid
		Vector3 epsilon = (offsetBoundingBox.boxEnd - offsetBoundingBox.boxOrigin) / 10001; // 之前是10001
		offsetBoundingBox.boxOrigin -= epsilon;
		offsetBoundingBox.boxEnd += epsilon;
		offsetBoundingBox.boxWidth = offsetBoundingBox.boxEnd - offsetBoundingBox.boxOrigin;
	}

	////////////////////////
	//  Core of Algorithm //
	////////////////////////
	int8_t
		OffsetModel::getFieldTypeInTet(const TetDis& dis) const {
		Scalar max_dis = donut::max(dis(0), dis(1), dis(2), dis(3));
		if (max_dis < abs_offsetDis) return -1;
		Scalar min_dis = donut::min(dis(0), dis(1), dis(2), dis(3));
		if (min_dis > abs_offsetDis) return 1;
		return 0;
	}

	void
		OffsetModel::createOctree(std::shared_ptr<OctreeNode>& node,
			const AABox<Vector3>& bbox, int depth) {
		node->depth = depth;
		node->bbox = bbox;
		node->center = bbox.center();
		node->width = bbox.boxWidth;
		node->circumRadius = 0.5 * bbox.boxWidth.norm();

		// set morton code
		int16_t x = node->xyzIdx(0);
		uint16_t y = node->xyzIdx(1);
		uint16_t z = node->xyzIdx(2);
		// order: zyx
		node->mortonCode = mortonEncode_LUT(x, y, z);

		// unsigned distance filter
		double circumRadius = node->circumRadius;
		node->centerDis = pointUDF(node->center);
		bool no_inter = (std::abs(node->centerDis - abs_offsetDis) > circumRadius);

		if (no_inter) {
			octree->depthNodes[depth].emplace_back(node);
			octree->allLeafNodes.emplace_back(node);
			return;
		}

#ifdef USE_SDF
        //double signed_centerDis = node->centerDis * pointSign(node->center); // use intersection test
		double signed_centerDis = pointPseudonormalSDF(node->center);
        /** or
         * double signed_centerDis = pointPseudonormalSDF(node->center);
         * double signed_centerDis = pointWnSDF(node->center);
         * */

        if (isOuterOffset) {
			if (signed_centerDis >= -abs_offsetDis - minNodeWidth - circumRadius &&
				signed_centerDis <= std::min(-abs_offsetDis + minNodeWidth + circumRadius,
					abs_offsetDis - minNodeWidth - circumRadius)) {
                octree->allLeafNodes.emplace_back(node);
                return;
            }
        } else {
            if (signed_centerDis >= std::max(abs_offsetDis - minNodeWidth - circumRadius,
				-abs_offsetDis + minNodeWidth + circumRadius) &&
				signed_centerDis <= abs_offsetDis + minNodeWidth + circumRadius) {
                octree->allLeafNodes.emplace_back(node);
				return;
            }
        }
#endif

		if (depth >= maxOctreeDepth - 1) {
			octree->depthNodes[depth].emplace_back(node);
			octree->validLeafNodes.emplace_back(node);
			octree->allLeafNodes.emplace_back(node);
			return;
		}

		Vector3 child_width = (node->width) * 0.5;
		node->is_leaf = false;
		for (int i = 0; i < 8; ++i) {
			node->child[i] = std::make_shared<OctreeNode>();
			node->child[i]->parent = node;
			node->child[i]->childIdx = i;

			//      0 - - - 2
			//     /       /
			//    /       /
			//   1 - - - 3

			Vector3 child_origin = bbox.boxOrigin;
			int local_x = (i & 1);
			int local_y = ((i >> 1) & 1);
			int local_z = ((i >> 2) & 1);
			child_origin.x() += child_width.x() * local_x;
			child_origin.y() += child_width.y() * local_y;
			child_origin.z() += child_width.z() * local_z;
			Vector3 child_end = child_origin + child_width;
			AABox<Vector3> child_bbox(child_origin, child_end);
			node->child[i]->xyzIdx = node->xyzIdx * 2 + Vector3i(local_x, local_y, local_z);

			createOctree(node->child[i], child_bbox, depth + 1);
		}
		octree->depthNodes[depth].emplace_back(node);
	}

	void
		OffsetModel::selectContributingTris(utils::timer::TimerInterface* algo_t) {
		size_t numValidLeafNodes = octree->validLeafNodes.size();
		LOG::INFO("[PCO] The number of valid leaf nodes = {}.", numValidLeafNodes);
		tetOffsetSurface.reserve(numValidLeafNodes * 5);
		double interRegion = std::sqrt(2.0) * abs_offsetDis;

		using namespace utils;

		timer::startTimer(&algo_t);
#pragma omp parallel for
		for (int i = 0; i < numValidLeafNodes; ++i) {
			auto& node = octree->validLeafNodes[i];
			node->setCorners();

			/// TODO: Optimize
			std::vector<CGAL_AABB::Primitive> intersection_primitives;
			AABox<Vector3> bbox = node->bbox;
			Point cubeMin = Point(bbox.boxOrigin(0) - interRegion, bbox.boxOrigin(1) - interRegion,
				bbox.boxOrigin(2) - interRegion);
			Point cubeMax = Point(bbox.boxEnd(0) + interRegion, bbox.boxEnd(1) + interRegion,
				bbox.boxEnd(2) + interRegion);
			Iso_cuboid cube(cubeMin, cubeMax);
			cgal_aabb.all_intersected_primitives(cube, std::back_inserter(intersection_primitives));
			LOG::DEBUG("[PCO] The number of intersected primitives = {}.", intersection_primitives.size());
			if (intersection_primitives.empty()) continue;

			int tetType = tetsTypeInNodeLUT[node->childIdx];
			for (int j = 0; j < 5; ++j) {
				OctreeTet octreeTet = node->getTet(j);
				Tet tet = octreeTet.coord;
				int tetLocalIdx = tetLocalIdxLUT[tetType][j];

				std::array<int, 4> verticesLocalId = octreeTet.vertsInLocalNode;
				std::array<size_t, 4> verticesGlobalId = octreeTet.globalVertsId;
				std::array<size_t, 6> edgesGlobalId = octreeTet.globalEdgesId;
				std::array<size_t, 4> facesGlobalId = octreeTet.globalFacesId;

				OffsetSurfaceBuilder builder(tetLocalIdx, tet, true,
					verticesGlobalId,
					edgesGlobalId,
					facesGlobalId);

				std::vector<std::pair<int, TetDis>> contriTrisAndDis;
				std::vector<std::pair<int, TetDis>> surviveTrisAndDis;

				size_t numPrimitives = intersection_primitives.size();
				size_t pri_i = 0;
				for (; pri_i < numPrimitives; ++pri_i) {
					const CGAL_AABB::Primitive& pri = intersection_primitives.at(pri_i);
					size_t triIdx = pri.id() - cgal_triangles.begin();

					TetDis dis = (utils::dishelper::simplexMultiSqrDis(tet, faceVertMats.at(triIdx)));
					dis = (dis.array().sqrt());

					int8_t disType = getFieldTypeInTet(dis);
					if (disType == -1) break;
					else {
						contriTrisAndDis.emplace_back(triIdx, dis);
						if (disType == 0) {
							builder.add_new_plane(dis.array() - abs_offsetDis, triIdx);

							if (builder.is_empty()) break;
							surviveTrisAndDis.emplace_back(triIdx, dis);
						}
					}
				}
				if (pri_i != numPrimitives) continue;

				builder.extract_offset_surface();
				if (builder.get_offset_surface().faces.empty()) continue;

#pragma omp critical
				{
					tetOffsetSurface.emplace_back(builder);

					int tetIdx = validTets.size();

					validTets.emplace_back(i, j, tet);
					validTets.back().tetType = tetType;

					for (const auto& tri_dis : contriTrisAndDis) {
						size_t triIdx = tri_dis.first;
						TetDis dis = tri_dis.second;

						LinearField lf(dis, tet);
						contriTetToTri[tetIdx][triIdx] = lf;
						contriTriToTet[triIdx][tetIdx] = lf;
					}

					for (const auto& face : builder.get_offset_surface().faces) {
						if (face.supporting_plane == OffsetSurface::None) continue;
						size_t planeIdx = face.supporting_plane - DIM - 1;
						size_t triIdx = surviveTrisAndDis[planeIdx].first;
						TetDis dis = surviveTrisAndDis[planeIdx].second;

						LinearField lf(dis, tet);
						surviveTetToTri[tetIdx][triIdx] = lf;
						surviveTriToTet[triIdx][tetIdx] = lf;
					}
				}
			}

		}

		timer::stopTimer(&algo_t);
		numValidTets = validTets.size();
	}

	void
		OffsetModel::buildCompatibleGraph(MDSU& mdsu,
			SparseGraphHandler*& G2Handler,
			const std::unordered_set<int>& updatedTets) {
		std::vector<SpGraphTriplet> G2Vec;

		std::unordered_set<int> visBaseFaTris;

		//#pragma omp parallel for
		for (int i = 0; i < numG2; ++i) {
			int baseFaTri = mdsu.find(G2IdxToTri.at(i));
			if (visBaseFaTris.count(baseFaTri)) continue;
			else visBaseFaTris.insert(baseFaTri);

			std::unordered_map<int, bool> isValidOtherTri;
			std::unordered_map<int, Scalar> alphaWithOtherTri;

			for (const auto& tet_field : surviveTriToTet.at(baseFaTri)) {
				int tetIdx = tet_field.first;

				if (!updatedTets.contains(tetIdx)) continue;
				if (!surviveTetToTri.contains(tetIdx)) continue;

				std::unordered_set<int> visOtherFaTris;
				const LinearField& baseField = tet_field.second;
				for (const auto& tri_field : surviveTetToTri.at(tetIdx)) {
					int otherFaTri = mdsu.find(tri_field.first);

					if (otherFaTri == baseFaTri) continue;

					if (visOtherFaTris.contains(otherFaTri)) continue;
					else visOtherFaTris.insert(otherFaTri);

					if (isValidOtherTri.contains(otherFaTri) && !isValidOtherTri.at(otherFaTri)) continue;

					const LinearField& other_field = tri_field.second;

					if (isIntersect(baseField, other_field, abs_offsetDis)) {
						TetDis mergeField = constMergeFunc(baseField.field, other_field.field);
						if (getFieldTypeInTet(mergeField) == 0) {
							Scalar sub_cos_alpha = gradientAngle(baseField, other_field);
							if (sub_cos_alpha <= COS_COMP_ANGLE) {
								isValidOtherTri[otherFaTri] = true;
								alphaWithOtherTri[otherFaTri] += (COS_COMP_ANGLE - sub_cos_alpha + 1e-6);
							}
							else isValidOtherTri[otherFaTri] = false;
						}
						else isValidOtherTri[otherFaTri] = false;
					}
					else isValidOtherTri[otherFaTri] = false;
				}
			}

			int baseFaTriIdx = triToG2Idx.at(baseFaTri);
			for (auto& tri_valid : isValidOtherTri) {
				int otherFaTri = tri_valid.first;
				if (tri_valid.second) {
					int otherFaTriIdx = triToG2Idx.at(otherFaTri);

					// 再检查contri中的四面体是否可以合并，保证拓扑正确
					for (const auto& tet_field : contriTriToTet.at(baseFaTri)) {
						// 得到 base_tri 贡献的所有四面体
						int tetIdx = tet_field.first; // 0~numValidTets, 注意并非tet的实际索引
						if (!surviveTetToTri.contains(tetIdx)) continue;
						if (!contriTetToTri.at(tetIdx).contains(otherFaTri)) continue;

						auto builder = tetOffsetSurface.at(tetIdx);
						const LinearField& baseField = tet_field.second;
						const LinearField& otherField = contriTetToTri.at(tetIdx).at(otherFaTri);
						TetDis mergeField = constMergeFunc(baseField.field, otherField.field);
						builder.add_new_plane(mergeField.array() - offsetDis, baseFaTri);
						if (builder.is_empty()) {
							tri_valid.second = false;
							break;
						}
					}
					if (!tri_valid.second) continue;

					Scalar weight = alphaWithOtherTri.at(otherFaTri);
					// TODO: optimize construction
					G2Vec.emplace_back(baseFaTriIdx, otherFaTriIdx, weight);
				}
			}
		}

		SpGraph G2(numG2, numG2);
		G2.setFromTriplets(G2Vec.begin(), G2Vec.end());
		G2.makeCompressed();

		G2Handler = new SparseGraphAdaptor(G2);
		G2Handler->updateEdges();
	}

	TetDis
		OffsetModel::getDSUField(int baseFaTri,
			int tetIdx,
			const MDSU& mdsu) {
		TetDis dsuField;
		if (!contriTriToTet.at(baseFaTri).contains(tetIdx)) {
			TetDis dis = utils::dishelper::simplexMultiSqrDis(validTets[tetIdx].tet,
				faceVertMats.at(baseFaTri));
			dis = (dis.array().sqrt());

			dsuField = dis;

			int cnt = 0;
			std::unordered_set<int> children = mdsu.getChildren(baseFaTri);
			for (int child : children) {
				if (child == baseFaTri) continue;
				++cnt;
				if (!contriTriToTet.at(child).contains(tetIdx)) {
					dis = utils::dishelper::simplexMultiSqrDis(validTets[tetIdx].tet,
						faceVertMats.at(child));
					dis = dis.array().sqrt();
					TetDis childField = dis;

					TetDis _dsuField = dsuField;
					mergeFunc(dsuField, childField);
				}
				else {
					mergeFunc(dsuField, contriTriToTet.at(child).at(tetIdx).field); // 可以不需要这行
				}
			}
		}
		else {
			dsuField = contriTriToTet.at(baseFaTri).at(tetIdx).field;
		}

		return dsuField;
	}

	std::unordered_set<int>
		OffsetModel::mergeNodesByClique(SparseGraphHandler* GHandler, MDSU& mdsu) {
		std::unordered_set<int> updatedTets;
		std::vector<std::list<int>> cliques = GHandler->run_maximalCliques("degeneracy");

		using PAIR_VAL = std::pair<int, Scalar>;
		std::priority_queue<IndexedValue<PAIR_VAL>, std::vector<IndexedValue<PAIR_VAL>>,
			decltype(&IndexedValue<PAIR_VAL>::asc_cmp)>
			maxHeap(IndexedValue<PAIR_VAL>::asc_cmp);
		std::vector<IndexedValue<PAIR_VAL>>
			nodePariVal;
		std::unordered_map<int, std::unordered_set<int>> nodeToClique;

		nodePariVal.reserve(cliques.size());
		for (int i = 0; i < cliques.size(); ++i) {
			Scalar edgeWeightSum = .0;
			std::list<int> clique = cliques[i];
			for (auto outer_it = clique.begin(); outer_it != clique.end(); ++outer_it) {
				nodeToClique[*outer_it].insert(i);
				auto inner_it = std::next(outer_it);
				for (; inner_it != clique.end(); ++inner_it) {
					edgeWeightSum += GHandler->G.coeff(*outer_it, *inner_it);
				}
			}
			nodePariVal.emplace_back(i, std::make_pair(clique.size(), -edgeWeightSum));
			maxHeap.emplace(i, std::make_pair(clique.size(), -edgeWeightSum));
		}

		std::set<IndexedValue<PAIR_VAL>> oldData;
		while (!maxHeap.empty()) {
			auto fro = maxHeap.top();
			maxHeap.pop();

			if (oldData.contains(fro)) continue;
			if (fro.value.first <= 1) continue;

			std::list<int> cur_clique = cliques[fro.index];
			std::unordered_set<int> curCliqueNodes;
			std::unordered_set<int> relatedCliques;

			auto baseNodeIter = cur_clique.begin();
			int baseNode = *baseNodeIter;
			int baseFaTri = mdsu.find(G2IdxToTri.at(baseNode));
			std::unordered_set<int> cliqueUpdatedTets;
			for (int node : cur_clique) {
				curCliqueNodes.insert(node);
				relatedCliques.insert(nodeToClique.at(node).begin(), nodeToClique.at(node).end());

				int otherFaTri = mdsu.find(G2IdxToTri.at(node));
				if (otherFaTri == baseFaTri) continue;

				mdsu.merge(baseFaTri, otherFaTri);
				for (const auto& tet_dis : contriTriToTet.at(G2IdxToTri.at(node))) {
					cliqueUpdatedTets.insert(tet_dis.first);
				}
			}
			updatedTets.insert(cliqueUpdatedTets.begin(), cliqueUpdatedTets.end());

			for (int tetIdx : cliqueUpdatedTets) {
				TetDis baseDSUField = getDSUField(baseFaTri, tetIdx, mdsu);

				for (auto otherNodeIter = std::next(baseNodeIter);
					otherNodeIter != cur_clique.end(); ++otherNodeIter) {
					int otherTri = G2IdxToTri.at(*otherNodeIter);

					TetDis otherDSUField = getDSUField(otherTri, tetIdx, mdsu);
					TetDis _baseDSUField = baseDSUField;
					/*if (getFieldTypeInTet(baseDSUField) == -1)
						LOG::ERROR("getFieldTypeInTet(baseDSUField) == -1");*/
					TetDis mergeDis = constMergeFunc(baseDSUField, otherDSUField);
					if (getFieldTypeInTet(mergeDis) == -1) continue;

					mergeFunc(baseDSUField, otherDSUField);
				}
				contriTriToTet[baseFaTri][tetIdx] = LinearField(baseDSUField, validTets[tetIdx].tet);
				contriTetToTri[tetIdx][baseFaTri] = LinearField(baseDSUField, validTets[tetIdx].tet);

				surviveTriToTet[baseFaTri][tetIdx] = LinearField(baseDSUField, validTets[tetIdx].tet);
				surviveTetToTri[tetIdx][baseFaTri] = LinearField(baseDSUField, validTets[tetIdx].tet);
			}

			for (int relatedClique : relatedCliques) {
				if (relatedClique == fro.index) continue;
				oldData.insert(nodePariVal[relatedClique]);

				std::list<int>& otherClique = cliques[relatedClique];

				for (int node : curCliqueNodes) {
					if (std::find(otherClique.begin(), otherClique.end(), node) !=
						otherClique.end()) {
						otherClique.remove(node);
					}
				}

				Scalar newEdgeWeightSum = .0;
				for (auto outer_it = otherClique.begin(); outer_it != otherClique.end(); ++outer_it) {
					auto inner_it = std::next(outer_it);
					for (; inner_it != otherClique.end(); ++inner_it) {
						newEdgeWeightSum += GHandler->G.coeff(*outer_it, *inner_it);
					}
				}

				nodePariVal[relatedClique].value.first = otherClique.size();
				nodePariVal[relatedClique].value.second = -newEdgeWeightSum;
				maxHeap.emplace(nodePariVal[relatedClique]);
			}

			GHandler->softDeleteMultiNodes(curCliqueNodes);
		}
		return updatedTets;
	}

	int
		OffsetModel::fieldsMerge(double COMP_ANGLE) {
		LOG::INFO("Compatible angle = {} degree", COMP_ANGLE);

		COS_COMP_ANGLE = std::cos(COMP_ANGLE * PI / 180.0);

		MDSU mdsu;
		std::unordered_set<int> updatedTets;

		triToG2Idx.reserve(surviveTriToTet.size());
		G2IdxToTri.reserve(surviveTriToTet.size());
		for (const auto& tri_tet : surviveTriToTet) {
			int tri = tri_tet.first;
			mdsu.fa[tri] = tri;
			triToG2Idx[tri] = numG2;
			G2IdxToTri[numG2++] = tri;
		}
		//std::cout << "numG2: " << numG2 << std::endl;
		for (const auto& tet_tri : surviveTetToTri)
			updatedTets.insert(tet_tri.first);

		SparseGraphHandler* G2Handler = new SparseGraphAdaptor();
		buildCompatibleGraph(mdsu, G2Handler, updatedTets);

		int maxMergeIter = 10;
		int mergeIter = 0;
		bool isNeedMerge = (G2Handler->numEdges() != 0);
		do {
			if (isNeedMerge) {
				++mergeIter;
				LOG::INFO("Merge Iter = {}.", mergeIter);
				LOG::INFO("-- The number of nodes in G2 = {}.", G2Handler->numNodes());
				LOG::INFO("-- The number of edges in G2 = {}.", G2Handler->numEdges());
				LOG::INFO("-- The degeneracy of G2 = {}.", G2Handler->getDegeneracy());

				updatedTets.clear();
				updatedTets = mergeNodesByClique(G2Handler, mdsu);

				for (int tetIdx : updatedTets) {
					if (!surviveTetToTri.contains(tetIdx)) continue;
					auto& builder = tetOffsetSurface[tetIdx];

					std::unordered_set<int> visFaTris;
					std::unordered_map<size_t, std::pair<int, TetDis>> candidateFields;

					for (const auto& tri_field : contriTetToTri.at(tetIdx)) {
						int fa_tri = mdsu.find(tri_field.first);
						if (fa_tri == -1) fa_tri = tri_field.first;
						if (visFaTris.count(fa_tri)) continue;
						visFaTris.insert(fa_tri);
						TetDis field = contriTetToTri.at(tetIdx).at(fa_tri).field;
						if (getFieldTypeInTet(field) == -1) {
							LOG::ERROR("error: getFieldTypeInTet(field) == -1!!\n");
						}

						size_t planeIdx = builder.add_new_plane(field.array() - abs_offsetDis, fa_tri);
						candidateFields[planeIdx] = std::make_pair(fa_tri, field);
					}

					builder.extract_offset_surface();
					if (builder.get_offset_surface().faces.empty()) {
						surviveTetToTri.erase(tetIdx);
						contriTetToTri.erase(tetIdx); // 是否要加？
						continue;
					}
					surviveTetToTri.at(tetIdx).clear();
					for (const auto& face : builder.get_offset_surface().faces) {
						if (face.supporting_plane == OffsetSurface::None) continue;
						int planeIdx = face.supporting_plane;

						if (!candidateFields.contains(planeIdx)) continue;
						int triIdx = candidateFields.at(planeIdx).first;

						TetDis dis = candidateFields[planeIdx].second;

						LinearField lf(dis, validTets[tetIdx].tet);
						if (mdsu.find(triIdx) == -1) {
							mdsu.fa[triIdx] = triIdx;
							triToG2Idx[triIdx] = numG2;
							G2IdxToTri[numG2++] = triIdx;
						}
						surviveTetToTri[tetIdx][triIdx] = lf;
						surviveTriToTet[triIdx][tetIdx] = lf;
					}
				}
			}

			buildCompatibleGraph(mdsu, G2Handler, updatedTets);

			isNeedMerge = (G2Handler->numEdges() != 0);
		} while (mergeIter < maxMergeIter && isNeedMerge);

		return mergeIter;
	}

	////////////////////////
	//   APIs for User    //
	////////////////////////
	void
		OffsetModel::launch(const Args& args) {
		using namespace utils;

		int step = 1;

		timer::TimerInterface* algo_t;
		timer::createTimer(&algo_t);

		timer::startTimer(&algo_t);
        minNodeWidth = offsetBoundingBox.diagLength() / (1 << (maxOctreeDepth - 1));
		createOctree(octree->root, offsetBoundingBox, 0);
		timer::stopTimer(&algo_t);
		LOG::INFO("[PCO] STEP[{}] - Create Octree spent {} s.",
			step++, timer::getElapsedTime(&algo_t) / 1000.0);

        timer::startTimer(&algo_t);
        octree->constructAllNodes();
        octree->constructNodePrimitives();
        octree->constructNodePrimitives();
        timer::stopTimer(&algo_t);
        LOG::INFO("[PCO] STEP[{}] - Create Octree Info spent {} s.",
                  step++, timer::getElapsedTime(&algo_t) / 1000.0);

		timer::startTimer(&algo_t);
		octree->constructTetPrimitives();
		timer::stopTimer(&algo_t);
		LOG::INFO("[PCO] STEP[{}] - Create Tets spent {} s.",
			step++, timer::getElapsedTime(&algo_t) / 1000.0);

		timer::startTimer(&algo_t);
		selectContributingTris(algo_t);
		timer::stopTimer(&algo_t);
		LOG::INFO("[PCO] STEP[{}] - Select contributing triangles and extract offset surface spent {} s.",
			step++, timer::getElapsedTime(&algo_t) / 1000.0);

		if (args.isMerge) {
			timer::startTimer(&algo_t);
			fieldsMerge(args.COMP_ANGLE);
			timer::stopTimer(&algo_t);
			LOG::INFO("[PCO] STEP[{}] - Fields merging spent {} s.",
				step++, timer::getElapsedTime(&algo_t) / 1000.0);
		}

		// use single floating type
		timer::startTimer(&algo_t);
		std::vector<Vector3f> vertVec;
		std::vector<SurfaceMesh::Polygon> faceVec;
		outputOffsetSurface(vertVec, faceVec);
		SurfaceMesh surfaceMesh(vertVec, faceVec);
		surfaceMesh.processOffsetMesh();
		timer::stopTimer(&algo_t);
		LOG::INFO("[PCO] STEP[{}] - Process the result spent {} s.",
			step++, timer::getElapsedTime(&algo_t) / 1000.0);

		if (!args.outFile.empty()) {
			timer::startTimer(&algo_t);
			surfaceMesh.writeMesh(args.outFile);
			timer::stopTimer(&algo_t);
			LOG::INFO("[PCO] STEP[{}] - Output result spent {} s.",
				step++, timer::getElapsedTime(&algo_t) / 1000.0);
		}

		LOG::INFO("[PCO] The algorithm spent {} s.", timer::getAllTimeValue(&algo_t) / 1000.0);
		deleteTimer(&algo_t);

#ifdef ENABLE_VIEWER
		if (args.enableView) {
			LOG::INFO("[PCO] Initialize viewer..."),
				viewer::MeshViewer(modelName, this,
					utils::file::getFileName(args.outFile), &surfaceMesh);
		}
#endif
	}

	////////////////////////
	//        Output      //
	////////////////////////
	void
		OffsetModel::outputOffsetSurface(std::vector<Vector3f>& vertVec,
			std::vector<std::vector<size_t>>& faceVec) const {
		using SPlanes = std::array<size_t, 3>; // global supporting planes
		using SPlanesEdge = std::array<size_t, 2>; // global supporting plane and edge idx

		size_t vertCnt = 0;

		/// Note that even if the three planes forming a vertex are distinct,
		/// the coordinates may still coincide due to numerical errors.
		/// So we perform an additional process (in 'launch()') to remove duplicate vertices and faces.

		std::unordered_map<SPlanes, size_t> faceVertToOutIdx; // 2 outer planes + 1 inner plane
		std::unordered_map<SPlanesEdge, size_t> edgeVertToOutIdx; // 1 outer plane + 2 inner planes
		std::unordered_map<size_t, size_t> tetVertToOutIdx; // 3 inner planes
		std::set<std::vector<size_t>> finalFace;

		for (auto offsetSurfaceBuilder : tetOffsetSurface) {
			const auto& complexCut = offsetSurfaceBuilder.get_complex();
			const auto& coplanar_planes = offsetSurfaceBuilder.get_coplanar_planes();
			const auto& global_extPlaneId = offsetSurfaceBuilder.get_global_extPlaneId();
			const auto& offsetSurface = offsetSurfaceBuilder.get_offset_surface();
			const auto& offsetSurfaceFaces = offsetSurface.faces;
			const auto& offsetSurfaceVertices = offsetSurface.vertices;
			const auto& offsetSurfaceVerticesCoord = offsetSurface.vertices_coord;

			std::unordered_map<size_t, size_t> innerVertToOutIdx; // first key is the local id of vertex (3 outter planes)
			for (const auto& f : offsetSurfaceFaces) {
				if (f.supporting_plane == OffsetSurface::None) continue;
				std::vector<size_t> face;

				for (int k = 0; k < f.vertices.size(); ++k) {
					size_t iv = f.vertices[k];
					const auto& v = offsetSurfaceVertices.at(iv);

					std::vector<size_t> local_oriPlaneId;
					std::vector<size_t> local_extPlaneId;
					// If the supporting plane is a face of the tetrahedron, it remains unchanged even if coinciding with an external plane.
					// The union-find parent node is always the ID of the last added coinciding external plane.
					(v[0] <= 3) ? local_oriPlaneId.push_back(v[0]) : local_extPlaneId.push_back(v[0]);
					(v[1] <= 3) ? local_oriPlaneId.push_back(v[1]) : local_extPlaneId.push_back(v[1]);
					(v[2] <= 3) ? local_oriPlaneId.push_back(v[2]) : local_extPlaneId.push_back(v[2]);
					int numOriPlane = local_oriPlaneId.size();

					if (numOriPlane == 0) {
						if (!innerVertToOutIdx.count(iv)) {
							innerVertToOutIdx[iv] = vertCnt++;
							vertVec.emplace_back(offsetSurfaceVerticesCoord.at(iv).cast<float>());
						}
						else {
						}
						face.emplace_back(innerVertToOutIdx.at(iv));
					}
					else if (numOriPlane == 1) {
						// First, obtain the global IDs of the two external faces, including the one coinciding with the supporting plane.
						size_t local_extFaceId_0 = coplanar_planes.find(local_extPlaneId[0]);
						size_t local_extFaceId_1 = coplanar_planes.find(local_extPlaneId[1]);
						size_t global_extFaceId_0 = global_extPlaneId.at(local_extFaceId_0 - 4);
						size_t global_extFaceId_1 = global_extPlaneId.at(local_extFaceId_1 - 4);

						// Next, retrieve the global ID of the internal face.
						size_t global_oriFaceId = complexCut.global_oriFacesId[local_oriPlaneId[0]];

						if (global_extFaceId_0 > global_extFaceId_1)
							std::swap(global_extFaceId_0, global_extFaceId_1);
						SPlanes support_planes = { global_extFaceId_0, global_extFaceId_1, global_oriFaceId };
						if (!faceVertToOutIdx.count(support_planes)) {
							faceVertToOutIdx[support_planes] = vertCnt++;
							vertVec.emplace_back(offsetSurfaceVerticesCoord.at(iv).cast<float>());
						}
						else {
						}
						face.emplace_back(faceVertToOutIdx.at(support_planes));
					}
					else if (numOriPlane == 2) {
						// First, obtain the global ID of the external face, including the one coinciding with the supporting plane.
						size_t local_extFaceId = coplanar_planes.find(local_extPlaneId[0]);
						size_t global_extFaceId = global_extPlaneId.at(local_extFaceId - 4);

						// Next, retrieve the local ID of the shared edge between the two internal faces.
						const auto& plane_edges_0 = complexCut.oriFacesEdges[local_oriPlaneId[0]];
						const auto& plane_edges_1 = complexCut.oriFacesEdges[local_oriPlaneId[1]];

						int local_edgeId = -1;
						for (int i = 0; i < 3; ++i) {
							for (int j = 0; j < 3; ++j) {
								if (plane_edges_0[i] == plane_edges_1[j]) {
									local_edgeId = plane_edges_0[i];
									break;
								}
							}
							if (local_edgeId != -1) break;
						}
						// Then, obtain the global ID of the shared edge.
						size_t global_edgeId = complexCut.global_oriEdgesId[local_edgeId];

						SPlanesEdge support_plane_edge = { global_extFaceId, global_edgeId };
						if (!edgeVertToOutIdx.count(support_plane_edge)) {
							edgeVertToOutIdx[support_plane_edge] = vertCnt++;
							vertVec.emplace_back(offsetSurfaceVerticesCoord.at(iv).cast<float>());
						}
						else {
						}
						face.emplace_back(edgeVertToOutIdx.at(support_plane_edge));
					}
					else {
						size_t global_verticesId = f.global_verticesId.at(k);
						if (!tetVertToOutIdx.count(global_verticesId)) {
							tetVertToOutIdx[global_verticesId] = vertCnt++;
							vertVec.emplace_back(offsetSurfaceVerticesCoord.at(iv).cast<float>());
						}
						else {
						}
						face.emplace_back(tetVertToOutIdx.at(global_verticesId));
					}
				}

				auto sorted_face = face;
				std::sort(sorted_face.begin(), sorted_face.end());
				if (!finalFace.contains(sorted_face)) {
					finalFace.insert(sorted_face);
#ifdef USE_SDF
					if (!isOuterOffset) std::reverse(face.begin(), face.end());
#endif // USE_SDF
					faceVec.emplace_back(face);
				}
			}
		}
	}

	void
		OffsetModel::outputOffsetSurfaceFast(std::vector<Vector3f>& vertVec,
			std::vector<std::vector<size_t>>& faceVec) const {
		for (auto offsetSurfaceBuilder : tetOffsetSurface) {
			const auto& offsetSurface = offsetSurfaceBuilder.get_offset_surface();
			const auto& offsetSurfaceFaces = offsetSurface.faces;
			const auto& offsetSurfaceVertices = offsetSurface.vertices;
			const auto& offsetSurfaceVerticesCoord = offsetSurface.vertices_coord;

			for (const auto& f : offsetSurfaceFaces) {
				if (f.supporting_plane == OffsetSurface::None) continue;
				std::vector<size_t> face;
				for (int k = 0; k < f.vertices.size(); ++k) {
					size_t iv = f.vertices[k];
					face.push_back(vertVec.size());
					vertVec.push_back(offsetSurfaceVerticesCoord.at(iv).cast<float>());
				}
				faceVec.push_back(face);
			}
		}
	}

} // namespace core
NAMESPACE_END(PCO)