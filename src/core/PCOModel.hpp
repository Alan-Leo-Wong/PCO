#pragma once

#include "Args.hpp"
#include "LinearField.hpp"
#include "mesh/TriMesh.hpp"
#include "utils/dsu/DSU.hpp"
#include "utils/Timer.hpp"
#include "detail/Octree.hpp"
#include "graph/GraphAdaptor.hpp"
#include "complex/OffsetSurfaceBuilder.hpp"

#include <set>
#include <utility>
#include <functional>
#include <unordered_set>
#include <unordered_map>

class TimerInterface;
NAMESPACE_BEGIN(PCO)
namespace core {
	using namespace utils::dsu;
	using namespace geometry;
	using namespace complex;
	using namespace graph;
	using namespace mesh;

	class OffsetModel : public TriMesh {
	public:
		using MergeFunc = std::function<void(TetDis&, const TetDis&)>;
		using ConstMergeFunc = std::function<TetDis(const TetDis&, const TetDis&)>;
		using TetCornerSign = Eigen::Array<int, 4, 1>;

		using SpGraph = SparseGraphHandler::SpGraph;
		using SpGraphTriplet = SparseGraphHandler::SpGraphTrip;
		using IdxToDis = std::unordered_map<size_t, LinearField>;

		struct PCOTet {
			int nodeIdx, localTetIdx;
			int tetType;
			Tet tet;

			PCOTet() = default;

			PCOTet(int _nodeIdx, int _localTetIdx, Tet _tet) :
				nodeIdx(_nodeIdx), localTetIdx(_localTetIdx), tet(std::move(_tet)) {}
		};

	public:
		/* Constructors and Destructors */
		OffsetModel() noexcept = default;

		OffsetModel(const std::string& filename, double offsetFactor_, int maxDepth_) noexcept;

		~OffsetModel() noexcept override = default;

	private:
		/**
		 * Set the bounding-box of offsetting surface.
		 */
		void setOffsetBoundingBox();

	public:
		/**
		 * [API] Set the offset distance by user.
		 * @param _offsetDis
		 */
		void setOffsetDis(double _offsetDis) noexcept {
			offsetDis = _offsetDis;
			abs_offsetDis = std::abs(_offsetDis);
			setOffsetBoundingBox();
			LOG::INFO("[PCO] Offset distance = ", offsetDis);
		}

		[[nodiscard]] double getOffsetDis() const noexcept { return offsetDis; }

		[[nodiscard]] AABox<Vector3> getOffsetBB() const noexcept { return offsetBoundingBox; }

	private:
		/**
		 * Compute type of the linear field in a tetrahedron.
		 * @param dis
		 * @return
		 */
		int8_t getFieldTypeInTet(const TetDis& dis) const;

	private:
		/**
		 * Build an adaptive octree, where each volume will be subdivided into 5 tets.
		 * @param node
		 * @param bbox
		 * @param depth
		 */
		void createOctree(std::shared_ptr<OctreeNode>& node, const AABox<Vector3>& bbox, int depth);

		/**
		 * Calculate the local distance field contributed by each tet.
		 */
		void selectContributingTris(utils::timer::TimerInterface* algo_t);

		/**
		 * Merge the undirected graph according to the minimum edge weight clique
		 * (cliques sorted by size from largest to smallest).
		 * @param GHandler
		 * @param mdsu
		 * @param G2IdxToTri
		 * @return
		 */
		std::unordered_set<int>
			mergeNodesByClique(SparseGraphHandler* GHandler, MDSU& mdsu);

		/**
		 * Obtain the local distance field for a tet belonging to the same disjoint-set collection.
		 * @param baseFaTri
		 * @param tetIdx
		 * @param mdsu
		 * @return
		 */
		TetDis getDSUField(int baseFaTri,
			int tetIdx,
			const MDSU& mdsu);

		/**
		 * Create the compatibility graph G2 between triangles.
		 * @param mdsu
		 * @param G2Handler
		 * @param updatedTets
		 */
		void buildCompatibleGraph(MDSU& mdsu,
			SparseGraphHandler*& G2Handler,
			const std::unordered_set<int>& updatedTets);

		/**
		 * Merge compatible local distance fields.
		 * @param COMP_ANGLE Angle threshold.
		 * @return Iterations of merging process.
		 */
		int fieldsMerge(double COMP_ANGLE);

	public:
		/**
		 * [API] Launch the algorithm for user.
		 * @param args
		 * @return
		 */
		void launch(const Args& args);

	public:

		/**
		 * [API] Output offset surface to two vectors.
		 * @param vertVec the output offset surface's vertices
		 * @param faceVec vertex index of each face
		 */
		void outputOffsetSurface(std::vector<Vector3f>& vertVec,
			std::vector<std::vector<size_t>>& faceVec) const;

        void outputOffsetSurfaceFast(std::vector<Vector3f> &vertVec,
            std::vector<std::vector<size_t>> &faceVec) const;

		/**
		 * [API] Output constructed offsetting surface in .obj.
		 * @param filename
		 * @return
		 */
		bool outputOffsetSurface(const std::string& filename) const;

	private:
        bool isOuterOffset;
		double offsetDis;
		double offsetFactor;
		double abs_offsetDis; // absolute value of offset

		AABox<Vector3> offsetBoundingBox; // Bounding-box of offsetting surface

	private:
		int maxOctreeDepth;
        double minNodeWidth;
		std::shared_ptr<Octree> octree;

	private:
		size_t numValidTets;
		std::vector<PCOTet> validTets;

		std::vector<OffsetSurfaceBuilder> tetOffsetSurface;

	private:
		MergeFunc mergeFunc;
		ConstMergeFunc constMergeFunc;
		double COS_COMP_ANGLE;

		// 关系图 G1
		std::unordered_map<size_t, IdxToDis> contriTetToTri; // 每个四面体索引(非id)到对其有贡献的三角形的映射, 值为距离
		std::unordered_map<size_t, IdxToDis> contriTriToTet; // 每个三角形到对其有贡献的四面体索引(非id)的映射, 值为距离

		size_t numG2 = 0;
		// 每次增量合并时所需要的数据结构
		std::unordered_map<int, int> triFa; // 并查集每个三角形(集合)的fa数组
		std::unordered_map<int, std::unordered_set<int>> triFaToTetSet;
		std::unordered_map<size_t, IdxToDis> surviveTetToTri; // 保存四面体索引(非id)到竞争成功的三角形(fa)的映射
		std::unordered_map<size_t, IdxToDis> surviveTriToTet; // 保存竞争成功的三角形(fa)到四面体索引(非id)的映射

		std::unordered_map<int, int> triToG2Idx;
		std::unordered_map<int, int> G2IdxToTri;

	private:
		std::string out_dir;

	private:
		/* Winding number query related */
		static constexpr int FWN_ACCURACY_SACLE = 6;
		static constexpr double FWN_INNER_EPSLION = 1e-5;
		static constexpr double FWN_OUTER_EPSLION = 1e-9;
	};

} // namespace core
NAMESPACE_END(PCO)