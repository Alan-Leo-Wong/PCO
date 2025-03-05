//
// Created by Lei on 3/31/2024.
//

#include "ComplexCut.hpp"
#include "utils/Log.hpp"
#include "detail/Handler.hpp"

#include <absl/container/flat_hash_map.h>
#include <implicit_predicates/implicit_predicates.h>

NAMESPACE_BEGIN(PCO)
namespace complex {
	template<typename T, typename KeepFunc>
	std::vector<size_t> shrink(std::vector<T>& c, const KeepFunc& keep) {
		const size_t s = c.size();
		std::vector<size_t> index_map(s, INVALID);
		size_t active_count = 0;
		for (size_t i = 0; i < s; i++) {
			if (!keep(i)) continue;

			if (i != active_count) {
				std::swap(c[active_count], c[i]);
			}
			index_map[i] = active_count;
			active_count++;
		}
		c.resize(active_count);

		return index_map;
	}

	int8_t signof(implicit_predicates::Orientation o) {
		PCO_ASSERT(o != implicit_predicates::INVALID);
		if (o > 0) return 1;
		else if (o < 0) return -1;
		else return 0;
	};

	int8_t ComplexCut::cut_0_face(size_t vid,
		size_t plane_index) {
		const auto& p = planes.get_plane(plane_index); // 当前添加的平面
		const auto& v = complex.vertices[vid];  // 当前考虑的顶点, 由 p0 p1 (p2) 构成
		const auto& p0 = planes.get_plane(v[0]);       // v[index]: 每个顶点使用其关联的平面 index 表示
		const auto& p1 = planes.get_plane(v[1]);
		const auto& p2 = planes.get_plane(v[2]);
#ifdef NON_ROBUST_COMPLEX_CUTTING
		return implicit_predicates::orient3d_nonrobust(p0.data(), p1.data(), p2.data(), p.data());
#else
		return implicit_predicates::orient3d(p0.data(), p1.data(), p2.data(), p.data());
#endif // NON_ROBUST_COMPLEX_CUTTING
	}

	std::array<size_t, 3> ComplexCut::cut_1_face(size_t eid,
		size_t plane_index,
		const std::vector<int8_t>& orientations) {
		auto& vertices = complex.vertices;
		auto& edges = complex.edges;
		vertices.reserve(vertices.size() + 1);
		edges.reserve(edges.size() + 2);

		size_t positive_subedge_id = INVALID;
		size_t negative_subedge_id = INVALID;
		size_t intersection_id = INVALID;

		const auto& e = edges[eid];
		const auto& end_points = e.vertices;
		const auto o0 = orientations[end_points[0]];
		const auto o1 = orientations[end_points[1]];
		LOG::DEBUG("edge {} has end points ({}, {}) with orientations {} {}",
			eid,
			end_points[0],
			end_points[1],
			o0,
			o1);

		if (o0 == 0) intersection_id = end_points[0];
		if (o1 == 0) intersection_id = end_points[1];

		auto compute_intersection_id = [&]() {
			size_t p0 = e.supporting_planes[0];
			size_t p1 = e.supporting_planes[1];
			LOG::DEBUG("Adding cut vertex: {}, {}, {}", p0, p1, plane_index);
			vertices.push_back({ p0, p1, plane_index });
			return vertices.size() - 1;
			};

		if (o0 == 0 && o1 == 0) {
			// Collinear!
			LOG::DEBUG("Edge {} is coplanar with plane {}", eid, plane_index);
		}
		else if (o0 >= 0 && o1 >= 0) {
			positive_subedge_id = eid;
		}
		else if (o0 <= 0 && o1 <= 0) {
			negative_subedge_id = eid;
		}
		else {
			PCO_ASSERT(intersection_id == INVALID);
			intersection_id = compute_intersection_id();
			Edge positive_subedge, negative_subedge;

			if (o0 > 0 && o1 < 0) {
				positive_subedge.vertices = { end_points[0], intersection_id };
				negative_subedge.vertices = { intersection_id, end_points[1] };
			}
			else {
				PCO_ASSERT(o0 < 0);
				PCO_ASSERT(o1 > 0);
				negative_subedge.vertices = { end_points[0], intersection_id };
				positive_subedge.vertices = { intersection_id, end_points[1] };
			}

			positive_subedge.supporting_planes = e.supporting_planes;
			negative_subedge.supporting_planes = e.supporting_planes;

			edges.push_back(std::move(positive_subedge));
			edges.push_back(std::move(negative_subedge));
			positive_subedge_id = edges.size() - 2;
			negative_subedge_id = edges.size() - 1;
			LOG::DEBUG("Adding positive subedge: {}", positive_subedge_id);
			LOG::DEBUG("Adding negative subedge: {}", negative_subedge_id);
		}
		LOG::DEBUG("Cut edge {}: {} {} with intersection {}",
			eid,
			positive_subedge_id,
			negative_subedge_id,
			intersection_id);
		return { positive_subedge_id, negative_subedge_id, intersection_id };
	}

	std::array<size_t, 3> ComplexCut::cut_2_face(size_t fid,
		size_t plane_index,
		const std::vector<int8_t>& orientations,
		const std::vector<std::array<size_t, 3>>& subedges) {
		auto& edges = complex.edges;
		auto& faces = complex.faces;
		edges.reserve(edges.size() + 1);
		faces.reserve(faces.size() + 2);

		auto& f = faces[fid];
		const size_t num_bd_edges = f.edges.size();
		LOG::DEBUG("face {} has {} bd edges", fid, num_bd_edges);

		std::vector<size_t> positive_subedges, negative_subedges;
		positive_subedges.reserve(num_bd_edges);
		negative_subedges.reserve(num_bd_edges);

		Edge cut_edge;
		size_t cut_edge_index = INVALID;
		bool face_is_coplanar = true;
		size_t cut_edge_positive_location = INVALID;
		size_t cut_edge_negative_location = INVALID;

		auto get_end_vertex = [&](size_t local_eid) {
			size_t curr_eid = f.edges[local_eid];
			size_t next_eid = f.edges[(local_eid + 1) % num_bd_edges];
			const auto& e0 = edges[curr_eid];
			const auto& e1 = edges[next_eid];
			if (e0.vertices[0] == e1.vertices[0] || e0.vertices[0] == e1.vertices[1]) {
				return e0.vertices[0];
			}
			else {
				PCO_ASSERT(e0.vertices[1] == e1.vertices[0] || e0.vertices[1] == e1.vertices[1]);
				return e0.vertices[1];
			}
			};

		for (size_t j = 0; j < num_bd_edges; j++) {
			const size_t eid = f.edges[j];

			bool last_positive = false;
			bool last_negative = false;
			size_t intersection_point = subedges[eid][2];
			size_t positive_subedge = subedges[eid][0];
			size_t negative_subedge = subedges[eid][1];
			LOG::DEBUG("{}: intersection point {}", eid, intersection_point);

			if (positive_subedge == INVALID && negative_subedge == INVALID) {
				// Edge is coplanar with the plane.
				cut_edge_index = eid;
			}
			else {
				if (positive_subedge != INVALID) {
					positive_subedges.push_back(positive_subedge);
					size_t end_vid = get_end_vertex(j);
					if (orientations[end_vid] <= 0) {
						// This edge is the last edge with positive subedge.
						cut_edge_positive_location = positive_subedges.size();
						last_positive = true;
					}
				}
				if (negative_subedge != INVALID) {
					negative_subedges.push_back(negative_subedge);
					size_t end_vid = get_end_vertex(j);
					if (orientations[end_vid] >= 0) {
						// This edge is the last edge with negative subedge.
						cut_edge_negative_location = negative_subedges.size();
						last_negative = true;
					}
				}
				face_is_coplanar = false;
			}
			if (intersection_point != INVALID) {
				if (last_positive) {
					cut_edge.vertices[0] = intersection_point;
				}
				else if (last_negative) {
					cut_edge.vertices[1] = intersection_point;
				}
			}
		}
		LOG::DEBUG("num positive subedges: {}", positive_subedges.size());
		LOG::DEBUG("num negative subedges: {}", negative_subedges.size());
		LOG::DEBUG("cut point: {}, {}", cut_edge.vertices[0], cut_edge.vertices[1]);

		if (face_is_coplanar) {
			LOG::DEBUG("Face {} is coplanar with cut plane", fid);
			return { INVALID, INVALID, INVALID };
		}

		if (positive_subedges.empty() || negative_subedges.empty()) {
			// No cut.
			if (positive_subedges.empty()) {
				PCO_ASSERT(!negative_subedges.empty());
				return { INVALID, fid, cut_edge_index };
			}
			else {
				PCO_ASSERT(!positive_subedges.empty());
				PCO_ASSERT(negative_subedges.empty());
				return { fid, INVALID, cut_edge_index };
			}
		}

		PCO_ASSERT(cut_edge_index == INVALID);
		{
			// Insert cut edge.
			cut_edge.supporting_planes = { f.supporting_plane, plane_index };
			cut_edge_index = edges.size();
			edges.emplace_back(std::move(cut_edge));
			LOG::DEBUG("Adding cut edge: {}", cut_edge_index);
		}

		// Create subfaces.
		Face positive_subface, negative_subface;
		positive_subface.supporting_plane = f.supporting_plane;
		negative_subface.supporting_plane = f.supporting_plane;
		positive_subface.positive_cell = f.positive_cell;
		positive_subface.negative_cell = f.negative_cell;
		negative_subface.positive_cell = f.positive_cell;
		negative_subface.negative_cell = f.negative_cell;

		PCO_ASSERT(cut_edge_index != INVALID);
		if (cut_edge_positive_location != positive_subedges.size()) {
			std::rotate(positive_subedges.begin(),
				positive_subedges.begin() + cut_edge_positive_location,
				positive_subedges.end());
		}
		if (cut_edge_negative_location != negative_subedges.size()) {
			std::rotate(negative_subedges.begin(),
				negative_subedges.begin() + cut_edge_negative_location,
				negative_subedges.end());
		}
		positive_subedges.push_back(cut_edge_index);
		positive_subface.edges = std::move(positive_subedges);
		PCO_ASSERT(positive_subface.edges.size() > 2);
		negative_subedges.push_back(cut_edge_index);
		negative_subface.edges = std::move(negative_subedges);
		PCO_ASSERT(negative_subface.edges.size() > 2);

		faces.emplace_back(std::move(positive_subface));
		faces.emplace_back(std::move(negative_subface));
		size_t positive_fid = faces.size() - 2;
		size_t negative_fid = faces.size() - 1;
		LOG::DEBUG("Adding positive subface: {}", positive_fid);
		LOG::DEBUG("Adding negative subface: {}", negative_fid);

		return { positive_fid, negative_fid, cut_edge_index };
	}

	std::array<size_t, 3> ComplexCut::cut_3_face(size_t cid,
		size_t plane_index,
		const std::vector<std::array<size_t, 3>>& subfaces) {
		auto& edges = complex.edges;
		auto& faces = complex.faces;
		auto& cells = complex.cells;

		const auto& c = cells[cid];
		const size_t num_bd_faces = c.faces.size();
		LOG::DEBUG("cell {} has {} bd faces", cid, num_bd_faces);

		size_t cut_face_id = INVALID;
		std::vector<size_t> positive_subfaces;
		std::vector<size_t> negative_subfaces;
		std::vector<size_t> cut_edges;
		std::vector<bool> cut_edge_orientations;
		positive_subfaces.reserve(num_bd_faces + 1);
		negative_subfaces.reserve(num_bd_faces + 1);
		cut_edges.reserve(num_bd_faces);
		cut_edge_orientations.reserve(num_bd_faces);

		auto compute_cut_edge_orientation = [&](size_t fid,
			const std::array<size_t, 3>& subface) -> bool {
				PCO_ASSERT(subface[2] != INVALID);
				const auto& f = faces[fid];
				bool s = c.signs[f.supporting_plane];

				if (subface[0] == INVALID || subface[1] == INVALID) {
					// Intersection edge is on the boundary of the face.
					auto itr = std::find(f.edges.begin(), f.edges.end(), subface[2]);
					PCO_ASSERT(itr != f.edges.end());
					size_t curr_i = itr - f.edges.begin();
					size_t next_i = (curr_i + 1) % f.edges.size();

					const auto& curr_e = edges[f.edges[curr_i]];
					const auto& next_e = edges[f.edges[next_i]];
					bool edge_is_consistent_with_face = (curr_e.vertices[1] == next_e.vertices[0] ||
						curr_e.vertices[1] == next_e.vertices[1]);

					bool on_positive_side = subface[0] != INVALID;
					LOG::DEBUG("cell/face: {}, face/edge: {}, face/cut plane: {}",
						s,
						edge_is_consistent_with_face,
						on_positive_side);

					uint8_t key = 0;
					if (s) key++;
					if (edge_is_consistent_with_face) key++;
					if (on_positive_side) key++;
					return key % 2 == 0;
				}
				else {
					// Intersection edge is a cross cut.
					//std::cout << "fid: " << fid << ", s: " <<std::boolalpha << s << std::endl;
					return !s;
				}
			};

		for (auto fid : c.faces) {
			const auto& subface = subfaces[fid];
			// plane与面重叠
			if (subface[0] == INVALID && subface[1] == INVALID) {
				cut_face_id = fid;
			}
			if (subface[0] != INVALID) {
				positive_subfaces.push_back(subface[0]);
			}
			if (subface[1] != INVALID) {
				negative_subfaces.push_back(subface[1]);
			}
			if (subface[2] != INVALID) {
				cut_edges.push_back(subface[2]);
				cut_edge_orientations.push_back(compute_cut_edge_orientation(fid, subface));
			}
		}

		if (positive_subfaces.empty() && negative_subfaces.empty()) {
			// The implicit function is identical over the whole cell.
			return { INVALID, INVALID, INVALID };
		}
		else if (positive_subfaces.empty()) {
			cells[cid].signs[plane_index] = false;
			return { INVALID, cid, cut_face_id };
		}
		else if (negative_subfaces.empty()) {
			cells[cid].signs[plane_index] = true;
			return { cid, INVALID, cut_face_id };
		}

		// Chain cut edges into a loop.
		{
			size_t num_cut_edges = cut_edges.size();
			PCO_ASSERT(num_cut_edges >= 3);
			absl::flat_hash_map<size_t, size_t> v2e;
			v2e.reserve(num_cut_edges);
			for (size_t i = 0; i < num_cut_edges; i++) {
				const auto eid = cut_edges[i];
				const auto& e = edges[eid];
				if (cut_edge_orientations[i]) {
					v2e[e.vertices[0]] = i;
				}
				else {
					v2e[e.vertices[1]] = i;
				}
			}
			std::vector<size_t> chained_cut_edges;
			chained_cut_edges.reserve(num_cut_edges);
			chained_cut_edges.push_back(0);
			while (chained_cut_edges.size() < num_cut_edges) {
				const size_t i = chained_cut_edges.back();
				const auto& e = edges[cut_edges[i]];
				const size_t vid = cut_edge_orientations[i] ? e.vertices[1] : e.vertices[0];
				const auto itr = v2e.find(vid);
				PCO_ASSERT(itr != v2e.end());
				const size_t next_i = itr->second;
				if (cut_edges[next_i] == cut_edges[chained_cut_edges.front()]) {
					break;
				}
				chained_cut_edges.push_back(next_i);
			}
			std::transform(chained_cut_edges.begin(),
				chained_cut_edges.end(),
				chained_cut_edges.begin(),
				[&](size_t i) { return cut_edges[i]; });
			std::swap(cut_edges, chained_cut_edges);
		}

		// Cross cut.
		PCO_ASSERT(!cut_edges.empty());
		Face cut_face;
		cut_face.edges = std::move(cut_edges);
		cut_face.supporting_plane = plane_index;
		faces.push_back(std::move(cut_face));
		cut_face_id = faces.size() - 1;
		LOG::DEBUG("Adding cut face: {}", cut_face_id);

		// Generate positive and negative subcell.
		Cell positive_cell, negative_cell;
		positive_cell.faces.reserve(positive_subfaces.size() + 1);
		negative_cell.faces.reserve(negative_subfaces.size() + 1);

		positive_subfaces.push_back(cut_face_id);
		positive_cell.faces = std::move(positive_subfaces);
		positive_cell.signs = c.signs;
		positive_cell.signs[plane_index] = true;

		negative_subfaces.push_back(cut_face_id);
		negative_cell.faces = std::move(negative_subfaces);
		negative_cell.signs = c.signs;
		negative_cell.signs[plane_index] = false;

		cells.emplace_back(std::move(positive_cell));
		cells.emplace_back(std::move(negative_cell));
		size_t positive_cell_id = cells.size() - 2;
		size_t negative_cell_id = cells.size() - 1;
		LOG::DEBUG("Adding positive subcell: {}", positive_cell_id);
		LOG::DEBUG("Adding negative subcell: {}", negative_cell_id);

		// Update cell id on each side of involved faces.
		{
			// cut face
			PCO_ASSERT(cut_face_id != INVALID);
			auto& cut_f = faces[cut_face_id];
			cut_f.positive_cell = positive_cell_id;
			cut_f.negative_cell = negative_cell_id;

			auto& positive_c = cells[positive_cell_id];
			auto& negative_c = cells[negative_cell_id];

			for (auto fid : positive_c.faces) {
				if (fid == cut_face_id) continue;
				auto& f = faces[fid];
				PCO_ASSERT(f.positive_cell == cid || f.negative_cell == cid);
				if (f.positive_cell == cid) {
					f.positive_cell = positive_cell_id;
				}
				else {
					f.negative_cell = positive_cell_id;
				}
			}
			for (auto fid : negative_c.faces) {
				if (fid == cut_face_id) continue;
				auto& f = faces[fid];
				PCO_ASSERT(f.positive_cell == cid || f.negative_cell == cid);
				if (f.positive_cell == cid) {
					f.positive_cell = negative_cell_id;
				}
				else {
					f.negative_cell = negative_cell_id;
				}
			}
		}

		return { positive_cell_id, negative_cell_id, cut_face_id };
	}

	void ComplexCut::consolidate() {
		// Shrink faces.
		{
			std::vector<bool> active_faces(complex.faces.size(), false);
			for (auto& c : complex.cells) {
				for (auto fid : c.faces) {
					active_faces[fid] = true;
				}
			}

			auto index_map = shrink(complex.faces, [&](size_t fid) { return active_faces[fid]; });

			for (auto& c : complex.cells) {
				std::transform(c.faces.begin(), c.faces.end(), c.faces.begin(), [&](size_t i) {
					PCO_ASSERT(index_map[i] != INVALID);
					return index_map[i];
					});
			}
		}

		// Shrink edges.
		{
			std::vector<bool> active_edges(complex.edges.size(), false);
			for (auto& f : complex.faces) {
				for (auto eid : f.edges) {
					active_edges[eid] = true;
				}
			}

			auto index_map = shrink(complex.edges, [&](size_t eid) { return active_edges[eid]; });

			for (auto& f : complex.faces) {
				std::transform(f.edges.begin(), f.edges.end(), f.edges.begin(), [&](size_t i) {
					PCO_ASSERT(index_map[i] != INVALID);
					return index_map[i];
					});
			}
		}

		// Shrink vertices.
		{
			std::vector<bool> active_vertices(complex.vertices.size(), false);
			for (auto& e : complex.edges) {
				for (auto vid : e.vertices) {
					PCO_ASSERT(vid != INVALID);
					active_vertices[vid] = true;
				}
			}

			auto index_map =
				shrink(complex.vertices, [&](size_t vid) { return active_vertices[vid]; });

			for (auto& e : complex.edges) {
				e.vertices[0] = index_map[e.vertices[0]];
				e.vertices[1] = index_map[e.vertices[1]];
			}
		}
	}


	size_t ComplexCut::add_plane(size_t plane_index,
		bool is_positive) {
		const size_t num_vertices = complex.vertices.size();
		const size_t num_edges = complex.edges.size();
		const size_t num_faces = complex.faces.size();
		const size_t num_cells = complex.cells.size();
		LOG::DEBUG("adding plane {}", plane_index);
		LOG::DEBUG("Before: {} {} {} {}", num_vertices, num_edges, num_faces, num_cells);

		auto& vertices = complex.vertices;
		auto& edges = complex.edges;
		auto& faces = complex.faces;
		auto& cells = complex.cells;

		// Reserve capacity. TODO: check if this helps.
		vertices.reserve(num_vertices + num_edges);

		// Step 1: handle 0-faces.
		std::vector<int8_t> orientations;
		orientations.reserve(num_vertices);
		bool valid_flag = false;
		for (size_t i = 0; i < num_vertices; i++) {
			const auto& v = vertices[i];
			orientations.push_back(cut_0_face(i, plane_index));
			if (!valid_flag) {
				if (is_positive && (orientations.back() == -1 || orientations.back() == 0)) valid_flag = true;
				else if (!is_positive && (orientations.back() == 1 || orientations.back() == 0)) valid_flag = true;
			}
		}
		if (!valid_flag) return INVALID;

		edges.reserve(num_edges * 2);
		faces.reserve(num_faces * 2);
		cells.reserve(num_cells * 2);

		// Step 2: handle 1-faces.
		std::vector<std::array<size_t, 3>> subedges;
		subedges.reserve(num_edges);
		for (size_t i = 0; i < num_edges; i++) {
			subedges.push_back(cut_1_face(i, plane_index, orientations));
		}

		// Step 3: handle 2-faces.
		std::vector<std::array<size_t, 3>> subfaces;
		subfaces.reserve(num_faces);
		for (size_t i = 0; i < num_faces; i++) {
			subfaces.push_back(cut_2_face(i, plane_index, orientations, subedges));
		}

		// Step 4: handle 3-faces.
		std::vector<std::array<size_t, 3>> subcells;
		subcells.reserve(num_cells);
		for (size_t i = 0; i < num_cells; i++) {
			// return: positive_cell, negative_cell, cut_face_id
			subcells.push_back(cut_3_face(i, plane_index, subfaces));
		}

		// Step 5: remove old cells and update cell indices
		{
			std::vector<bool> to_keep(cells.size(), false);
			for (const auto& subcell : subcells) {
				if (is_positive) {
					if (subcell[0] != INVALID) to_keep[subcell[0]] = true;
				}
				else {
					if (subcell[1] != INVALID) to_keep[subcell[1]] = true;
				}
			}

			auto index_map = shrink(cells, [&](size_t i) { return to_keep[i]; });

			// Update cell indices in faces.
			for (auto& f : faces) {
				if (is_positive) {
					if (f.positive_cell != INVALID) f.positive_cell = index_map[f.positive_cell];
				}
				else {
					if (f.negative_cell != INVALID) f.negative_cell = index_map[f.negative_cell];
				}
			}
		}

		// Step 6: check for coplanar planes.
		size_t coplanar_plane = INVALID;
		for (size_t i = 0; i < num_faces; i++) {
			const auto& subface = subfaces[i];
			if (subface[0] == INVALID && subface[1] == INVALID) {
				const auto& f = faces[i];
				coplanar_plane = f.supporting_plane;
			}
		}

		// Step 7: consolidate.
		consolidate();
		LOG::DEBUG("After: {} {} {} {}",
			complex.vertices.size(),
			complex.edges.size(),
			complex.faces.size(),
			complex.cells.size());

		return coplanar_plane;
	}

	size_t ComplexCut::add_plane(const Plane& plane, bool is_positive) {
		complex.cells[0].signs.push_back(false);
		planes.add_new_plane(plane);
		return add_plane(planes.get_num_planes() + DIM, is_positive);
	}

} // namespace complex
NAMESPACE_END(PCO)