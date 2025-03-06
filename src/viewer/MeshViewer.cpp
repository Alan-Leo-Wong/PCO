#include "MeshViewer.hpp"
#include "FileDialog.hpp"

#include <omp.h>
#include <polyscope/pick.h>
#include <polyscope/volume_grid.h>
#include <polyscope/point_cloud.h>
#include <polyscope/implicit_helpers.h>
#include <polyscope/surface_vector_quantity.h>

#ifdef ERROR
#   undef ERROR
#endif

NAMESPACE_BEGIN(PCO)
namespace viewer {

	namespace internal {
		using pSurfaceMesh = MeshViewer::pSurfaceMesh;
		using pCurveNetwork = MeshViewer::pCurveNetwork;

		TriMesh* ourTriMesh = nullptr;
		SurfaceMesh* ourSurfaceMesh = nullptr;

		pSurfaceMesh* psTriMesh = nullptr, * psSurfaceMesh = nullptr;
		pCurveNetwork* triMeshCurveNetwork = nullptr, * surfaceMeshCurveNetwork = nullptr;

		MeshData triMeshData, surfaceMeshData;

		void setVertexScalarQuantity(MeshData& meshData, pSurfaceMesh* psMesh, bool isReg) {
			const size_t nVertices = meshData.nVertices;
			std::vector<double> valX(nVertices);
			std::vector<double> valY(nVertices);
			std::vector<double> valZ(nVertices);
			std::vector<double> valMag(nVertices);
			std::vector<std::array<double, 3>> randColor(nVertices);
			for (size_t iV = 0; iV < nVertices; iV++) {
				valX[iV] = meshData.vertexPositionsGLM[iV].x / 10000;
				valY[iV] = meshData.vertexPositionsGLM[iV].y;
				valZ[iV] = meshData.vertexPositionsGLM[iV].z;
				valMag[iV] = glm::length(meshData.vertexPositionsGLM[iV]);

				randColor[iV] = { {polyscope::randomUnit(), polyscope::randomUnit(), polyscope::randomUnit()} };
			}

			if (isReg) {
				meshData.vertexScalarQuantities["cX"] = psMesh->addVertexScalarQuantity("cX", valX);
				meshData.vertexScalarQuantities["cY"] = psMesh->addVertexScalarQuantity("cY", valY);
				meshData.vertexScalarQuantities["cZ"] = psMesh->addVertexScalarQuantity("cZ", valZ);

				meshData.vertexColorQuantities["cZ"] = psMesh->addVertexColorQuantity("vColor", randColor);

				meshData.vertexScalarQuantities["cY_sym"] = psMesh->addVertexScalarQuantity("cY_sym", valY,
					polyscope::DataType::SYMMETRIC);
				meshData.vertexScalarQuantities["cNorm"] = psMesh->addVertexScalarQuantity("cNorm", valMag,
					polyscope::DataType::MAGNITUDE);

				meshData.vertexScalarQuantities["cY_dist"] = psMesh->addVertexDistanceQuantity("cY_dist", valY);
				meshData.vertexScalarQuantities["cY_signeddist"] = psMesh->addVertexSignedDistanceQuantity(
					"cY_signeddist", valY);
			}
			else {
				meshData.vertexScalarQuantities["cX"]->updateData(valX);
				meshData.vertexScalarQuantities["cY"]->updateData(valY);
				meshData.vertexScalarQuantities["cZ"]->updateData(valZ);

				meshData.vertexColorQuantities["cZ"]->updateData(randColor);

				meshData.vertexScalarQuantities["cY_sym"]->updateData(valY);
				meshData.vertexScalarQuantities["cNorm"]->updateData(valMag);

				meshData.vertexScalarQuantities["cY_dist"]->updateData(valY);
				meshData.vertexScalarQuantities["cY_signeddist"]->updateData(valY);
			}
		}

		void setFaceScalarQuantity(MeshData& meshData, pSurfaceMesh* psMesh, bool isReg) {
			const size_t nFaces = meshData.nFaces;
			std::vector<double> fArea(nFaces);
			std::vector<std::array<double, 3>> fColor(nFaces);
			for (size_t iF = 0; iF < nFaces; iF++) {
				const std::vector<size_t>& face = meshData.faceIndices[iF];

				// Compute something like area
				double area = 0;
				for (size_t iV = 1; iV < face.size() - 1; iV++) {
					glm::vec3 p0 = meshData.vertexPositionsGLM[face[0]];
					glm::vec3 p1 = meshData.vertexPositionsGLM[face[iV]];
					glm::vec3 p2 = meshData.vertexPositionsGLM[face[iV + 1]];
					area += 0.5f * glm::length(glm::cross(p1 - p0, p2 - p0));
				}
				fArea[iF] = area;
				fColor[iF] = { {polyscope::randomUnit(), polyscope::randomUnit(), polyscope::randomUnit()} };
			}
			if (isReg) {
				meshData.faceScalarQuantities["face area"] = psMesh->addFaceScalarQuantity("face area", fArea,
					polyscope::DataType::MAGNITUDE);
				meshData.faceColorQuantities["fColor"] = psMesh->addFaceColorQuantity("fColor", fColor);
			}
			else {
				meshData.faceScalarQuantities["face area"]->updateData(fArea);
				meshData.faceColorQuantities["fColor"]->updateData(fColor);
			}
		}

		void setEdgeScalarQuantity(MeshData& meshData, pSurfaceMesh* psMesh, bool isReg) {
			const size_t nFaces = meshData.nFaces;
			std::vector<double> eLen;
			std::vector<double> heLen;
			std::vector<double> cAngle;
			std::vector<uint32_t> edgeOrdering;
			for (size_t iF = 0; iF < nFaces; iF++) {
				const std::vector<size_t>& face = meshData.faceIndices[iF];

				for (size_t iV = 0; iV < face.size(); iV++) {
					size_t i0 = face[iV];
					size_t i1 = face[(iV + 1) % face.size()];
					size_t im1 = face[(iV + face.size() - 1) % face.size()];
					glm::vec3 p0 = meshData.vertexPositionsGLM[i0];
					glm::vec3 p1 = meshData.vertexPositionsGLM[i1];
					glm::vec3 pm1 = meshData.vertexPositionsGLM[im1];

					double len = glm::length(p0 - p1);

					double angle = glm::acos(glm::dot(glm::normalize(p1 - p0), glm::normalize(pm1 - p0)));

					size_t iMin = std::min(i0, i1);
					size_t iMax = std::max(i0, i1);

					auto p = std::make_pair(iMin, iMax);
					if (meshData.seenEdges.find(p) == meshData.seenEdges.end()) {
						eLen.push_back(len);
						edgeOrdering.push_back(
							edgeOrdering.size()); // totally coincidentally, this is the trivial ordering
						meshData.seenEdges.insert(p);
					}
					heLen.push_back(len);
					cAngle.push_back(angle);
				}
			}
			meshData.nEdges = edgeOrdering.size();

			if (isReg) {
				psMesh->setEdgePermutation(edgeOrdering);
				meshData.edgeScalarQuantities["edge length"] = psMesh->addEdgeScalarQuantity("edge length", eLen);
				meshData.halfEdgeScalarQuantities["halfedge length"] = psMesh->addHalfedgeScalarQuantity(
					"halfedge length", heLen);
				meshData.cornerScalarQuantities["corner angle"] = psMesh->addCornerScalarQuantity("corner angle",
					cAngle);
			}
			else {
				meshData.halfEdgeScalarQuantities["halfedge length"]->updateData(heLen);
				meshData.cornerScalarQuantities["corner angle"]->updateData(cAngle);
			}

		}

		void setNormalQuantity(MeshData& meshData, pSurfaceMesh* psMesh, bool isReg) {
			const size_t nVertices = meshData.nVertices;
			const auto& fNormals = meshData.fNormals;
			const auto& vNormals = meshData.vNormals;
			std::vector<int> vertexNeighbors(nVertices, 0);
			std::vector<double> normalDeviation(nVertices, .0);

			for (size_t iF = 0; iF < meshData.faceIndices.size(); ++iF) {
				const auto& fNormal = fNormals.at(iF);
				for (size_t v = 0; v < meshData.faceIndices.at(iF).size(); ++v) {
					size_t iV = meshData.faceIndices.at(iF).at(v);
					const auto& vNormal = vNormals.at(iV);
					normalDeviation[iV] += glm::dot(vNormal, fNormal);
					vertexNeighbors[iV]++;
				}
			}

			std::vector<glm::vec3> vNormalsRand(nVertices, glm::vec3{ 0., 0., 0. });
			std::vector<glm::vec3> toZero(nVertices, glm::vec3{ 0., 0., 0. });
			for (size_t iV = 0; iV < nVertices; iV++) {
				vNormalsRand[iV] = vNormals[iV] * (float)polyscope::randomUnit() * 5000.f;
				toZero[iV] = -meshData.vertexPositionsGLM[iV];
				normalDeviation[iV] /= vertexNeighbors[iV];
			}

			if (isReg) {
				meshData.faceVectorQuantities["face normals"] = psMesh->addFaceVectorQuantity("face normals", fNormals);
				meshData.vertexVectorQuantities["area vertex normals"] = psMesh->addVertexVectorQuantity(
					"area vertex normals", vNormals);
				meshData.vertexVectorQuantities["rand length normals"] = psMesh->addVertexVectorQuantity(
					"rand length normals", vNormalsRand);
				meshData.vertexVectorQuantities["toZero"] = psMesh->addVertexVectorQuantity("toZero", toZero,
					polyscope::VectorType::AMBIENT);
				meshData.vertexScalarQuantities["normal deviation"] = psMesh->addVertexScalarQuantity("normal deviation", normalDeviation);
			}
			else {
				meshData.vertexVectorQuantities["area vertex normals"]->updateData(vNormals);
				meshData.vertexVectorQuantities["rand length normals"]->updateData(vNormalsRand);
				meshData.vertexVectorQuantities["toZero"]->updateData(toZero);
				meshData.vertexScalarQuantities["normal deviation"]->updateData(normalDeviation);
			}
		}

		void setIntrisicVectorQuantity(MeshData& meshData, pSurfaceMesh* psMesh, bool isReg) {
			const size_t nVertices = meshData.nVertices;
			const size_t nFaces = meshData.nFaces;
			const auto& fNormals = meshData.fNormals;
			const auto& fCenters = meshData.fCenters;
			const auto& vNormals = meshData.vNormals;

			// Some kind of intrinsic vector field
			// Project this weird swirly field on to the surface (the ABC flow)
			auto spatialFunc = [&](glm::vec3 p) {
				float A = 1.;
				float B = 1.;
				float C = 1.;
				float xComp = A * std::sin(p.z) + C * std::cos(p.y);
				float yComp = B * std::sin(p.x) + A * std::cos(p.z);
				float zComp = C * std::sin(p.y) + B * std::cos(p.x);
				return glm::vec3{ xComp, yComp, zComp };
				};

			auto constructBasis = [&](glm::vec3 unitNormal) -> std::tuple<glm::vec3, glm::vec3> {
				glm::vec3 basisX{ 1., 0., 0. };
				basisX -= dot(basisX, unitNormal) * unitNormal;
				if (std::abs(basisX.x) < 0.1) {
					basisX = glm::vec3{ 0., 1., 0. };
					basisX -= glm::dot(basisX, unitNormal) * unitNormal;
				}
				basisX = glm::normalize(basisX);
				glm::vec3 basisY = glm::normalize(glm::cross(unitNormal, basisX));
				return std::make_tuple(basisX, basisY);
				};

			// vertex tangent bases
			std::vector<glm::vec3> vertexBasisX(nVertices);
			std::vector<glm::vec3> vertexBasisY(nVertices);
			for (size_t i = 0; i < nVertices; i++) {
				std::tie(vertexBasisX[i], vertexBasisY[i]) = constructBasis(vNormals[i]);
			}

			// face tangent bases
			std::vector<glm::vec3> faceBasisX(nFaces);
			std::vector<glm::vec3> faceBasisY(nFaces);
			for (size_t i = 0; i < nFaces; i++) {
				std::tie(faceBasisX[i], faceBasisY[i]) = constructBasis(fNormals[i]);
			}

			// At vertices
			std::vector<glm::vec2> vertexTangentVec(nVertices, glm::vec3{ 0., 0., 0. });
			for (size_t iV = 0; iV < nVertices; iV++) {
				glm::vec3 pos = meshData.vertexPositionsGLM[iV];
				glm::vec3 basisX = vertexBasisX[iV];
				glm::vec3 basisY = vertexBasisY[iV];

				glm::vec3 v = spatialFunc(pos);
				glm::vec2 vTangent{ glm::dot(v, basisX), glm::dot(v, basisY) };
				vertexTangentVec[iV] = vTangent;
			}
			if (isReg) {
				meshData.vertexTangentVectorQuantities["tangent vertex vec"] = psMesh->addVertexTangentVectorQuantity(
					"tangent vertex vec", vertexTangentVec, vertexBasisX, vertexBasisY);
				meshData.vertexTangentVectorQuantities["tangent vertex vec line"] = psMesh->addVertexTangentVectorQuantity(
					"tangent vertex vec line", vertexTangentVec, vertexBasisX,
					vertexBasisY, 2);
			}
			else {
				meshData.vertexTangentVectorQuantities["tangent vertex vec"]->updateData(vertexTangentVec);
				meshData.vertexTangentVectorQuantities["tangent vertex vec"]->tangentBasisX.data = vertexBasisX;
				meshData.vertexTangentVectorQuantities["tangent vertex vec"]->tangentBasisY.data = vertexBasisY;

				meshData.vertexTangentVectorQuantities["tangent vertex vec line"]->updateData(vertexTangentVec);
				meshData.vertexTangentVectorQuantities["tangent vertex vec line"]->tangentBasisX.data = vertexBasisX;
				meshData.vertexTangentVectorQuantities["tangent vertex vec line"]->tangentBasisY.data = vertexBasisY;
			}


			// At faces
			std::vector<glm::vec2> faceTangentVec(nFaces, glm::vec3{ 0., 0., 0. });
			for (size_t iF = 0; iF < nFaces; iF++) {

				glm::vec3 pos = fCenters[iF];
				glm::vec3 basisX = faceBasisX[iF];
				glm::vec3 basisY = faceBasisY[iF];

				glm::vec3 v = spatialFunc(pos);
				glm::vec2 vTangent{ glm::dot(v, basisX), glm::dot(v, basisY) };
				faceTangentVec[iF] = vTangent;
			}
			if (isReg) {
				meshData.faceTangentVectorQuantities["tangent face vec"] = psMesh->addFaceTangentVectorQuantity(
					"tangent face vec", faceTangentVec, faceBasisX, faceBasisY);
				meshData.faceTangentVectorQuantities["tangent face vec cross"] = psMesh->addFaceTangentVectorQuantity(
					"tangent face vec cross", faceTangentVec, faceBasisX, faceBasisY,
					4);
			}
			else {
				meshData.faceTangentVectorQuantities["tangent face vec"]->updateData(vertexTangentVec);
				meshData.faceTangentVectorQuantities["tangent face vec"]->tangentBasisX.data = faceBasisX;
				meshData.faceTangentVectorQuantities["tangent face vec"]->tangentBasisY.data = faceBasisY;

				meshData.faceTangentVectorQuantities["tangent face vec cross"]->updateData(vertexTangentVec);
				meshData.faceTangentVectorQuantities["tangent face vec cross"]->tangentBasisX.data = faceBasisX;
				meshData.faceTangentVectorQuantities["tangent face vec cross"]->tangentBasisY.data = faceBasisY;
			}

			// 1-form
			const size_t nEdges = meshData.nEdges;
			auto& seenEdges = meshData.seenEdges;
			std::vector<float> edgeForm(nEdges, 0.);
			std::vector<char> edgeOrient(nEdges, false);
			bool isTriangle = true;
			size_t iEdge = 0;
			seenEdges.clear();
			for (size_t iF = 0; iF < nFaces; iF++) {
				const std::vector<size_t>& face = meshData.faceIndices[iF];

				if (face.size() != 3) {
					isTriangle = false;
					break;
				}

				glm::vec3 pos = fCenters[iF];

				for (size_t j = 0; j < face.size(); j++) {
					size_t vA = face[j];
					size_t vB = face[(j + 1) % face.size()];
					size_t iMin = std::min(vA, vB);
					size_t iMax = std::max(vA, vB);
					auto p = std::make_pair(iMin, iMax);
					if (seenEdges.find(p) ==
						seenEdges.end()) { // use the hashset again to iterate over edges in order
						glm::vec3 v = spatialFunc(pos);
						glm::vec3 edgeVec = meshData.vertexPositionsGLM[vB] - meshData.vertexPositionsGLM[vA];
						edgeForm[iEdge] = glm::dot(edgeVec, v);
						edgeOrient[iEdge] = (vB > vA);
						seenEdges.insert(p);
						iEdge++;
					}
				}
			}
			if (isTriangle) {
				if (isReg) {
					meshData.oneFormTangentVectorQuantities["intrinsic 1-form"] = psMesh->addOneFormTangentVectorQuantity(
						"intrinsic 1-form", edgeForm, edgeOrient);
				}
				else {
					meshData.oneFormTangentVectorQuantities["intrinsic 1-form"]->oneForm = edgeForm;
					meshData.oneFormTangentVectorQuantities["intrinsic 1-form"]->canonicalOrientation = edgeOrient;
				}
			}
		}

		void constructCurveNetwork(MeshData& meshData, pCurveNetwork* curveNetwork, bool isReg) {
			const auto& edges = meshData.edges;
			const auto& vertexPositionsGLM = meshData.vertexPositionsGLM;

			// Useful data
			size_t nNodes = vertexPositionsGLM.size();
			size_t nEdges = edges.size();

			{ // Add some node values
				std::vector<double> valX(nNodes);
				std::vector<double> valXabs(nNodes);
				std::vector<std::array<double, 3>> randColor(nNodes);
				std::vector<glm::vec3> randVec(nNodes);
				for (size_t iN = 0; iN < nNodes; iN++) {
					valX[iN] = vertexPositionsGLM[iN].x;
					valXabs[iN] = std::fabs(vertexPositionsGLM[iN].x);
					randColor[iN] = { {polyscope::randomUnit(), polyscope::randomUnit(), polyscope::randomUnit()} };
					randVec[iN] = glm::vec3{ polyscope::randomUnit() - .5, polyscope::randomUnit() - .5,
											polyscope::randomUnit() - .5 };
				}
				if (isReg) {
					meshData.curveNetworkNodeScalarQuantities["nX"] = curveNetwork->addNodeScalarQuantity("nX",
						valX);
					meshData.curveNetworkNodeScalarQuantities["nXabs"] = curveNetwork->addNodeScalarQuantity(
						"nXabs", valXabs);

					meshData.curveNetworkNodeColorQuantities["nColor"] = curveNetwork->addNodeColorQuantity(
						"nColor", randColor);

					meshData.curveNetworkNodeVectorQuantities["randVecN"] = curveNetwork->addNodeVectorQuantity(
						"randVecN", randVec);
				}
				else {
					meshData.curveNetworkNodeScalarQuantities["nX"]->updateData(valX);
					meshData.curveNetworkNodeScalarQuantities["nXabs"]->updateData(valXabs);

					meshData.curveNetworkNodeColorQuantities["nColor"]->updateData(randColor);

					meshData.curveNetworkNodeVectorQuantities["randVecN"]->updateData(randVec);
				}
			}

			{ // Add some edge values
				std::vector<double> edgeLen(nEdges);
				std::vector<std::array<double, 3>> randColor(nEdges);
				std::vector<glm::vec3> randVec(nEdges);
				for (size_t iE = 0; iE < nEdges; iE++) {
					auto edge = edges[iE];
					size_t nA = std::get<0>(edge);
					size_t nB = std::get<1>(edge);

					edgeLen[iE] = glm::length(vertexPositionsGLM[nA] - vertexPositionsGLM[nB]);
					randColor[iE] = { {polyscope::randomUnit(), polyscope::randomUnit(), polyscope::randomUnit()} };
					randVec[iE] = glm::vec3{ polyscope::randomUnit() - .5, polyscope::randomUnit() - .5,
											polyscope::randomUnit() - .5 };
				}
				if (isReg) {
					meshData.curveNetworkEdgeScalarQuantities["edge len"] = curveNetwork->addEdgeScalarQuantity(
						"edge len", edgeLen, polyscope::DataType::MAGNITUDE);

					meshData.curveNetworkEdgeColorQuantities["eColor"] = curveNetwork->addEdgeColorQuantity(
						"eColor", randColor);

					meshData.curveNetworkEdgeVectorQuantities["randVecE"] = curveNetwork->addEdgeVectorQuantity(
						"randVecE", randVec);
				}
				else {
					meshData.curveNetworkEdgeScalarQuantities["edge len"]->updateData(edgeLen);

					meshData.curveNetworkEdgeColorQuantities["eColor"]->updateData(randColor);

					meshData.curveNetworkEdgeVectorQuantities["randVecE"]->updateData(randVec);
				}

			}

			// set a node radius quantity from above
			curveNetwork->setNodeRadiusQuantity("nXabs");
		}

		void setTriMeshInternalData(bool isReg = false) {
			// triangle mesh
			if (triMeshData.meshName.empty()) return;
			auto& vertexPositionsGLM = triMeshData.vertexPositionsGLM;
			auto& faceIndices = triMeshData.faceIndices;

			const size_t nVertices = ourTriMesh->numVertices();
			triMeshData.nVertices = nVertices;
			vertexPositionsGLM.resize(triMeshData.nVertices);
			const size_t nFaces = ourTriMesh->numFaces();
			triMeshData.nFaces = nFaces;
			faceIndices.resize(triMeshData.nFaces);

			const auto& vertMat = ourTriMesh->getVertMat();
#pragma omp parallel for
			for (int i = 0; i < nVertices; ++i) {
				const auto& v = vertMat.row(i);
				vertexPositionsGLM[i] = glm::vec3{
						static_cast<float>(v(0)),
						static_cast<float>(v(1)),
						static_cast<float>(v(2)),
				};
			}

			const auto& faceMat = ourTriMesh->getFaceMat();
#pragma omp parallel for
			for (int i = 0; i < triMeshData.nFaces; ++i) {
				const auto& f = faceMat.row(i);
				faceIndices[i].emplace_back(f(0));
				faceIndices[i].emplace_back(f(1));
				faceIndices[i].emplace_back(f(2));
			}

			if (isReg) {
				psTriMesh = polyscope::registerSurfaceMesh(triMeshData.meshName, vertexPositionsGLM,
					faceIndices);
			}

			/* add some quantities */
			// scalar
			internal::setVertexScalarQuantity(triMeshData, psTriMesh, isReg);
			internal::setVertexScalarQuantity(triMeshData, psTriMesh, isReg);
			internal::setFaceScalarQuantity(triMeshData, psTriMesh, isReg);
			internal::setEdgeScalarQuantity(triMeshData, psTriMesh, isReg);

			// vector
			// Face & vertex normals
			auto& fNormals = triMeshData.fNormals;
			auto& fCenters = triMeshData.fCenters;
			auto& vNormals = triMeshData.vNormals;
			fNormals.clear(), fNormals.resize(nFaces);
			fCenters.clear(), fCenters.resize(nFaces);
			vNormals.clear(), vNormals.resize(nVertices, glm::vec3{ 0., 0., 0. });
			const auto& faceNormalMat = ourTriMesh->getFaceNormalMat();
			const auto& vertNormalMat = ourTriMesh->getVertNormalMat();
#pragma omp parallel for
			for (int iF = 0; iF < nFaces; ++iF) {
				const std::vector<size_t>& face = faceIndices[iF];

				// Compute a center (used follow)
				glm::vec3 C = { 0., 0., 0. };
				for (size_t iV : face) {
					C += vertexPositionsGLM[iV];
				}
				C /= face.size();
				fCenters[iF] = C;

				fNormals[iF] = glm::vec3{ (float)faceNormalMat(iF, 0),
										 (float)faceNormalMat(iF, 1),
										 (float)faceNormalMat(iF, 2) };
			}
#pragma omp parallel for
			for (int iV = 0; iV < nVertices; ++iV) {
				vNormals[iV] = glm::vec3{ (float)vertNormalMat(iV, 0),
										 (float)vertNormalMat(iV, 1),
										 (float)vertNormalMat(iV, 2) };
			}
			internal::setNormalQuantity(triMeshData, psTriMesh, isReg);
			internal::setIntrisicVectorQuantity(triMeshData, psTriMesh, isReg);

			// Add a curve network from the edges
			triMeshData.edges.clear();
			for (size_t iF = 0; iF < nFaces; iF++) {
				std::vector<size_t>& face = faceIndices[iF];
				for (size_t iV = 0; iV < face.size(); iV++) {
					size_t i0 = face[iV];
					size_t i1 = face[(iV + 1) % face.size()];
					if (i0 < i1) {
						triMeshData.edges.push_back({ i0, i1 });
					}
				}
			}
			// Add the curve
			if (!triMeshData.edges.empty()) {
				if (isReg) {
					triMeshData.curveName = triMeshData.meshName + " curves";
					triMeshCurveNetwork = polyscope::registerCurveNetwork(triMeshData.curveName,
						vertexPositionsGLM,
						triMeshData.edges);
				}
				internal::constructCurveNetwork(triMeshData, triMeshCurveNetwork, isReg);
			}
		}

		void setSurfaceMeshInternalData(bool isReg = false) {
			// surface mesh
			if (surfaceMeshData.meshName.empty()) return;
			auto& vertexPositionsGLM = surfaceMeshData.vertexPositionsGLM;
			auto& faceIndices = surfaceMeshData.faceIndices;

			const size_t nVertices = ourSurfaceMesh->numVertices();
			surfaceMeshData.nVertices = nVertices;
			if (nVertices == 0) {
				LOG::WARN("[PCO][Viewer] The number of surface mesh's vertices is zero, which will not be rendered.");
				return;
			}
			vertexPositionsGLM.resize(nVertices);
			const size_t nFaces = ourSurfaceMesh->numFaces();
			surfaceMeshData.nFaces = nFaces;
			faceIndices.resize(nFaces);

			const auto& vertVec = ourSurfaceMesh->getVertVec();
#pragma omp parallel for
			for (int i = 0; i < nVertices; ++i) {
				const auto& v = vertVec.at(i);
				vertexPositionsGLM[i] = glm::vec3{
						static_cast<float>(v(0)),
						static_cast<float>(v(1)),
						static_cast<float>(v(2)),
				};
			}

			const auto& faceVec = ourSurfaceMesh->getFaceVec();
#pragma omp parallel for
			for (int i = 0; i < nFaces; ++i) {
				const auto& f = faceVec.at(i);
				faceIndices[i] = f;
			}

			if (isReg) {
				psSurfaceMesh = polyscope::registerSurfaceMesh(surfaceMeshData.meshName, vertexPositionsGLM,
					faceIndices);
			}

			/* add some quantities */
			// scalar
			internal::setVertexScalarQuantity(surfaceMeshData, psSurfaceMesh, isReg);
			internal::setFaceScalarQuantity(surfaceMeshData, psSurfaceMesh, isReg);

			// vector
			// Face & vertex normals
			auto& fNormals = surfaceMeshData.fNormals;
			auto& fCenters = surfaceMeshData.fCenters;
			auto& vNormals = surfaceMeshData.vNormals;
			fNormals.clear(), fNormals.resize(nFaces);
			fCenters.clear(), fCenters.resize(nFaces);
			vNormals.clear(), vNormals.resize(nVertices, glm::vec3{ 0., 0., 0. });
			const auto& faceNormalMat = ourSurfaceMesh->getFaceNormalMat();
			const auto& vertNormalMat = ourSurfaceMesh->getVertNormalMat();
#pragma omp parallel for
			for (int iF = 0; iF < nFaces; ++iF) {
				const std::vector<size_t>& face = faceIndices[iF];

				// Compute a center (used follow)
				glm::vec3 C = { 0., 0., 0. };
				for (size_t iV : face) {
					C += vertexPositionsGLM[iV];
				}
				C /= face.size();
				fCenters[iF] = C;

				fNormals[iF] = glm::vec3{ (float)faceNormalMat(iF, 0),
										 (float)faceNormalMat(iF, 1),
										 (float)faceNormalMat(iF, 2) };
			}
#pragma omp parallel for
			for (int iV = 0; iV < nVertices; ++iV) {
				vNormals[iV] = glm::vec3{ (float)vertNormalMat(iV, 0),
										 (float)vertNormalMat(iV, 1),
										 (float)vertNormalMat(iV, 2) };
			}
			internal::setNormalQuantity(surfaceMeshData, psSurfaceMesh, isReg);

			// Add a curve network from the edges
			surfaceMeshData.edges.clear();
			for (size_t iF = 0; iF < nFaces; iF++) {
				std::vector<size_t>& face = faceIndices[iF];
				for (size_t iV = 0; iV < face.size(); iV++) {
					size_t i0 = face[iV];
					size_t i1 = face[(iV + 1) % face.size()];
					if (i0 < i1) {
						surfaceMeshData.edges.push_back({ i0, i1 });
					}
				}
			}

			// Add the curve
			if (!surfaceMeshData.edges.empty()) {
				if (isReg) {
					surfaceMeshData.curveName = surfaceMeshData.meshName + " curves";
					surfaceMeshCurveNetwork = polyscope::registerCurveNetwork(surfaceMeshData.curveName,
						vertexPositionsGLM,
						surfaceMeshData.edges);

				}
				internal::constructCurveNetwork(surfaceMeshData, surfaceMeshCurveNetwork, isReg);
			}
		}

		/* callback */
		void callback() {
			ImGui::PushItemWidth(100);

			/// TODO
			/*if (ImGui::CollapsingHeader("Mesh", ImGuiTreeNodeFlags_DefaultOpen)) {
				float w = ImGui::GetContentRegionAvail().x;
				float p = ImGui::GetStyle().FramePadding.x;
				if (ImGui::Button("Load##Mesh", ImVec2((w - p) / 2.f, 0))) {
					fileDialogOpen();
				}
				ImGui::SameLine(0, p);
				if (ImGui::Button("Save##Mesh", ImVec2((w - p) / 2.f, 0))) {
					fileDialogSave();
				}
			}

			if (psTriMesh != nullptr) {
				if (ImGui::CollapsingHeader("Triangle Mesh manipulation", ImGuiTreeNodeFlags_DefaultOpen)) {
					float w = ImGui::GetContentRegionAvail().x;
					float p = ImGui::GetStyle().FramePadding.x;

					static float noisePer = 0.5;
					static float minVal = 0;
					static float maxVal;
					ImGui::InputFloat("noise percentage", &noisePer);
					ImGui::InputFloat("min value", &minVal);
					ImGui::InputFloat("max value", &maxVal);
					ImGui::SameLine(0, p);
					if (ImGui::Button("Add noise", ImVec2((w - p) / 2.f, 0))) {
						ourTriMesh->addNoise(noisePer, minVal, maxVal);
						setTriMeshInternalData();
					}
				}
			}

			if (psSurfaceMesh != nullptr) {
				if (ImGui::CollapsingHeader("Surface Mesh manipulation", ImGuiTreeNodeFlags_DefaultOpen)) {
					float w = ImGui::GetContentRegionAvail().x;
					float p = ImGui::GetStyle().FramePadding.x;
					if (ImGui::Button("Stitch borders", ImVec2((w - p) / 2.f, 0))) {
						if (psSurfaceMesh->isEnabled()) {
							ourSurfaceMesh->stitch();
							setSurfaceMeshInternalData();
						}
					}
					ImGui::SameLine(0, p);
					if (ImGui::Button("Stitch and orient", ImVec2((w - p) / 2.f, 0))) {
						if (psSurfaceMesh->isEnabled()) {
							ourSurfaceMesh->stitchAndOrientPolygonSoups();
							setSurfaceMeshInternalData();
						}
					}

					if (ImGui::Button("Repair", ImVec2((w - p) / 2.f, 0))) {
						if (psSurfaceMesh->isEnabled()) {
							ourSurfaceMesh->repair();
							setSurfaceMeshInternalData();
						}
					}
					ImGui::SameLine(0, p);
					if (ImGui::Button("Reverse orientation", ImVec2((w - p) / 2.f, 0))) {
						if (psSurfaceMesh->isEnabled()) {
							ourSurfaceMesh->reverseOrientation();
							setSurfaceMeshInternalData();
						}
					}
				}
			}*/

			// ImGui::SameLine();
		}
	}

	MeshViewer::MeshViewer(const std::string& _meshName, SurfaceMesh* _surfaceMesh) {
		// Options
		// polyscope::options::autocenterStructures = true;
		// polyscope::view::windowWidth = 600;
		// polyscope::view::windowHeight = 800;
		// polyscope::options::maxFPS = -1;
		polyscope::options::verbosity = 100;
		polyscope::options::enableRenderErrorChecks = true;

		// Initialize polyscope
		polyscope::init();

		internal::ourSurfaceMesh = _surfaceMesh;
		internal::surfaceMeshData.meshName = _meshName;

		internal::setSurfaceMeshInternalData(true);
		render();
	}

	MeshViewer::MeshViewer(const std::string& _triMeshName, TriMesh* _triMesh,
		const std::string& _surfaceMeshName, SurfaceMesh* _surfaceMesh) {
		if (_triMeshName == _surfaceMeshName) {
			LOG::ERROR("[Viewer] Duplicate mesh name! Polyscope terminated.");
		}
		else {
			// Options
			// polyscope::options::autocenterStructures = true;
			// polyscope::view::windowWidth = 600;
			// polyscope::view::windowHeight = 800;
			// polyscope::options::maxFPS = -1;
			polyscope::options::verbosity = 100;
			polyscope::options::enableRenderErrorChecks = true;

			// Initialize polyscope
			polyscope::init();

			internal::ourTriMesh = _triMesh;
			internal::ourSurfaceMesh = _surfaceMesh;
			internal::triMeshData.meshName = _triMeshName;
			internal::surfaceMeshData.meshName = _surfaceMeshName;

			internal::setTriMeshInternalData(true);
			internal::setSurfaceMeshInternalData(true);
			render();
		}
	}

	void MeshViewer::render() {

		// Add a few gui elements
		polyscope::state::userCallback = internal::callback;

		// Show the gui
		polyscope::show();
	}

} // namespace viewer
NAMESPACE_END(PCO)
