#pragma once

#include "Config.hpp"
#include "utils/Log.hpp"
#include "utils/File.hpp"
#include "detail/BasicDataType.hpp"

#include <string>
#include <iomanip>
#include <fstream>
#include <random>
#include <omp.h>

NAMESPACE_BEGIN(PCO)

//#define CUBE_WRITE

namespace geometry {
	// forward declare
	template<typename T>
	struct AABox;
} // namespace geometry

namespace geometry::gvis {
	// Helper function to write single vertex to OBJ file
	static void write_vertex(std::ostream& output, const Vector3& v) {
		output << "v " << v.x() << " " << v.y() << " " << v.z() << std::endl;
	}

	// Helper function to write single vertex to OBJ file
	static void write_vertex(std::ostream& output, const Vector3& v, const Vector3& rgb) {
		output << "v " << v.x() << " " << v.y() << " " << v.z() << " " << rgb.x() << " " << rgb.y() << " "
			<< rgb.z() << std::endl;
	}

	// Helper function to write single vertex to OBJ file
	static void write_vertex_to_xyz(std::ostream& output, const Vector3& v) {
		output << v.x() << " " << v.y() << " " << v.z() << std::endl;
	}

	// Helper function to write face
	static void write_face(std::ostream& output, const Vector3i& f) {
		output << "f " << f.x() << " " << f.y() << " " << f.z() << std::endl;
	}

	// Helper function to write face
	static void write_face(std::ostream& output, const Eigen::Vector4i& f) {
		output << "f " << f.x() << " " << f.y() << " " << f.z() << " " << f.w() << std::endl;
	}

	// Helper function to write line
	static void write_line(std::ostream& output, const Eigen::Vector2i& l) {
		output << "l " << l.x() << " " << l.y() << std::endl;
	}

	// Helper function to write full cube (using relative vertex positions in the OBJ file - support for this should be widespread by now)
	PCO_INLINE void
		writeCube(const Vector3& nodeOrigin, const Vector3& unit, const size_t& vertBegIdx, std::ostream& output) {
		//	   2-------1
		//	  /|      /|
		//	 / |     / |
		//	7--|----8  |
		//	|  4----|--3
		//	| /     | /
		//	5-------6
		// Create vertices
		Vector3 v1 = nodeOrigin + Vector3(0, unit.y(), unit.z());
		Vector3 v2 = nodeOrigin + Vector3(0, 0, unit.z());
		Vector3 v3 = nodeOrigin + Vector3(0, unit.y(), 0);
		const Vector3& v4 = nodeOrigin;
		Vector3 v5 = nodeOrigin + Vector3(unit.x(), 0, 0);
		Vector3 v6 = nodeOrigin + Vector3(unit.x(), unit.y(), 0);
		Vector3 v7 = nodeOrigin + Vector3(unit.x(), 0, unit.z());
		Vector3 v8 = nodeOrigin + Vector3(unit.x(), unit.y(), unit.z());

		// write them in reverse order, so relative position is -i for v_i
		write_vertex(output, v1);
		write_vertex(output, v2);
		write_vertex(output, v3);
		write_vertex(output, v4);
		write_vertex(output, v5);
		write_vertex(output, v6);
		write_vertex(output, v7);
		write_vertex(output, v8);

		// create faces
#if defined(MESH_WRITE)
			// back
		write_face(output, Eigen::Vector3i(vertBegIdx + 1, vertBegIdx + 3, vertBegIdx + 4));
		write_face(output, Eigen::Vector3i(vertBegIdx + 1, vertBegIdx + 4, vertBegIdx + 2));
		// bottom
		write_face(output, Eigen::Vector3i(vertBegIdx + 4, vertBegIdx + 3, vertBegIdx + 6));
		write_face(output, Eigen::Vector3i(vertBegIdx + 4, vertBegIdx + 6, vertBegIdx + 5));
		// right
		write_face(output, Eigen::Vector3i(vertBegIdx + 3, vertBegIdx + 1, vertBegIdx + 8));
		write_face(output, Eigen::Vector3i(vertBegIdx + 3, vertBegIdx + 8, vertBegIdx + 6));
		// top
		write_face(output, Eigen::Vector3i(vertBegIdx + 1, vertBegIdx + 2, vertBegIdx + 7));
		write_face(output, Eigen::Vector3i(vertBegIdx + 1, vertBegIdx + 7, vertBegIdx + 8));
		// left
		write_face(output, Eigen::Vector3i(vertBegIdx + 2, vertBegIdx + 4, vertBegIdx + 5));
		write_face(output, Eigen::Vector3i(vertBegIdx + 2, vertBegIdx + 5, vertBegIdx + 7));
		// front
		write_face(output, Eigen::Vector3i(vertBegIdx + 5, vertBegIdx + 6, vertBegIdx + 8));
		write_face(output, Eigen::Vector3i(vertBegIdx + 5, vertBegIdx + 8, vertBegIdx + 7));
#elif defined(CUBE_WRITE)
        // back
		write_face(output, Eigen::Vector4i(vertBegIdx + 3, vertBegIdx + 4, vertBegIdx + 2, vertBegIdx + 1));
		// bottom
		write_face(output, Eigen::Vector4i(vertBegIdx + 6, vertBegIdx + 5, vertBegIdx + 4, vertBegIdx + 3));
		// right
		write_face(output, Eigen::Vector4i(vertBegIdx + 1, vertBegIdx + 8, vertBegIdx + 6, vertBegIdx + 3));
		// top
		write_face(output, Eigen::Vector4i(vertBegIdx + 1, vertBegIdx + 2, vertBegIdx + 7, vertBegIdx + 8));
		// left
		write_face(output, Eigen::Vector4i(vertBegIdx + 4, vertBegIdx + 5, vertBegIdx + 7, vertBegIdx + 2));
		// front
		write_face(output, Eigen::Vector4i(vertBegIdx + 8, vertBegIdx + 7, vertBegIdx + 5, vertBegIdx + 6));
#else
		write_line(output, Eigen::Vector2i(vertBegIdx + 1, vertBegIdx + 2));
		write_line(output, Eigen::Vector2i(vertBegIdx + 2, vertBegIdx + 7));
		write_line(output, Eigen::Vector2i(vertBegIdx + 7, vertBegIdx + 8));
		write_line(output, Eigen::Vector2i(vertBegIdx + 8, vertBegIdx + 1));

		write_line(output, Eigen::Vector2i(vertBegIdx + 3, vertBegIdx + 4));
		write_line(output, Eigen::Vector2i(vertBegIdx + 4, vertBegIdx + 5));
		write_line(output, Eigen::Vector2i(vertBegIdx + 5, vertBegIdx + 6));
		write_line(output, Eigen::Vector2i(vertBegIdx + 6, vertBegIdx + 3));

		write_line(output, Eigen::Vector2i(vertBegIdx + 3, vertBegIdx + 1));
		write_line(output, Eigen::Vector2i(vertBegIdx + 4, vertBegIdx + 2));
		write_line(output, Eigen::Vector2i(vertBegIdx + 5, vertBegIdx + 7));
		write_line(output, Eigen::Vector2i(vertBegIdx + 6, vertBegIdx + 8));
#endif

		//vertBegIdx += 8;
	}

	PCO_INLINE void
		writeRectangle(const Eigen::Vector2d& nodeOrigin, const Eigen::Vector2d& unit, const size_t& vertBegIdx,
			std::ostream& output) {
		//  2-------1
		//  |       |
		//  |       |
		//  |       |
		//  4-------3

		// Create vertices
		Eigen::Vector3d o(nodeOrigin.x(), nodeOrigin.y(), 0);
		Eigen::Vector3d v1 = o + Eigen::Vector3d(unit.x(), unit.y(), 0);
		Eigen::Vector3d v2 = o + Eigen::Vector3d(0, unit.y(), 0);
		Eigen::Vector3d v3 = o + Eigen::Vector3d(unit.x(), 0, 0);
		const Eigen::Vector3d& v4 = o;

		// write them in reverse order, so relative position is -i for v_i
		write_vertex(output, v1);
		write_vertex(output, v2);
		write_vertex(output, v3);
		write_vertex(output, v4);

		// create faces
#if defined(MESH_WRITE)
		write_face(output, Eigen::Vector3i(vertBegIdx + 1, vertBegIdx + 2, vertBegIdx + 3));
		write_face(output, Eigen::Vector3i(vertBegIdx + 3, vertBegIdx + 2, vertBegIdx + 4));
#elif defined(CUBE_WRITE)
		write_face(output, Eigen::Vector4i(vertBegIdx + 1, vertBegIdx + 2, vertBegIdx + 4, vertBegIdx + 3));
#else
		write_line(output, Eigen::Vector2i(vertBegIdx + 1, vertBegIdx + 2));
		write_line(output, Eigen::Vector2i(vertBegIdx + 2, vertBegIdx + 4));
		write_line(output, Eigen::Vector2i(vertBegIdx + 4, vertBegIdx + 3));
		write_line(output, Eigen::Vector2i(vertBegIdx + 3, vertBegIdx + 1));
#endif
	}

	//        PCO_INLINE void write_tri(const std::vector<Vector3> &triVerts, int &vertBegIdx, std::ofstream &output) {
	//            if (triVerts.size() != 3) {
	//                LOG::qpWarn("Invalid number of triangle vertex!");
	//                return;
	//            }
	//            if (vertBegIdx <= 0) vertBegIdx = 1;
	//            for (int i = 0; i < 3; ++i) {
	//                output << "v " << triVerts[i].transpose() << "\n";
	//            }
	//            output << "f " << vertBegIdx << " " << vertBegIdx + 1 << " " << vertBegIdx + 2 << "\n";
	//            vertBegIdx += 3;
	//        }

	PCO_INLINE void
		write_tri(const Vector3& triVert_0, const Vector3& triVert_1, const Vector3& triVert_2,
			int& vertBegIdx, std::ostream& output) {
		if (vertBegIdx <= 0) vertBegIdx = 1;

		output << "v " << triVert_0.transpose() << "\n";
		output << "v " << triVert_1.transpose() << "\n";
		output << "v " << triVert_2.transpose() << "\n";
		output << "f " << vertBegIdx << " " << vertBegIdx + 1 << " " << vertBegIdx + 2 << "\n";

		vertBegIdx += 3;
	}

	PCO_INLINE void writePointCloud(const std::vector<Vector3>& points, std::ostream& output) {
		for (const auto& point : points)
			write_vertex(output, point);
	}

	PCO_INLINE void writePointCloud_xyz(const std::vector<Vector3>& points, std::ostream& output) {
		for (const auto& point : points)
			write_vertex_to_xyz(output, point);
	}

	PCO_INLINE void
		writePointCloud(const std::vector<Vector3>& points, const std::vector<Vector3>& rgbs, std::ostream& output) {
		for (size_t i = 0; i < points.size(); ++i)
			write_vertex(output, points[i], rgbs[i]);
	}

	PCO_INLINE void writePointCloud(const MatrixX& points, std::ostream& output) {
		for (size_t i = 0; i < points.rows(); ++i)
			write_vertex(output, points.row(i));
	}

	PCO_INLINE void writePointCloud_xyz(const MatrixX& points, std::ostream& output) {
		for (size_t i = 0; i < points.rows(); ++i)
			write_vertex_to_xyz(output, points.row(i));
	}

	PCO_INLINE void writePointCloud(const MatrixX& points, const std::vector<Vector3>& rgbs, std::ostream& output) {
		for (size_t i = 0; i < points.rows(); ++i)
			write_vertex(output, points.row(i), rgbs[i]);
	}

	PCO_INLINE void writePointCloud(const Vector3& point, const Vector3& rgb, std::ostream& output) {
		write_vertex(output, point, rgb);
	}
} // namespace geometry::gvis

namespace geometry::pointgen {

	template<typename Real>
	struct GaussianSampleGen {
	private:
		unsigned int numSamples;

	public:
		explicit GaussianSampleGen(const unsigned int& _numSamples) : numSamples(_numSamples) {}

		std::vector<Real> run(const AABox<Real>& box) {
			static std::random_device rd;
			static std::mt19937 gen(rd());
			std::normal_distribution<double> dist(0.0, 1.0);

			const int dim = box.boxOrigin.rows();
			Real mean = (box.boxEnd + box.boxOrigin) / 2.0;
			Real stddev = (box.boxEnd - box.boxOrigin) / 6.0;

			std::vector<Real> res;
			res.resize(numSamples);

			// std::ofstream out(".\\gaussian_points.xyz");

#pragma omp parallel for
			for (int i = 0; i < numSamples; ++i) {
				// 生成高斯样本
				Real sample;
				for (int j = 0; j < dim; ++j) {
					sample(j) = mean(j) + stddev(j) * dist(gen);
				}
				res[i] = sample;
				// 打印样本
				// std::cout << "Sample " << i + 1 << ": (" << sample(0) << ", " << sample(1) << ", " << sample(2) << ")" << std::endl;
				// out << sample(0) << " " << sample(1) << " " << sample(2) << std::endl;
			}
			return res;
		}
	};

	template<typename Real>
	struct UniformSampleGen {
	private:
		unsigned int numSamples;

		constexpr static const int Base[3] = {
				2, // X轴上的基数
				3, // Y轴上的基数
				5 // Z轴上的基数
		};

		// 生成Halton序列的第index个值
		double haltonSequence(int index, int base) const {
			double result = 0.0;
			double f = 1.0 / base;
			int i = index;

			while (i > 0) {
				result += f * (i % base);
				i = std::floor(i / base);
				f /= base;
			}

			return result;
		}

		// 将Halton序列值映射到[min, max]范围内
		double mapToRange(double value, double min, double max) const {
			return min + value * (max - min);
		}

		// 在[minArea, maxArea]范围内进行蓝噪声采样
		std::vector<Real> generateUniformSamples(const Real& minArea,
			const Real& maxArea) const {
			std::vector<Real> samples(numSamples);

			const int dim = minArea.rows(); // TODO: 使用Real内部的Rows, 实现constexpr
#pragma omp parallel for
			for (int i = 0; i < numSamples; ++i) {
				Real val;
				for (int j = 0; j < dim; ++j) {
					val(j) = mapToRange(haltonSequence(i, Base[j]), minArea(j), maxArea(j));
				}
				samples[i] = val;
				//                    std::cout << "val = " << val.transpose() << std::endl;
			}

			return samples;
		}

	public:
		explicit UniformSampleGen(const unsigned int& _numSamples) : numSamples(_numSamples) {}

		// Generate uniform samples in the given 3D space using Sobol sequence
		std::vector<Real>
			run(const AABox<Real>& box) const {
			return generateUniformSamples(box.boxOrigin, box.boxEnd);
		}
	};

	template<typename Real>
	struct UniformSampleGen_Sobol {
	private:
		unsigned int numSamples;

		// Sobol direction numbers for up to 21201 dimensions
		constexpr static const unsigned int sobolDirectionNumbers[][50] = {
			// Direction numbers for 1 dimension
			{1u, 0u},
			// Direction numbers for 2 dimensions
			{1u, 1u, 0u},
			// Direction numbers for 3 dimensions
			{1u, 3u, 1u, 0u},
			// Direction numbers for 4 dimensions
			{1u, 7u, 5u, 1u, 0u},
			// ... direction numbers for higher dimensions
		};

		// Helper function to calculate the number of bits required to represent a value
		unsigned int countBits(unsigned int value) {
			unsigned int count = 0;
			while (value > 0) {
				value >>= 1;
				count++;
			}
			return count;
		}

		// Helper function to calculate the ith component of the Sobol sequence
		double sobolSample(unsigned int i, unsigned int n) {
			unsigned int bits = countBits(n);
			unsigned int directionCount =
				sizeof(sobolDirectionNumbers) / sizeof(sobolDirectionNumbers[0]);

			if (i >= directionCount) {
				std::cerr << "Sobol sequence not supported for " << i << " dimensions."
					<< std::endl;
				return 0.0f;
			}

			const unsigned int* directionNumbers = sobolDirectionNumbers[i];
			unsigned int directionCountBits = countBits(directionNumbers[0]);

			if (bits > directionCountBits) {
				std::cerr << "Sobol sequence not supported for " << bits << " bits."
					<< std::endl;
				return 0.0f;
			}

			unsigned int result = 0;
			for (unsigned int j = 1; j <= bits; ++j) {
				if ((n & (1u << (bits - j))) != 0) {
					result ^= directionNumbers[j];
				}
			}

			return static_cast<double>(result) / static_cast<double>(1u << bits);
		}

	public:
		explicit UniformSampleGen_Sobol(const unsigned int& _numSamples) : numSamples(_numSamples) {}

		// Generate uniform samples in the given 3D space using Sobol sequence
		std::vector<Real>
			run(const AABox<Real>& box) {
			const int dim = box.boxEnd.rows();
			unsigned int currentSample = 0;

			std::vector<Real> samples;
			samples.resize(numSamples);

#pragma omp parallel for
			for (int i = 0; i < numSamples; ++i) {
				Real sample;
				for (unsigned int k = 0; k < dim; ++k) {
					sample(k) = box.boxOrigin(k) + sobolSample(k, i);
					sample(k) *= box.boxEnd(k) - box.boxOrigin(k);
				}
				samples[i] = sample;
			}

			return samples;
		}
	};
}

namespace geometry {
	template<typename Real>
	struct Edge {
		using type = std::pair<Real, Real>;
	};

	// An Axis Aligned Box (AAB) of a certain Real - to be initialized with a boxOrigin and boxEnd
	template<typename T>
	struct AABox {
		//using T = typename Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;
		T boxOrigin;
		T boxEnd;
		T boxWidth;

		AABox() = default;

		AABox(const T& _boxOrigin, const T& _boxEnd) : boxOrigin(_boxOrigin),
			boxEnd(_boxEnd),
			boxWidth(_boxEnd - _boxOrigin) {}

		void scaleAndTranslate(const double& scale_factor = 1.0, const T& translation = T::Zero()) {
			if (scale_factor == 1.0 && translation == T::Zero()) return;
			if (scale_factor <= .0) {
				std::cerr << "Invalid scale factor!\n";
				return;
			}

			// Calculate centroid point
			const T boxCenter = center();

			// Zoom and translation based on the centroid point
			const T scaled_min_point = (boxOrigin - boxCenter) * scale_factor + boxCenter + translation;
			const T scaled_max_point = (boxEnd - boxCenter) * scale_factor + boxCenter + translation;

			// Update coordinate
			update(scaled_min_point, scaled_max_point);
		}

		AABox(const AABox<T>& _box) {
			boxOrigin = _box.boxOrigin;
			boxEnd = _box.boxEnd;
			boxWidth = _box.boxWidth;
		}

		AABox& operator=(const AABox<T>& _box) {
			boxOrigin = _box.boxOrigin;
			boxEnd = _box.boxEnd;
			boxWidth = _box.boxWidth;

			return *this;
		}

		double diagLength() const { return boxWidth.norm(); }

		T center() const { return (boxOrigin + boxEnd) / 2.0; }

		T corner(int i) const {
            if constexpr (T::RowsAtCompileTime == 3) {
                if (i == 0) return boxOrigin;
                if (i == 1) return boxOrigin + T(boxWidth.x(), 0, 0);
                if (i == 2) return boxOrigin + T(0, boxWidth.y(), 0);
                if (i == 3) return boxOrigin + T(boxWidth.x(), boxWidth.y(), 0);

                if (i == 4) return boxOrigin + T(0, 0, boxWidth.z());
                if (i == 5) return boxOrigin + T(boxWidth.x(), 0, boxWidth.z());
                if (i == 6) return boxOrigin + T(0, boxWidth.y(), boxWidth.z());
                if (i == 7) return boxOrigin + T(boxWidth.x(), boxWidth.y(), boxWidth.z());
            }
            else if constexpr (T::RowsAtCompileTime == 2) {
                if (i == 0) return boxOrigin;
                if (i == 1) return boxOrigin + T(boxWidth.x(), 0);
                if (i == 2) return boxOrigin + T(0, boxWidth.y());
                if (i == 3) return boxOrigin + T(boxWidth.x(), boxWidth.y());
            }
		}

		std::vector<T> corner() const {
			std::vector<T> corners;
			/*for (int i = 0; i < 4; ++i) {
				T p = boxOrigin;
				p.x() += boxWidth.x() * (i & 1);
				p.y() += boxWidth.y() * ((i >> 1) & 1);
				corners.emplace_back(p);
			}*/

            if constexpr (T::RowsAtCompileTime == 3) {
                corners.resize(8);

                corners[0] = (boxOrigin);
                corners[1] = (boxOrigin + T(boxWidth.x(), 0, 0));
                corners[2] = (boxOrigin + T(0, boxWidth.y(), 0));
                corners[3] = (boxOrigin + T(boxWidth.x(), boxWidth.y(), 0));
                for (int i = 0; i < 4; ++i) {
                    T p = corners[i];
                    p.z() += boxWidth.z();
                    corners[i + 4] = p;
                }
            }
            else if constexpr (T::RowsAtCompileTime == 2) {
                corners.resize(4);

                corners[0] = (boxOrigin);
                corners[1] = (boxOrigin + T(boxWidth.x(), 0));
                corners[2] = (boxOrigin + T(0, boxWidth.y()));
                corners[3] = (boxOrigin + T(boxWidth.x(), boxWidth.y()));
            }
			return corners;
		}

		bool isInBox(const T& queryPoint) const {
			int dim = boxOrigin.rows();
			for (int i = 0; i < dim; ++i)
				if (!(boxOrigin(i) <= queryPoint(i) && queryPoint(i) <= boxEnd(i))) return false;
			return true;
		}

		void output(const std::string& filename) const {
			utils::file::checkDir(filename);
			std::ofstream out(filename);
			if (!out) {
				LOG::ERROR("[PCO][Geometry][I/O]: File \"{}\" could not be opened!", filename);
				return;
			}
			LOG::INFO("[PCO][Geometry] Output bounding-box to \"{}\" ...", filename);

			gvis::writeCube(boxOrigin, boxWidth, 0, out);

			out.close();
		}

		void output(std::ofstream& out, int& vert_beg_idx) const {
			if (vert_beg_idx <= 0) vert_beg_idx = 0;
			if constexpr (T::RowsAtCompileTime == 3) {
				gvis::writeCube(boxOrigin, boxWidth, vert_beg_idx, out);
				vert_beg_idx += 8;
			}
			else if constexpr (T::RowsAtCompileTime == 2) {
				gvis::writeRectangle(boxOrigin, boxWidth, vert_beg_idx, out);
				vert_beg_idx += 4;
			}
		}

		friend std::ostream& operator<<(std::ostream& out, const AABox<T>& bbox) {
			//if constexpr (std::is_base_of_v<Eigen::MatrixBase<T>, T>)
			if constexpr (T::RowsAtCompileTime == 3)
				gvis::writeCube(bbox.boxOrigin, bbox.boxWidth, 0, out);
			else if constexpr (T::RowsAtCompileTime == 2)
				gvis::writeRectangle(bbox.boxOrigin, bbox.boxWidth, 0, out);
			return out;
		}

	private:
		void update(const T& _boxOrigin, const T& _boxEnd) {
			boxOrigin = _boxOrigin;
			boxEnd = _boxEnd;
			boxWidth = boxEnd - boxOrigin;
		}
	};

	template<typename T>
	struct Triangle {
		T p1, p2, p3;
		T normal;

		Triangle() noexcept = default;

		Triangle(const T& _p1, const T& _p2, const T& _p3) noexcept : p1(_p1), p2(_p2), p3(_p3) {
			normal = (p2 - p1).cross(p3 - p1).normalized();
		}
	};

	//        template<typename Real = Scalar, class T = Vector<Real, 3>>
	//        struct Tet_ {
	//            typedef T Point;
	//            typedef std::array<Point, 4> Tetrahedron;
	//            typedef std::array<int, 4> TetNeighbor;
	//        };
	template<typename REAL>
	using Tet_ = Matrix<REAL, 4, 3>;
	using Tet = Tet_<Scalar>;
    using TetCoord = Tet;

} // namespace geometry

NAMESPACE_END(PCO)