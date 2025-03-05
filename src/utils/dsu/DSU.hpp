//
// Created by Lei on 12/15/2023.
//

#ifndef LINEARSHARPOFFSET_DSU_HPP
#define LINEARSHARPOFFSET_DSU_HPP

#include "Config.hpp"
#include "utils/Log.hpp"
#include "detail/BasicDataType.hpp"

#include <iostream>
#include <utility>

NAMESPACE_BEGIN(PCO)
namespace utils::dsu {

	// TODO: REPLACED WITH "META TEMPLATE"
	class DSUBase {
	public:
		using DATA_TYPE = int;
		using FA_TYPE = std::vector<DATA_TYPE>;
		using MFA_TYPE = std::unordered_map<DATA_TYPE, DATA_TYPE>;
		using CHILDREN_TYPE = std::unordered_set<DATA_TYPE>;

	protected:
		virtual void init(DATA_TYPE n) = 0;

		// TODO: Polymorphisms at Compile time
		// template<typename DATA>
		// void init(const DATA &data);

//            virtual void init(std::vector<_DATA_TYPE> const &in_data) = 0;

		virtual DATA_TYPE find(DATA_TYPE x) const = 0;

		virtual DATA_TYPE merge(DATA_TYPE x, DATA_TYPE y) = 0;
	};

	/////////////////////////////////
	//!  DSU using normal vector
	/////////////////////////////////
	class DSU : public DSUBase {
	public:
		FA_TYPE fa;

	public:
		DSU() noexcept = default;

		explicit DSU(FA_TYPE _fa) noexcept : fa(std::move(_fa)) {}

		virtual ~DSU() noexcept = default;

	public:
		void init(DATA_TYPE n) override {
			fa.resize(n);
			std::iota(fa.begin(), fa.end(), 0);
		}

		DATA_TYPE find(DATA_TYPE x) const override {
			return fa[x] == x ? x : find(fa[x]);
		}

		DATA_TYPE merge(DATA_TYPE x, DATA_TYPE y) override {
			fa[find(y)] = find(x);
			return find(x);
		}
	};

	/////////////////////////////////
	//!  DSU using unordered_map
	/////////////////////////////////
	class MDSU : public DSUBase {
	public:
		MFA_TYPE fa;
		std::unordered_map<DATA_TYPE, CHILDREN_TYPE> children;

	public:
		MDSU() noexcept = default;

		explicit MDSU(MFA_TYPE _fa) noexcept : fa(std::move(_fa)) {}

		virtual ~MDSU() noexcept = default;

	public:
		void init(DATA_TYPE n) override {
			for (DATA_TYPE i = 0; i < n; ++i) {
				fa[i] = i;
			}
		}

		DATA_TYPE find(DATA_TYPE x) const override {
			if (fa.find(x) == fa.end()) {
				//std::cout << "error x " << x << std::endl;
				return -1;
			}
			return (fa.at(x) == x) ? x : find(fa.at(x));
		}

		DATA_TYPE merge(DATA_TYPE x, DATA_TYPE y) override {
			DATA_TYPE fa_x = find(x);
			DATA_TYPE fa_y = find(y);
			fa[fa_y] = fa_x;

			children[fa_x].insert(fa_y);
			children[fa_x].insert(children[fa_y].begin(), children[fa_y].end());

			return fa_x;
		}

		CHILDREN_TYPE getChildren(DATA_TYPE x) const {
			if (children.empty()) return MDSU::CHILDREN_TYPE();
			DATA_TYPE fa_x = find(x);
			if (!children.empty() && !children.contains(fa_x)) {
				LOG::ERROR("[DSU] Invalid query in mdsu!");
				return MDSU::CHILDREN_TYPE();
			}
			return children.at(fa_x);
		}
	};

} // namespace utils::dsu
NAMESPACE_END(PCO)

#endif //LINEARSHARPOFFSET_DSU_HPP
