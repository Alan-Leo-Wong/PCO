//
// Created by Lei on 4/1/2024.
//

#ifndef PCO_ARGS_HPP
#define PCO_ARGS_HPP

#include "Config.hpp"

#include <string>

NAMESPACE_BEGIN(PCO)
namespace core {

	struct Args {
		std::string inFile;
		std::string outFile;

		int maxOctreeDepth = -1;

		double offsetFactor;

		bool isMerge = false;
		double COMP_ANGLE = -1.0;

		bool postProcessing = false;
		bool enableView = false;
	};

} // namespace core
NAMESPACE_END(PCO)

#endif //PCO_ARGS_HPP
