#include "Config.hpp"
#include "core/PCOModel.hpp"

#include <string>
#include <filesystem>

#include <CLI/CLI.hpp>

using namespace PCO;
using namespace core;
using namespace utils;
using std::string;

int main(int argc, char **argv) {
    std::cout << "/////////////////////////////////////////////////////////\n";
    std::cout << "//                                                     //\n";
    std::cout << "//     PCO: Precision-Controllable Offset Surfaces     //\n";
    std::cout << "//                with Sharp Features                  //\n";
    std::cout << "//                                                     //\n";
    std::cout << "/////////////////////////////////////////////////////////\n\n";

    ////////////////////////////////////////////////////////

    CLI::App app("PCO");

    Args args;
    app.add_option("-f,--file", args.inFile,
                   "Input model's name with extension.")->required();
    app.add_option("-F,--File", args.outFile,
                   "Output path of offsetting surface.");

    app.add_option("-o,--offset", args.offsetFactor,
                   "Specify the offset factor.")->required();

    app.add_option("-d,--depth", args.maxOctreeDepth,
                   "Specify the maximum octree depth.")->required()
            ->check(CLI::Validator(
                    std::function < std::string(std::string & ) > ([&args](std::string &val) -> std::string {
                        if (args.isMerge && std::stoi(val) < 0)
                            throw CLI::ValidationError(
                                    "The maximum octree depth must be provided and larger than zero");
                        return "";
                    }), "octree check"));

    app.add_flag("-m,--merge", args.isMerge);
    app.add_option("-c,--comp", args.COMP_ANGLE,
                   "Specify the compatible angel threshold in merging fields.")
            ->check(CLI::Validator(
                    std::function < std::string(std::string & ) > ([&args](std::string &val) -> std::string {
                        if (args.isMerge && (std::stod(val) < .0 || std::stod(val) > 180.0))
                            throw CLI::ValidationError(
                                    "If merging distance fields is enabled, the compatible angel threshold must be provided and larger than zero");
                        return "";
                    }), "merge check"));

    app.add_flag("--pmp", args.postProcessing);

    app.add_flag("--view", args.enableView);

    app.add_flag("--fast", args.fastCompute);

    try {
        app.parse(argc, argv);
    }
    catch (const CLI::ValidationError &e) {
        return app.exit(e);
    }
    catch (const CLI::ParseError &e) {
        return app.exit(e);
    }

    ////////////////////////////////////////////////////////

    OffsetModel model(args.inFile, args.offsetFactor, args.maxOctreeDepth);

    static constexpr bool enableOut = true; // for debug and visualization
    std::string outDir;
    const std::string &modelName = model.modelName;
    if (args.outFile.empty()) {
        outDir = file::concatFilePath(VIS_DIR, modelName,
                                      std::to_string(args.offsetFactor),
                                      std::to_string(args.maxOctreeDepth));
        if (enableOut) {
            if (args.isMerge)
                args.outFile = file::concatFilePath(outDir, "offset_" + std::to_string(args.COMP_ANGLE) + ".obj");
            else
                args.outFile = file::concatFilePath(outDir, "offset.obj");
        }
    } else {
        size_t lastSeparatorPos = args.outFile.find_last_of("/\\");
        outDir = args.outFile.substr(0, lastSeparatorPos);
    }

    model.launch(args);

    return 0;
}