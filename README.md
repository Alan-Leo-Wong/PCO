
# Precision-Controllable Offset Surfaces with Sharp Features (SIGGRAPH Asia 2024)

[![Project Page](https://img.shields.io/badge/Project-Page-blue?style=flat&logo=google-chrome&logoColor=white)](https://alan-leo-wong.github.io/SIGASIA24-PCO-ProjectPage/)
[![License: GPL](https://img.shields.io/badge/License-GPLv3.0-yellow.svg)](https://opensource.org/licenses/gpl-3.0)

<div style="width:100%; margin:0 auto; text-align:center;">
  <img src="https://raw.githubusercontent.com/Alan-Leo-Wong/SIGASIA24-PCO-ProjectPage/main/src/assets/gallery.png" 
       style="max-width:100%; min-width:600px; height:auto; border:1px solid #eee;">
</div>
This repository contains the official implementation of the paper **"PCO: Precision-Controllable Offset Surfaces with Sharp Features"** accepted to *ACM Transactions on Graphics (Proceedings of SIGGRAPH ASIA 2024)*. Our method generates high-fidelity offset surfaces with preserved sharp features while maintaining user-controllable precision.

## Installation

### Prerequisites
- C++17 compiler
- CMake 3.20+
- OpenMP (recommended for parallel computation)
- [Eigen](https://gitlab.com/libeigen/eigen/-/releases/3.4.0) 3.4+
- [spdlog](https://github.com/gabime/spdlog)
- [libigl](https://github.com/libigl/libigl)
- [CGAL](https://github.com/CGAL/cgal)
- [implicit predicates](https://github.com/qnzhou/implicit_predicates)
- [quick-cliques](https://github.com/darrenstrash/quick-cliques) (GPLv3 licensed)
- [polyscope](https://github.com/nmwsharp/polyscope) (optional, enable with `ENABLE_VIEWER`)

### Build Instructions
```bash
git clone --recursive https://github.com/alan-leo-wong/PCO.git
cd PCO

mkdir build
cd build

# Basic configuration
cmake .. -DCMAKE_BUILD_TYPE=Release
# or
# Optional features
cmake .. -DENABLE_VIEWER=ON -DUSE_SDF=ON  # Enable viewer and SDF features

cmake --build . -j your_core_num
```

## Useage

### Full Parameter Specification
| Parameter     | Description                             | Validation & Notes                                                      |
|---------------|-----------------------------------------|-------------------------------------------------------------------------|
| `-f,--file`   | Input model (OBJ/PLY/OFF)               | **Required**                                                            |
| `-F,--File`   | Output path                             | **Required**                                                            |
| `-o,--offset` | Offset distance (% of bbox diagonal)    | **Required**<br>• Positive: outward offset<br>• Negative: inward offset |
| `-d,--depth`  | Maximum octree depth                    | **Required**<br>• Must be ≥1 (recommend 8)                              |
| `-m,--merge`  | Enable distance field merging           | **Required when merging**                                               |
| `-c,--comp`   | Compatible angle threshold              | **Required with `-m`**<br>• 0°-180°                                     |
| `--pmp`       | Use CGAL for basic mesh post-processing | **Optional**                                                            |
| `--view`      | Enable visualization                    | **Requires `-DENABLE_VIEWER=ON`**                                           |


### Basic Example
```bash
# Basic offset generation
./PCO -f input.obj -F output.obj -o 2 -d 8 --pmp

# With field merging and angle threshold
./PCO -f input.obj -F output.obj -o -2 -d 8 -m -c 170.0 --pmp
```

## Citation
```bibtex
@article{wang2024pco,
    author = {Wang, Lei and Wang, Xudong and Wang, Pengfei and Chen, Shuangmin and Xin, Shiqing and Guo, Jiong and Wang, Wenping and Tu, Chenghe},
    title = {PCO: Precision-Controllable Offset Surfaces with Sharp Features},
    journal = {ACM Trans. Graph.},
    volume = {43},
    number = {6},
    doi = {10.1145/3687920},
    year = {2024},
    publisher = {ACM}
}
```

## License
Currently, this project is licensed under **GPLv3** due to dependencies on GPL-licensed components. See [LICENSE](LICENSE) for full terms.


## 🐛 Known Issues
**Current Constraints:** The quick-cliques component may occasionally fail on *complex meshes* or with *small compatible angles*. 
We have tested the implementation, though rare edge cases may still exist. Reproducible bug reports via [GitHub Issues](https://github.com/Alan-Leo-Wong/PCO/issues) are greatly appreciated to improve robustness.