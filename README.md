
# Precision-Controllable Offset Surfaces with Sharp Features (SIGGRAPH Asia 2024)

[![Project Page](https://img.shields.io/badge/Project-Page-blue?style=flat&logo=google-chrome&logoColor=white)](https://alan-leo-wong.github.io/SIGASIA24-PCO-ProjectPage/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

<img src="https://raw.githubusercontent.com/Alan-Leo-Wong/SIGASIA24-PCO-ProjectPage/main/src/assets/gallery.png" width="800" alt="PCO Github Repo Teaser Image">

This repository contains the official implementation of the paper **"PCO: Precision-Controllable Offset Surfaces with Sharp Features"** accepted to _ACM Transactions on Graphics (Proceedings of SIGGRAPH ASIA 2024)_. Our method generates high-fidelity offset surfaces with preserved sharp features while maintaining user-controllable precision.

## Installation

### Prerequisites
- C++17 compiler
- CMake 3.20+
- Eigen 3.4+
- OpenMP (recommended for parallel computation)
- [libigl](https://github.com/libigl/libigl)
- [CGAL](https://github.com/CGAL/cgal)
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
./PCO -f input.obj -F output.obj -o 0.1 -d 8 --pmp

# With field merging and angle threshold
./PCO -f input.obj -F output.obj -o -0.1 -d 8 -m -c 170.0 --pmp
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
This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

**Note:** The [quick-cliques](https://github.com/darrenstrash/quick-cliques) dependency is separately licensed under GPLv3. Users are responsible for complying with both licenses when using the combined work.