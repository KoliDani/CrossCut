# CrossCut

**CrossCut** is a preliminary C++ polygon geometry library focused on precise curve intersection detection. It supports segment, arc, and Bézier curve intersections using a sweep line algorithm for efficient spatial processing. Designed as a foundational building block for geometry pipelines, polygon clipping, and computational geometry tooling.

> ⚠️ **Work in progress** — Boolean polygon operations are currently under development and not yet fully functional.

---

## Features

- **Segment intersection** — precise line segment crossing detection
- **Arc intersection** — circular arc intersection support
- **Bézier intersection** — cubic and higher-order Bézier curve intersection via subdivision
- **Sweep line algorithm** — efficient spatial event processing for handling multiple curve intersections
- **Composable geometry primitives** — `Point`, `Line`, `Segment`, `ArcSegment`, `Bezier`, `BSpline`, `Path`, `Polygon` types built around a shared polymorphic base
- **Boolean operations** *(in progress)* — polygon union, intersection, difference and xor

---

## Project Structure

```
CrossCut/
├── CrossCut.cpp              # Entry point / usage example
├── CMakeLists.txt
├── data_structures/          # Core geometry types (Point, Path, Line, etc.)
├── bezier_intersector/       # Bézier intersection logic and sweep line implementation
├── boolean/                  # Boolean polygon operations (WIP)
└── utils/                    # Shared math utilities
```

---

## Getting Started

### Prerequisites

- CMake 3.15+
- C++17 compatible compiler (Clang or GCC)
- macOS: recommended to use GCC via Homebrew for OpenMP support

```bash
brew install cmake gcc
```

### Build

```bash
git clone https://github.com/KoliDani/CrossCut.git
cd CrossCut
mkdir build && cd build
cmake ..
cmake --build .
```

### Run

```bash
./CrossCut
```

---

## Roadmap

- [x] Segment intersection
- [x] Arc intersection
- [x] Bézier intersection
- [x] Sweep line event processing
- [ ] Boolean union
- [ ] Boolean intersection
- [ ] Boolean difference
- [ ] Python bindings

---

## License
MIT — see [LICENSE](LICENSE) for details.