# TIleThePlain

Graphics / image-processing experiment for loading a tile image, transforming it in matrix form, and previewing the result.

The repository mixes C++, OpenCV, Eigen, OpenGL / GLFW, and small Python helpers to explore tile preprocessing and geometric transforms.

## Highlights

- Loads tile imagery with OpenCV
- Converts image data between `cv::Mat` and Eigen matrices
- Applies rotation / transformation operations in matrix form
- Uses OpenGL / GLFW / ImGui dependencies for interactive visualization
- Includes Python helpers for preprocessing and viewing assets

## Tech Stack

- C++
- OpenCV
- Eigen3
- OpenGL
- GLFW
- ImGui
- Python
- CMake

## Repository Layout

```text
src/
├── main.cpp            # Main transform / preview experiment
├── HelperFunctions.*   # Conversion and transform helpers
├── preProcess.py       # Asset preprocessing helper
└── view.py             # Small viewing helper

data/                   # Input images and generated assets
HW3_data.zip            # Coursework data bundle
```

## Building

### Prerequisites

- CMake
- OpenCV
- Eigen3
- OpenGL
- GLFW

### Build

```bash
cmake -S . -B build
cmake --build build
```

The build generates the `main` target defined in `CMakeLists.txt`.

## Notes

- This repository reads like a focused graphics / coursework experiment rather than a productized application.
- Some paths in the current source are hard-coded for the original development environment and may need local adjustment before running.
