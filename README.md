# Adaptive Sampling Ray Tracer

## Overview

This project implements adaptive sampling techniques for Monte Carlo path tracing, building upon the "Ray Tracing in One Weekend" series. Instead of using uniform sampling (allocating the same number of samples to every pixel), this renderer intelligently distributes computational resources by detecting when pixels have converged to stable colors.

The implementation compares two adaptive sampling methods:
- **Variance-based sampling**: Uses standard error of the mean to detect convergence
- **Contrast-based sampling**: Uses per-channel contrast metric to detect convergence

### Key Results
On the "bouncing_spheres" test scene:
- **Variance method**: 59.3% reduction in ray count (40.7% of uniform baseline)
- **Contrast method**: 25.8% reduction in ray count (74.2% of uniform baseline)

## Building the Project

### Prerequisites
- CMake 3.1 or higher
- C++ compiler with C++11 support

### Build Instructions

On Windows:

```shell
$ cmake -B build
$ cmake --build build --config Release  # Create release binaries in `build\Release`
```

On Linux / macOS:

```shell
# Configure and build release binaries under `build/Release`
$ cmake -B build/Release -DCMAKE_BUILD_TYPE=Release
$ cmake --build build/Release
```

### Running The Programs

You can run the programs by executing the binaries placed in the build directory:

To render a scene and save the output:

```powershell
.\build\Release\theNextWeek.exe > output.ppm
```

    $ build\Release\inOneWeekend > image.ppm

The renderer outputs a PPM image file and logs performance metrics to stderr (console).
The generated PPM file can be viewed directly as a regular computer image, if your operating system
supports this image type. If your system doesn't handle PPM files, then you should be able to find
PPM file viewers online. We used PPM Image Preview extension for Visual Studio Code.


## Running the Project

### Basic Usage



### Selecting a Scene

Edit `src/TheNextWeek/main.cc` and change the switch value in the `main()` function (line ~408):

```cpp
switch(1) {  // Change this number to select a scene
    case 1:  bouncing_spheres();          break;
    case 2:  checkered_spheres();         break;
    // ... other scenes
}
```

For the results in the report, use `case 1` (bouncing_spheres scene).

## Configuring Adaptive Sampling

All adaptive sampling parameters are configured in `src/TheNextWeek/camera.h`.

### Key Parameters

```cpp
// In the camera class:
int    samples_per_pixel = 500;        // Maximum samples per pixel
double convergence_threshold = 0.01;   // Threshold for convergence detection
int    batch_size = 16;                // Check convergence every N samples
```

### Selecting Convergence Method

In `camera.h`, locate the `render()` function (around line 133) and change the `convergence_method` variable:

```cpp
// Method: 0 = Variance, 1 = Luminance Contrast, 2 = Per-Channel Contrast, 3 = Per-Channel Mitchel Contrast
int convergence_method = 0;  // Change this value
```

- **Method 0 (Variance)**: Standard error of the mean convergence
- **Method 3 (Per-Channel Mitchell)**: Per-channel contrast using Mitchell's thresholds (R=0.4, G=0.3, B=0.6)

### Reproducing Report Results

#### Uniform Sampling (Baseline)
```cpp
double convergence_threshold = 1e10;  // Set very high to disable early termination
int    samples_per_pixel = 500;
```

#### Variance-Based Adaptive Sampling
```cpp
int    convergence_method = 0;
double convergence_threshold = 0.01;
int    samples_per_pixel = 500;
int    batch_size = 16;
```

#### Contrast-Based Adaptive Sampling
```cpp
int    convergence_method = 3;  // Use Mitchell's per-channel thresholds
double convergence_threshold = 0.01;  // Not used for method 3
int    samples_per_pixel = 500;
int    batch_size = 16;
```

### Generating Heatmaps

To visualize sample distribution, enable debug mode in `camera.h` (around line 215):

```cpp
bool debug_view = true;  // TRUE for HEATMAP, FALSE for NORMAL IMAGE
```

Heatmap colors:
- **Green**: Few samples (easy regions that converged quickly)
- **Red**: Many samples (complex regions requiring full budget)

## Performance Metrics

After rendering, the console displays:
- Total rays fired
- Average samples per pixel
- Number and percentage of converged pixels
- Sample distribution (low/medium/high effort pixels)

Example output:
```
Done.
Total Rays Fired: 73335244
Average Samples Per Pixel: 203.709
Converged Pixels: 273629 / 360000  76.0081%

Sample Distribution:
  Low Effort Pixels  (< 166 samples): 205517 (57.0881%)
  Medium Effort Pixels: 40512 (11.2533%)
  High Effort Pixels (> 333 samples): 113971 (31.6586%)
```

## Scene Configuration

For the bouncing_spheres scene, camera settings are in `main.cc`:

```cpp
cam.aspect_ratio      = 16.0 / 9.0;
cam.image_width       = 800;
cam.samples_per_pixel = 500;
cam.max_depth         = 50;
cam.vfov     = 20;
cam.lookfrom = point3(13,2,3);
cam.lookat   = point3(0,0,0);
cam.defocus_angle = 0.6;
cam.focus_dist    = 10.0;
```

Adjust `image_width` to change resolution (height is computed from aspect ratio).

## AI Disclosure

AI tools were used to assist in the development of this project in the following ways:

- **Debugging**: AI assisted in identifying and resolving bugs during development
- **Code optimization**: General improvements to code, such as optimizing online variance calculation for better performance, and adding guards for division by zero in contrast metrics
- **Code cleanup**: Refactoring and formatting improvements for better readability
- **Documentation**: AI assistance in writing and structuring this README documentation

All core algorithmic decisions, implementation choices, and experimental design were made manually.

## References

- **[Kajiya1986]**: Kajiya, J. T. "The Rendering Equation." SIGGRAPH 1986.
- **[Mitchell1987]**: Mitchell, D. P. "Generating Antialiased Images at Low Sampling Densities." SIGGRAPH 1987.
- **[Sik2013]**: Šik, M., et al. "A Survey of Adaptive Sampling in Realistic Image Synthesis." 2013.
- **[Ray Tracing in One Weekend: The Next Week]**: Shirley, P. "Ray Tracing in One Weekend Book Series." https://raytracing.github.io/






[book1]:           books/RayTracingInOneWeekend.html
[book2]:           books/RayTracingTheNextWeek.html
[book3]:           books/RayTracingTheRestOfYourLife.html
[CONTRIBUTING]:    CONTRIBUTING.md
[cover1]:          images/cover/CoverRTW1-small.jpg
[cover2]:          images/cover/CoverRTW2-small.jpg
[cover3]:          images/cover/CoverRTW3-small.jpg
[discussions]:     https://github.com/RayTracing/raytracing.github.io/discussions/
[GitHub home]:     https://github.com/RayTracing/raytracing.github.io/
[ImageMagick]:     https://imagemagick.org/
[implementations]: https://github.com/RayTracing/raytracing.github.io/wiki/Implementations
[issues]:          https://github.com/RayTracing/raytracing.github.io/issues/
[PRINTING.md]:     PRINTING.md
[web1]:            https://raytracing.github.io/books/RayTracingInOneWeekend.html
[web2]:            https://raytracing.github.io/books/RayTracingTheNextWeek.html
[web3]:            https://raytracing.github.io/books/RayTracingTheRestOfYourLife.html
