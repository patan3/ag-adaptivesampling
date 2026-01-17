#ifndef CAMERA_H
#define CAMERA_H
//==============================================================================================
// Originally written in 2016 by Peter Shirley <ptrshrl@gmail.com>
//
// To the extent possible under law, the author(s) have dedicated all copyright and related and
// neighboring rights to this software to the public domain worldwide. This software is
// distributed without any warranty.
//
// You should have received a copy (see file COPYING.txt) of the CC0 Public Domain Dedication
// along with this software. If not, see <http://creativecommons.org/publicdomain/zero/1.0/>.
//==============================================================================================

#include "hittable.h"
#include "material.h"


class camera {
  public:
    double aspect_ratio      = 1.0;  // Ratio of image width over height
    int    image_width       = 100;  // Rendered image width in pixel count
    int    samples_per_pixel = 500;  // Count of max random samples for each pixel
    int    max_depth         = 10;   // Maximum number of ray bounces into scene
    color  background;               // Scene background color

    double vfov     = 90;              // Vertical view angle (field of view)
    point3 lookfrom = point3(0,0,0);   // Point camera is looking from
    point3 lookat   = point3(0,0,-1);  // Point camera is looking at
    vec3   vup      = vec3(0,1,0);     // Camera-relative "up" direction

    double defocus_angle = 0;  // Variation angle of rays through each pixel
    double focus_dist = 10;    // Distance from camera lookfrom point to plane of perfect focus

    double convergence_threshold = 0.01; // Threshold for adaptive sampling convergenge
    int    batch_size = 16;  // Every 16 samples, check convergence

    long long total_samples_shotted = 0;  // Counter for total rays fired
    long long converged_pixels = 0;         // Counter for total converged pixels

    // Helper for calculating luminance from RGB color
    static double get_luminance(const color& c) {
        return 0.2126 * c.x() + 0.7152 * c.y() + 0.0722 * c.z();
    }

    // Helper for checking if pixel has converged based on standard error
    static bool check_convergence(long N, double sum_L, double sum_L2, double threshold) {
        if (N < 2) return false;  // At least 2 samples for variance
        
        double mean = sum_L / N;
        double variance = (sum_L2 / N) - (mean * mean);
        
        if (variance < 0) variance = 0; // Ensure variance is positive .-.
        
        double std_error = std::sqrt(variance) / std::sqrt(N);
        return std_error < threshold;
    }

    // Helper for checking if pixel has converged based on contrast (luminance)
    static bool check_convergence_contrast(double min_L, double max_L, double threshold) {
        double sum = max_L + min_L;
        
        // Avoid division by zero
        if (std::abs(sum) < 1e-10) {
            return true;  // Converged (i.e., no significant luminance variation)
        }
        
        double contrast = (max_L - min_L) / sum;
        return contrast < threshold;
    }

    // Helper for checking if pixel has converged based on contrast metric (per-channel)
    static bool check_convergence_channels(double min_R, double max_R, 
                                          double min_G, double max_G,
                                          double min_B, double max_B,
                                          double threshold) {
        // Calculate contrast for each channel separately
        double contrast_R = 0.0, contrast_G = 0.0, contrast_B = 0.0;
        
        double sum_R = max_R + min_R;
        if (std::abs(sum_R) > 1e-10) {
            contrast_R = (max_R - min_R) / sum_R;
        }
        
        double sum_G = max_G + min_G;
        if (std::abs(sum_G) > 1e-10) {
            contrast_G = (max_G - min_G) / sum_G;
        }
        
        double sum_B = max_B + min_B;
        if (std::abs(sum_B) > 1e-10) {
            contrast_B = (max_B - min_B) / sum_B;
        }
        
        // Weight by human eye sensitivity (same weights as luminance)
        // Red: 0.2126, Green: 0.7152, Blue: 0.0722
        double weighted_contrast = 0.2126 * contrast_R + 0.7152 * contrast_G + 0.0722 * contrast_B;
        
        return weighted_contrast < threshold;
    }


    // Helper for checking convergence based on per-channel contrast (paper method)
    static bool check_convergence_channels_paper(double min_R, double max_R, 
                                                double min_G, double max_G,
                                                double min_B, double max_B) {
        // Separate thresholds per channel (from paper)
        const double threshold_R = 0.4;
        const double threshold_G = 0.3; 
        const double threshold_B = 0.6; 
        
        // Calculate contrast for each channel separately
        double contrast_R = 0.0, contrast_G = 0.0, contrast_B = 0.0;
        
        double sum_R = max_R + min_R;
        if (std::abs(sum_R) > 1e-10) {
            contrast_R = (max_R - min_R) / sum_R;
        }
        
        double sum_G = max_G + min_G;
        if (std::abs(sum_G) > 1e-10) {
            contrast_G = (max_G - min_G) / sum_G;
        }
        
        double sum_B = max_B + min_B;
        if (std::abs(sum_B) > 1e-10) {
            contrast_B = (max_B - min_B) / sum_B;
        }
        
        // Converged if ALL channels are below their respective thresholds
        return (contrast_R < threshold_R) && 
            (contrast_G < threshold_G) && 
            (contrast_B < threshold_B);
    }

    void render(const hittable& world) {
        initialize();

        std::cout << "P3\n" << image_width << ' ' << image_height << "\n255\n";

        // Sample distribution counters
        int count_low = 0;
        int count_mid = 0;
        int count_high = 0;

        for (int j = 0; j < image_height; j++) {
            std::clog << "\rScanlines remaining: " << (image_height - j) << ' ' << std::flush;
            for (int i = 0; i < image_width; i++) {
                color pixel_color(0,0,0);
                
                // Method: 0 = Variance, 1 = Luminance Contrast, 2 = Per-Channel Contrast, 3 = Per-Channel Mitchel Contrast
                int convergence_method = 0;
                
                // Track statistics for variance method
                double sum_L = 0.0;   // Sum of luminance values
                double sum_L2 = 0.0;  // Sum of squared luminance values
                
                // Track statistics for luminance contrast method
                double min_L = infinity;  // Minimum luminance
                double max_L = -infinity; // Maximum luminance
                
                // Track statistics for per-channel contrast method
                double min_R = infinity, max_R = -infinity;
                double min_G = infinity, max_G = -infinity;
                double min_B = infinity, max_B = -infinity;
                
                long current_samples = 0;
                
                // Continue until convergence OR max samples
                while (current_samples < samples_per_pixel) {
                    ray r = get_ray(i, j);
                    color sample_color = ray_color(r, max_depth, world);
                    pixel_color += sample_color;
                    
                    // Track statistics for all methods
                    double L = get_luminance(sample_color);
                    sum_L += L;
                    sum_L2 += L * L;
                    
                    // Luminance min/max for simplified contrast
                    if (L < min_L) min_L = L;
                    if (L > max_L) max_L = L;
                    
                    // Per-channel min/max for per-channel contrast
                    double R = sample_color.x();
                    double G = sample_color.y();
                    double B = sample_color.z();
                    if (R < min_R) min_R = R;
                    if (R > max_R) max_R = R;
                    if (G < min_G) min_G = G;
                    if (G > max_G) max_G = G;
                    if (B < min_B) min_B = B;
                    if (B > max_B) max_B = B;
                    
                    current_samples++;
                    
                    // Every batch_size samples, check convergence
                    if (convergence_threshold < 1e9 && current_samples % batch_size == 0) {
                        bool converged = false;
                        
                        if (convergence_method == 0) {
                            // Variance method 
                            converged = check_convergence(current_samples, sum_L, sum_L2, convergence_threshold);
                        } else if (convergence_method == 1) {
                            // Luminance contrast 
                            converged = check_convergence_contrast(min_L, max_L, convergence_threshold);
                        } else if (convergence_method == 2) {
                            // Per-channel contrast  
                            converged = check_convergence_channels(min_R, max_R, min_G, max_G, min_B, max_B, convergence_threshold);
                        } else if (convergence_method == 3) {
                            // Per-channel contrast (paper method)
                            converged = check_convergence_channels_paper(min_R, max_R, min_G, max_G, min_B, max_B);
                        }
                        
                        if (converged) {
                            converged_pixels++;
                            break;  // Pixel has converged :D
                        }
                    }
                }
                
                // Track total samples shot
                total_samples_shotted += current_samples;
                
                // Classify pixel difficulty
                if (current_samples < samples_per_pixel / 3) {
                    count_low++;
                } else if (current_samples > 2 * samples_per_pixel / 3) {
                    count_high++;
                } else {
                    count_mid++;
                }
                
                // Scale by actual number of samples taken
                double scale = 1.0 / current_samples;

                bool debug_view = false; // TRUE for HEATMAP, FALSE FOR NORMAL IMAGE

                if (debug_view) {
                    // Green = Low Samples, Red = High Samples
                    double t = double(current_samples - batch_size) / (samples_per_pixel - batch_size);
                    if (t < 0) t = 0;
                    if (t > 1) t = 1;
                    
                    // Simple interpolation from Green (0,1,0) to Red (1,0,0)
                    color heatmap(t, 1.0 - t, 0.0); 
                    
                    // Write heatmap color (already in [0,1] range)
                    write_color(std::cout, heatmap);
                } else {
                    // Normal Render Logic
                    write_color(std::cout, scale * pixel_color);
                }
            }
        }

        int total_pixels = image_width * image_height;

        std::clog << "\rDone.                 \n";
        std::clog << "Total Rays Fired: " << total_samples_shotted << "\n";
        std::clog << "Average Samples Per Pixel: " 
                  << (double)total_samples_shotted / (total_pixels) << "\n";
        std::clog << "Converged Pixels: " << converged_pixels << " / "
                  << (total_pixels) << "  " << converged_pixels * 100.0 / (total_pixels) << "%\n";
        
        // Sample distribution analysis
        std::clog << "\nSample Distribution:\n";
        std::clog << "  Low Effort Pixels  (< " << samples_per_pixel / 3 << " samples): " 
                  << count_low << " (" << (count_low * 100.0 / total_pixels) << "%)\n";
        std::clog << "  Medium Effort Pixels: " 
                  << count_mid << " (" << (count_mid * 100.0 / total_pixels) << "%)\n";
        std::clog << "  High Effort Pixels (> " << 2 * samples_per_pixel / 3 << " samples): " 
                  << count_high << " (" << (count_high * 100.0 / total_pixels) << "%)\n";
    }

  private:
    int    image_height;         // Rendered image height
    double pixel_samples_scale;  // Color scale factor for a sum of pixel samples
    point3 center;               // Camera center
    point3 pixel00_loc;          // Location of pixel 0, 0
    vec3   pixel_delta_u;        // Offset to pixel to the right
    vec3   pixel_delta_v;        // Offset to pixel below
    vec3   u, v, w;              // Camera frame basis vectors
    vec3   defocus_disk_u;       // Defocus disk horizontal radius
    vec3   defocus_disk_v;       // Defocus disk vertical radius

    void initialize() {
        image_height = int(image_width / aspect_ratio);
        image_height = (image_height < 1) ? 1 : image_height;

        pixel_samples_scale = 1.0 / samples_per_pixel;

        center = lookfrom;

        // Determine viewport dimensions.
        auto theta = degrees_to_radians(vfov);
        auto h = std::tan(theta/2);
        auto viewport_height = 2 * h * focus_dist;
        auto viewport_width = viewport_height * (double(image_width)/image_height);

        // Calculate the u,v,w unit basis vectors for the camera coordinate frame.
        w = unit_vector(lookfrom - lookat);
        u = unit_vector(cross(vup, w));
        v = cross(w, u);

        // Calculate the vectors across the horizontal and down the vertical viewport edges.
        vec3 viewport_u = viewport_width * u;    // Vector across viewport horizontal edge
        vec3 viewport_v = viewport_height * -v;  // Vector down viewport vertical edge

        // Calculate the horizontal and vertical delta vectors from pixel to pixel.
        pixel_delta_u = viewport_u / image_width;
        pixel_delta_v = viewport_v / image_height;

        // Calculate the location of the upper left pixel.
        auto viewport_upper_left = center - (focus_dist * w) - viewport_u/2 - viewport_v/2;
        pixel00_loc = viewport_upper_left + 0.5 * (pixel_delta_u + pixel_delta_v);

        // Calculate the camera defocus disk basis vectors.
        auto defocus_radius = focus_dist * std::tan(degrees_to_radians(defocus_angle / 2));
        defocus_disk_u = u * defocus_radius;
        defocus_disk_v = v * defocus_radius;
    }

    ray get_ray(int i, int j) const {
        // Construct a camera ray originating from the defocus disk and directed at a randomly
        // sampled point around the pixel location i, j.

        auto offset = sample_square();
        auto pixel_sample = pixel00_loc
                          + ((i + offset.x()) * pixel_delta_u)
                          + ((j + offset.y()) * pixel_delta_v);

        auto ray_origin = (defocus_angle <= 0) ? center : defocus_disk_sample();
        auto ray_direction = pixel_sample - ray_origin;
        auto ray_time = random_double();

        return ray(ray_origin, ray_direction, ray_time);
    }

    vec3 sample_square() const {
        // Returns the vector to a random point in the [-.5,-.5]-[+.5,+.5] unit square.
        return vec3(random_double() - 0.5, random_double() - 0.5, 0);
    }

    vec3 sample_disk(double radius) const {
        // Returns a random point in the unit (radius 0.5) disk centered at the origin.
        return radius * random_in_unit_disk();
    }

    point3 defocus_disk_sample() const {
        // Returns a random point in the camera defocus disk.
        auto p = random_in_unit_disk();
        return center + (p[0] * defocus_disk_u) + (p[1] * defocus_disk_v);
    }

    color ray_color(const ray& r, int depth, const hittable& world) const {
        // If we've exceeded the ray bounce limit, no more light is gathered.
        if (depth <= 0)
            return color(0,0,0);

        hit_record rec;

        // If the ray hits nothing, return the background color.
        if (!world.hit(r, interval(0.001, infinity), rec))
            return background;

        ray scattered;
        color attenuation;
        color color_from_emission = rec.mat->emitted(rec.u, rec.v, rec.p);

        if (!rec.mat->scatter(r, rec, attenuation, scattered))
            return color_from_emission;

        color color_from_scatter = attenuation * ray_color(scattered, depth-1, world);

        return color_from_emission + color_from_scatter;
    }
};


#endif
