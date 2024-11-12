#ifndef GENERATOR_H
#define GENERATOR_H

#include <algorithm> // For std::sort
#include <cmath>     // For atan2 and hypot
#include <vector>
#include <random>
#include <iostream>
#include <glm/glm.hpp>
#include <array>     // For std::array
#include <unordered_map> // Include the unordered_map header

#include <delaunator/delaunator.hpp>

#include <mapgen/data.h>
#include <mapgen/drawer.h>
#include <mapgen/config.h>

#include <mapgen/pass_water.h>
#include <mapgen/pass_biomes.h>
#include <mapgen/pass_coast.h>
#include <mapgen/pass_elevation.h>
#include <mapgen/pass_moisture.h>
#include <mapgen/pass_noisy_edges.h>
#include <mapgen/pass_rivers.h>
#include <mapgen/pass_roads.h>
#include <mapgen/pass_towns.h>

class Generator {
public:
    Generator(Config& config)
        : config_(config) {}

    void loadFromConfig() {
        generateMap();
    }
    void generateMap() {
        std::cout << "Start of generating map" << std::endl;
        generatePoints(config_.seed);
        triangulatePoints();
        borderCheck();
        executePasses();
        std::cout << "End of generating map" << std::endl;
    }
    void executePasses() {
        if (config_.isWater) {
            MapPassWater waterPass;
            waterPass.executePass(centers, corners, edges, 
                    config_.noiseIsland, config_.thresholdWater, config_.thresholdWaterCount);
        }
        if (config_.isCoast) {
            MapPassCoast coastPass;
            coastPass.executePass(centers, corners, edges);
        }
        if (config_.isElevation) {
            MapPassElevation elevationPass;
            elevationPass.executePass(centers, corners, edges);
        }
        if (config_.isRivers) {
            MapPassRivers riversPass;
            riversPass.executePass(centers, corners, edges, &config_);
        }
        if (config_.isMoisture) {
            MapPassMoisture moisturePass;
            moisturePass.executePass(centers, corners, edges, &config_);
        }
        if (config_.isBiomes) {
            MapPassBiomes biomesPass;
            biomesPass.executePass(centers, corners, edges);
        }
        if (config_.isRoads) {
            MapPassRoads roadsPass;
            roadsPass.executePass(centers, corners, edges);
        }
        if (config_.isTowns) {
            MapPassTowns townsPass;
            townsPass.executePass(centers, corners, edges);
        }
        if (config_.isNoisyEdges) {
            MapPassNoisyEdges edgesPass;
            edgesPass.executePass(centers, corners, edges);
        }
    }
    void saveToPNG(const char* filename) {
        MapDrawer::saveToPNG(filename, 
            config_.gridSize, config_.imageSize,
            centers, &config_);
    }

private:
    Config& config_;
    std::vector<MapPoint> points_;
    std::vector<std::shared_ptr<MapCenter>> centers;
    std::vector<std::shared_ptr<MapCorner>> corners;
    std::vector<std::shared_ptr<MapEdge>> edges;

    Config getDefaultConfig() {
        Config defaultConfig;
        // ... set any specific default values if needed ...
        return defaultConfig;
    }

    // Generate points with optional jitter
    void generatePoints(int seed) {
        std::mt19937 gen(seed);
        std::uniform_real_distribution<> dist(-1.0, 1.0);
        points_.clear();

        for (unsigned int x = 0; x <= config_.gridSize; x++) {
            for (unsigned int y = 0; y <= config_.gridSize; y++) {
                double jitterX = x + config_.jitter * (dist(gen) - 0.5);
                double jitterY = y + config_.jitter * (dist(gen) - 0.5);
                points_.emplace_back(MapPoint{jitterX, jitterY});

                // std::cout << "emplacing point at: " << jitterX << ", " << jitterY << std::endl;
            }
        }

         // Adding boundary points to extend coverage beyond the grid
        double extra = 1.0;  // Extend by one unit beyond the grid size
        for (unsigned int x = 0; x <= config_.gridSize; x++) {
            // Top and bottom boundary
            points_.push_back({x - extra, -extra}); // Top outside
            points_.push_back({x + extra, config_.gridSize + extra}); // Bottom outside
        }
        for (unsigned int y = 0; y <= config_.gridSize; y++) {
            // Left and right boundary
            points_.push_back({-extra, y - extra}); // Left outside
            points_.push_back({config_.gridSize + extra, y + extra}); // Right outside
        }

        // Add corners if not already included
        points_.push_back({-extra, -extra});
        points_.push_back({-extra, config_.gridSize + extra});
        points_.push_back({config_.gridSize + extra, -extra});
        points_.push_back({config_.gridSize + extra, config_.gridSize + extra});
    }
    // Perform Delaunay triangulation and prepare for Voronoi construction
    void triangulatePoints() {
        std::vector<double> coords;
        for (const auto& pt : points_) {
            coords.push_back(pt.x);
            coords.push_back(pt.y);
        }
        delaunator::Delaunator d(coords);
        collectDelaunayData(d);
    }
    // Helper function to calculate the centroid of a polygon formed by points
    MapPoint calculateCentroid(const std::vector<std::shared_ptr<MapCorner>>& corners) {
        double centroidX = 0, centroidY = 0;
        for (const auto& corner : corners) {
            centroidX += corner->point.x;
            centroidY += corner->point.y;
        }
        size_t count = corners.size();
        return {centroidX / count, centroidY / count};
    }
    // Function to sort corners radially around their centroid
    void sortCorners(std::vector<std::shared_ptr<MapCorner>>& corners) {
        MapPoint centroid = calculateCentroid(corners);
        std::sort(corners.begin(), corners.end(), [centroid](const std::shared_ptr<MapCorner>& a, const std::shared_ptr<MapCorner>& b) {
            double angleA = std::atan2(a->point.y - centroid.y, a->point.x - centroid.x);
            double angleB = std::atan2(b->point.y - centroid.y, b->point.x - centroid.x);
            return angleA < angleB;
        });
    }
    void borderCheck() {
        for (const auto& center : centers) {
            for (const auto& corner : center->corners) {
                if (corner->point.x < 0 || corner->point.x > config_.gridSize ||
                    corner->point.y < 0 || corner->point.y > config_.gridSize) {
                    corner->border = true;
                    center->border = true;
                    break;
                }
            }
        }
    }
    // Collect data from Delaunay triangulation to construct both Delaunay and Voronoi diagrams
    void collectDelaunayData(const delaunator::Delaunator& delaunay) {
        centers.resize(points_.size());
        corners.clear();
        edges.clear();
        int edgeIndex = 0;
        int cornerIndex = 0;

        // Initialize centers
        for (size_t i = 0; i < points_.size(); ++i) {
            auto center = std::make_shared<MapCenter>();
            center->index = i;
            center->point = points_[i];
            centers[i] = center;
        }

        // Map to store circumcenters by triangle index to link half-edges later
        std::unordered_map<std::size_t, std::shared_ptr<MapCorner>> triangleCircumcenters;

        // Calculate corners and link to centers
        for (std::size_t i = 0; i < delaunay.triangles.size(); i += 3) {
            std::size_t i0 = delaunay.triangles[i];
            std::size_t i1 = delaunay.triangles[i + 1];
            std::size_t i2 = delaunay.triangles[i + 2];

            auto circumcenter = calculateCircumcenter(
                points_[i0].x, points_[i0].y,
                points_[i1].x, points_[i1].y,
                points_[i2].x, points_[i2].y
            );

            auto corner = std::make_shared<MapCorner>();
            corner->point = {circumcenter[0], circumcenter[1]};
            corner->index = cornerIndex++;
            corners.push_back(corner);
            triangleCircumcenters[i / 3] = corner;

            // Assign corner to each center involved in this triangle
            centers[i0]->corners.push_back(corner);
            centers[i1]->corners.push_back(corner);
            centers[i2]->corners.push_back(corner);
            
            corner->touches.push_back(centers[i0]);
            corner->touches.push_back(centers[i1]);
            corner->touches.push_back(centers[i2]);
        }

        for (auto& center : centers) {
            //sortCorners(center->corners);
        }

        // Create half-edge information and connect corners
        for (size_t e = 0; e < delaunay.halfedges.size(); ++e) {
            if (delaunay.halfedges[e] == -1) continue; // Skip edge if it does not have a pair

            size_t e_pair = delaunay.halfedges[e];
            auto edge = std::make_shared<MapEdge>();
            edge->index = edgeIndex++;

            edge->v0 = triangleCircumcenters[e / 3]; // From current triangle
            edge->v1 = triangleCircumcenters[e_pair / 3]; // To paired triangle

            edge->midpoint = {(edge->v0->point.x + edge->v1->point.x) / 2, (edge->v0->point.y + edge->v1->point.y) / 2};

            edges.push_back(edge);

            // Connect edges to centers and corners
            centers[delaunay.triangles[e_pair]]->borders.push_back(edge);
            centers[delaunay.triangles[e]]->borders.push_back(edge);

            edge->v0->protrudes.push_back(edge);
            edge->v1->protrudes.push_back(edge);

            edge->v0->adjacent.push_back(edge->v1);
            edge->v1->adjacent.push_back(edge->v0);

            edge->d0 = centers[delaunay.triangles[e]];
            edge->d1 = centers[delaunay.triangles[e_pair]];
            
            centers[delaunay.triangles[e]]->neighbors.push_back(centers[delaunay.triangles[e_pair]]);
            centers[delaunay.triangles[e_pair]]->neighbors.push_back(centers[delaunay.triangles[e]]);
        }
    }
    std::array<double, 2> calculateCircumcenter(double x1, double y1, double x2, double y2, double x3, double y3) {
        double D = 2 * (x1*(y2 - y3) + x2*(y3 - y1) + x3*(y1 - y2));
        double Ux = ((x1*x1 + y1*y1) * (y2 - y3) + (x2*x2 + y2*y2) * (y3 - y1) + (x3*x3 + y3*y3) * (y1 - y2)) / D;
        double Uy = ((x1*x1 + y1*y1) * (x3 - x2) + (x2*x2 + y2*y2) * (x1 - x3) + (x3*x3 + y3*y3) * (x2 - x1)) / D;
        return {Ux, Uy};
    }
};

#endif // GENERATOR_H
