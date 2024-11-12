#ifndef MAP_PASS_ROADS_H
#define MAP_PASS_ROADS_H

#include <mapgen/data.h>
#include <mapgen/config.h>
#include <snorri/log.h>
#include <snorri/noise.h>
#include <vector>
#include <queue>
#include <unordered_set>
#include <random>  // Include for std::random_device and std::mt19937

class MapPassRoads {
public:
    MapPassRoads() {}

    void executePass(std::vector<std::shared_ptr<MapCenter>>& centers,
                                    std::vector<std::shared_ptr<MapCorner>>& corners,
                                    std::vector<std::shared_ptr<MapEdge>>& edges) {
        Log::console("starting road creation process!");

        std::vector<int> centerContour(centers.size(), 999);  // Initialize all contours to a high value
        std::vector<int> cornerContour(corners.size(), 999);  // Initialize corner contours similarly
        std::queue<std::shared_ptr<MapCenter>> queue;

        const std::vector<double> elevationThresholds = {0, 0.05, 0.54, 0.74};

        // Initial assignment for coastal and ocean centers
        for (auto& center : centers) {
            if (center->coast || center->ocean) {
                centerContour[center->index] = 1;
                queue.push(center);
            }
        }

        // Flood fill to assign contour levels
        while (!queue.empty()) {
            auto p = queue.front();
            queue.pop();
            int currentLevel = centerContour[p->index];

            for (auto& neighbor : p->neighbors) {
                int newLevel = currentLevel;
                while (neighbor->elevation > elevationThresholds[newLevel] && newLevel < elevationThresholds.size() - 1) {
                    if (neighbor->water) {
                        // Skip increasing contour level if the neighbor is water
                        break;
                    }
                    newLevel++;
                }
                if (newLevel < centerContour[neighbor->index]) {
                    centerContour[neighbor->index] = newLevel;
                    queue.push(neighbor);
                }
            }
        }

        // Assign contour levels to corners based on the minimum of their adjacent centers
        for (auto& center : centers) {
            for (auto& corner : center->corners) {
                if (cornerContour[corner->index] > centerContour[center->index]) {
                    cornerContour[corner->index] = centerContour[center->index];
                }
            }
        }

        // Mark edges as roads where corners have different contour levels
        for (auto& edge : edges) {
            if (edge->v0 && edge->v1 && cornerContour[edge->v0->index] != cornerContour[edge->v1->index]) {
                edge->road = true;  // Set the road flag to true for the edge
                Log::console("Road created on edge: " + std::to_string(edge->index));
            }
        }

        Log::console("road creation process completed!");
    }
};

#endif // MAP_PASS_ROADS_H
