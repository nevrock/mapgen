#ifndef MAP_PASS_ELEVATION_H
#define MAP_PASS_ELEVATION_H

#include <mapgen/data.h>
#include <mapgen/config.h>
#include <snorri/log.h>
#include <queue>
#include <limits>
#include <cmath>
#include <unordered_map> // Add this include directive

class MapPassElevation {
public:
    MapPassElevation() {}

    void executePass(std::vector<std::shared_ptr<MapCenter>>& centers,
                 std::vector<std::shared_ptr<MapCorner>>& corners,
                 std::vector<std::shared_ptr<MapEdge>>& edges) {
        std::queue<std::shared_ptr<MapCorner>> queue;

        Log::console("graph pass elevation!");

        // Initialize corners' elevation based on border status
        for (auto& corner : corners) {
            if (corner->border) {
                corner->elevation = 0.0;
                queue.push(corner);
            } else {
                corner->elevation = std::numeric_limits<double>::infinity();
            }
        }

        // Assign elevations based on distance from borders
        while (!queue.empty()) {
            auto current = queue.front();
            queue.pop();

            for (auto& neighbor : current->adjacent) {
                double newElevation = 0.01 + current->elevation;
                if (!current->water && !neighbor->water) {
                    newElevation += 1.0;
                    //if (needsMoreRandomness) {
                    //newElevation += config->getRandomDouble();  // Assuming Dict has a method to provide a random double
                    //}
                }

                if (newElevation < neighbor->elevation) {
                    neighbor->elevation = newElevation;
                    queue.push(neighbor);
                }
            }
        }

        redistributeElevations(corners);
        calculateDownslopes(corners);
        calculateCenterElevations(centers);
    }
    void redistributeElevations(std::vector<std::shared_ptr<MapCorner>>& corners) {
        const double SCALE_FACTOR = 1.1;
        std::vector<double> sortedElevations(corners.size());
        std::vector<double> newElevations(corners.size());

        // Copy elevations to a separate vector for sorting
        for (size_t i = 0; i < corners.size(); ++i) {
            sortedElevations[i] = corners[i]->elevation;
        }

        // Sort the temporary vector of elevations
        std::sort(sortedElevations.begin(), sortedElevations.end());

        // Calculate new elevations based on sorted values
        for (size_t i = 0; i < sortedElevations.size(); ++i) {
            double y = static_cast<double>(i) / (sortedElevations.size() - 1);
            double x = std::sqrt(SCALE_FACTOR) - std::sqrt(SCALE_FACTOR * (1 - y));
            if (x > 1.0) x = 1.0;  // Clamp the value to not exceed 1.0
            newElevations[i] = x;
        }

        // Map new elevations back to the original corners without altering their order
        std::unordered_map<double, std::queue<double>> elevationMap;
        for (double elevation : sortedElevations) {
            elevationMap[elevation].push(newElevations.front());
            newElevations.erase(newElevations.begin());
        }

        for (auto& corner : corners) {
            double originalElevation = corner->elevation;
            if (!elevationMap[originalElevation].empty()) {
                corner->elevation = elevationMap[originalElevation].front();
                elevationMap[originalElevation].pop();
            }
        }
    }
    void calculateCenterElevations(std::vector<std::shared_ptr<MapCenter>>& centers) {
        // Assign elevations to centers based on the average of corner elevations
        for (auto& center : centers) {
            double sumElevation = 0.0;
            for (auto& corner : center->corners) {
                sumElevation += corner->elevation;
            }
            center->elevation = sumElevation / center->corners.size();
        }
    }
    void calculateDownslopes(std::vector<std::shared_ptr<MapCorner>>& corners) {
        // Iterate over each corner to determine the downslope corner
        for (auto& q : corners) {
            std::shared_ptr<MapCorner> r = q; // Start by assuming the current corner is the lowest
            // Check all adjacent corners to find the one with the lowest elevation
            for (auto& s : q->adjacent) {
                if (s->elevation <= r->elevation) {
                    r = s;
                }
            }
            // Set the downslope corner
            q->downslope = r;
        }
    }
};

#endif // MAP_PASS_ELEVATION_H
