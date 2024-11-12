#ifndef PASS_WATER_H
#define PASS_WATER_H

#include <mapgen/data.h>

#include <snorri/log.h>
#include <snorri/noise.h>

#include <vector>
#include <queue>
#include <unordered_set>

class MapPassWater {
public:
    MapPassWater() {}

    void executePass(std::vector<std::shared_ptr<MapCenter>>& centers,
                 std::vector<std::shared_ptr<MapCorner>>& corners,
                 std::vector<std::shared_ptr<MapEdge>>& edges, 
                 NoiseConfig& noiseIsland, 
                 float thresholdWater, 
                 int thresholdWaterCount) {
        Log::console("graph pass water!");

        // Initialize water properties based on noise or border status
        for (const auto& corner : corners) {
            double nx = corner->point.x;
            double ny = corner->point.y;
            Noise noise(&noiseIsland); 
            float noiseVal = noise.getNoise(nx, ny);

            corner->water = (noiseVal > thresholdWater) || corner->border;
        }

        // Set water properties for centers based on adjacent corners
        for (const auto& center : centers) {
            if (center->border) {
                center->water = true;
                continue;
            }
            int waterCount = 0;
            for (const auto& corner : center->corners) {
                if (corner->water || corner->border) {
                    waterCount++;
                }
            }
            // Log::console("graph water pass, center water count: " + std::to_string(waterCount));
            if (waterCount > thresholdWaterCount) {
                center->water = true;
                center->ocean = false;  // Assume it's a lake initially
            }
        }

        // Flood fill to identify ocean centers starting from border water centers
        std::queue<std::shared_ptr<MapCenter>> queue;
        for (const auto& center : centers) {
            if (center->border && center->water) {
                center->ocean = true;
                queue.push(center);
            }
        }

        while (!queue.empty()) {
            auto current = queue.front();
            queue.pop();
            for (const auto& neighbor : current->neighbors) {
                if (neighbor->water && !neighbor->ocean) {
                    neighbor->ocean = true;
                    queue.push(neighbor);
                }
            }
        }
    }


private:
   
};

#endif // PASS_WATER_H
