#ifndef MAP_PASS_RIVERS_H
#define MAP_PASS_RIVERS_H

#include <mapgen/data.h>
#include <mapgen/config.h>
#include <snorri/log.h>
#include <snorri/noise.h>
#include <vector>
#include <queue>
#include <unordered_set>
#include <random>  // Include for std::random_device and std::mt19937
#include <functional> // Add this include directive

class MapPassRivers {
public:
    MapPassRivers() {}

    void executePass(std::vector<std::shared_ptr<MapCenter>>& centers,
                    std::vector<std::shared_ptr<MapCorner>>& corners,
                    std::vector<std::shared_ptr<MapEdge>>& edges, Config* config) {
        Log::console("graph pass rivers!");

        // Setup random number generator for the corners
        int seed = config->seed;
        std::mt19937 gen(seed);
        std::uniform_int_distribution<> dis(0, corners.size() - 1);

        // Lambda to generate random index within range
        auto randomFunc = [&gen, &dis]() -> int {
            return dis(gen); // This generates a random index within the specified range.
        };

        int riverCount = config->riverCount;

        createRivers(corners, edges, riverCount, randomFunc);
    }

    std::shared_ptr<MapEdge> lookupEdgeFromCorner(std::vector<std::shared_ptr<MapEdge>>& edges, std::shared_ptr<MapCorner> c1, std::shared_ptr<MapCorner> c2) {
        // Find the edge that connects two corners
        for (auto& edge : edges) {
            if ((edge->v0 == c1 && edge->v1 == c2) || (edge->v0 == c2 && edge->v1 == c1)) {
                return edge;
            }
        }
        return nullptr;  // Return null if no edge is found
    }

    void createRivers(std::vector<std::shared_ptr<MapCorner>>& corners, 
    std::vector<std::shared_ptr<MapEdge>>& edges, int size, std::function<int()> randomFunc) {
        // Loop for a number of times based on half the size of the corner list
        for (int i = 0; i < size; i++) {
            auto q = corners[randomFunc()];  // Random corner based on provided random function

            // Continue if the corner is in the ocean or the elevation is not suitable for starting a river
            if (q->ocean || q->elevation < 0.3 || q->elevation > 0.9) {
                i-=1;
                continue;
            }

            // Flow the river downhill from corner to corner
            while (!q->coast) {
                if (q == q->downslope) {
                    break;  // Stop if there's no lower elevation to flow to
                }

                auto edge = lookupEdgeFromCorner(edges, q, q->downslope);  // Retrieve the edge between q and its downslope
                if (!edge) continue;  // If no edge found, skip to the next iteration

                edge->river += 1;  // Increment the river count on the edge
                q->river += 1;     // Increment the river count on the current corner
                q->downslope->river += 1;  // Increment the river count on the downslope corner

                q = q->downslope;  // Move to the next downslope corner
            }
        }
    }
};

#endif // MAP_PASS_RIVERS_H
