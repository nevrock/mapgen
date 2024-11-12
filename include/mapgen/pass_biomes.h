#ifndef MAP_PASS_BIOMES_H
#define MAP_PASS_BIOMES_H

#include <mapgen/data.h>
#include <mapgen/config.h>

#include <snorri/log.h>
#include <snorri/noise.h>

#include <vector>
#include <queue>
#include <unordered_set>
#include <random>  // Include for std::random_device and std::mt19937

class MapPassBiomes {
public:
    MapPassBiomes() {}

    void executePass(std::vector<std::shared_ptr<MapCenter>>& centers,
                    std::vector<std::shared_ptr<MapCorner>>& corners,
                    std::vector<std::shared_ptr<MapEdge>>& edges) {
        Log::console("graph pass biomes!");

        assignBiomes(centers);
    }
    void assignBiomes(std::vector<std::shared_ptr<MapCenter>>& centers) {
        for (auto& p : centers) {
            p->biome = getBiome(p);
        }
    }
    static std::string getBiome(const std::shared_ptr<MapCenter>& p) {
        if (p->ocean) {
            return "ocean";
        } else if (p->water) {
            if (p->elevation < 0.1) return "marsh";
            if (p->elevation > 0.8) return "ice";
            return "lake";
        } else if (p->coast) {
            return "beach";
        } else if (p->elevation > 0.8) {
            if (p->moisture > 0.50) return "snow";
            else if (p->moisture > 0.33) return "tundra";
            else if (p->moisture > 0.16) return "bare";
            else return "scorched";
        } else if (p->elevation > 0.6) {
            if (p->moisture > 0.66) return "taiga";
            else if (p->moisture > 0.33) return "shrubland";
            else return "temperate_desert";
        } else if (p->elevation > 0.3) {
            if (p->moisture > 0.83) return "temperate_rain_forest";
            else if (p->moisture > 0.50) return "temperate_decidious_forest";
            else if (p->moisture > 0.16) return "grassland";
            else return "temperate_desert";
        } else {
            if (p->moisture > 0.66) return "tropical_rain_forest";
            else if (p->moisture > 0.33) return "tropical_seasonal_forest";
            else if (p->moisture > 0.16) return "grassland";
            else return "subtropical_desert";
        }
    }
};

#endif // MAP_PASS_BIOMES_H
