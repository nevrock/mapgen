#ifndef MAP_PASS_TOWNS_H
#define MAP_PASS_TOWNS_H

#include <mapgen/data.h>
#include <mapgen/config.h>
#include <snorri/log.h>
#include <snorri/noise.h>
#include <vector>
#include <queue>
#include <unordered_set>
#include <random>  // Include for std::random_device and std::mt19937

class MapPassTowns {
public:
    MapPassTowns() {}

    void executePass(std::vector<std::shared_ptr<MapCenter>>& centers,
                                    std::vector<std::shared_ptr<MapCorner>>& corners,
                                    std::vector<std::shared_ptr<MapEdge>>& edges) {
        Log::console("map pass towns!");
        
    }

private:
    /*
    std::vector<Vector2> packBuildings() {
        for (const Building& building : buildings) {
            bool placed = false;
            for (float y = 0; y < 100 && !placed; y += 1.0f) {
                for (float x = 0; x < 100 && !placed; x += 1.0f) {
                    Vector2 center = {x, y};
                    if (canPlaceBuilding(building, center)) {
                        packedPositions.push_back(center);
                        placed = true;
                    }
                }
            }
        }
        return packedPositions;
    }
    bool isPointInsidePolygon(const Vector2& point) const {
        int intersections = 0;
        for (size_t i = 0; i < polygon.size(); i++) {
            Vector2 p1 = polygon[i];
            Vector2 p2 = polygon[(i + 1) % polygon.size()];

            if ((point.y > std::min(p1.y, p2.y)) && (point.y <= std::max(p1.y, p2.y)) && 
                (point.x <= std::max(p1.x, p2.x)) && (p1.y != p2.y)) {
                float xinters = (point.y - p1.y) * (p2.x - p1.x) / (p2.y - p1.y) + p1.x;
                if (p1.x == p2.x || point.x <= xinters) {
                    intersections++;
                }
            }
        }
        return intersections % 2 != 0;
    }
    bool canPlaceBuilding(const Building& building, const Vector2& center) const {
        Vector2 halfSize = {building.width / 2, building.height / 2};
        Vector2 bottomLeft = {center.x - halfSize.x, center.y - halfSize.y};
        Vector2 bottomRight = {center.x + halfSize.x, center.y - halfSize.y};
        Vector2 topLeft = {center.x - halfSize.x, center.y + halfSize.y};
        Vector2 topRight = {center.x + halfSize.x, center.y + halfSize.y};

        return isPointInsidePolygon(bottomLeft) && isPointInsidePolygon(bottomRight) &&
               isPointInsidePolygon(topLeft) && isPointInsidePolygon(topRight);
    }
    */
};

#endif // MAP_PASS_TOWNS_H
