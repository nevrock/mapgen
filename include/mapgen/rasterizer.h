#ifndef DRAWER_BIOMES_H
#define DRAWER_BIOMES_H

#include <vector>
#include <random>
#include <iostream>
#include <glm/glm.hpp>
#include <mapgen/data.h>
#include <mapgen/config.h>
#include <mapgen/drawer.h>

struct Pixel {
    glm::vec3 color;
    bool hasRiver;
    bool hasRoad;

    Pixel() : color(255.0f, 255.0f, 255.0f), hasRiver(false), hasRoad(false) {} // white background
};

class MapRasterizer {
public:
    MapRasterizer() {}

    static glm::vec3 getColorFromConfig(const std::string& colorName, const Config& config) {
        if (colorName == "ocean") return config.ocean;
        if (colorName == "lake") return config.lake;
        if (colorName == "coast") return config.coast;
        if (colorName == "ground") return config.ground;
        if (colorName == "marsh") return config.marsh;
        if (colorName == "ice") return config.ice;
        if (colorName == "beach") return config.beach;
        if (colorName == "snow") return config.snow;
        if (colorName == "tundra") return config.tundra;
        if (colorName == "bare") return config.bare;
        if (colorName == "scorched") return config.scorched;
        if (colorName == "taiga") return config.taiga;
        if (colorName == "shrubland") return config.shrubland;
        if (colorName == "temperate_desert") return config.temperateDesert;
        if (colorName == "temperate_rain_forest") return config.temperateRainForest;
        if (colorName == "temperate_decidious_forest") return config.temperateDecidiousForest;
        if (colorName == "grassland") return config.grassland;
        if (colorName == "tropical_rain_forest") return config.tropicalRainForest;
        if (colorName == "tropical_seasonal_forest") return config.tropicalSeasonalForest;
        if (colorName == "subtropical_desert") return config.subtropicalDesert;

        return glm::vec3(0.0); 
    }

    static void rasterize(unsigned int gridSize, unsigned int rasterSize, 
                    unsigned int startX, unsigned int startY,
                    std::vector<std::shared_ptr<MapCenter>>& centers, Config* config,
                    std::vector<std::vector<Pixel>>& pixels) {
        int imageWidth = rasterSize;
        int imageHeight = rasterSize;

        std::mt19937 rng(std::random_device{}());
        std::uniform_real_distribution<> dist(0, 1);

        std::vector<glm::vec3> centerColors;
        for (const auto& center : centers) {
            glm::vec3 color = getColorFromConfig(center->biome, *config);
            centerColors.push_back(color);
        }

        int i = 0;
        for (auto& center : centers) {
            std::vector<MapPoint> vertices;
            for (auto& edge : center->borders) {
                for (auto& p : edge->noisyPoints0) {
                    double normalizedX = (p.x / gridSize) * rasterSize;
                    double normalizedY = (p.y / gridSize) * rasterSize;
                    vertices.push_back({normalizedX, normalizedY});
                }
            }

            sortPoints(vertices);
            fillPolygon(vertices, centerColors[i], pixels, imageWidth, imageHeight, gridSize, center);
            i++;
        }

        for (auto& center : centers) {
            for (auto& edge : center->borders) {
                if (edge->v1->ocean) {
                    continue;
                }
                if (edge->river > 0) {
                    drawRiver(edge, pixels, imageWidth, imageHeight, edge->river + config->riverFactor, gridSize, config);
                }
            }
        }

        for (auto& center : centers) {
            for (auto& edge : center->borders) {
                if (edge->road) {
                    drawRoad(edge, pixels, imageWidth, imageHeight, config->roadSize, gridSize, config);
                }
            }
        }
    }

private:
    static MapPoint calculateCentroid(const std::vector<MapPoint>& points) {
        double centroidX = 0, centroidY = 0;
        for (const auto& point : points) {
            centroidX += point.x;
            centroidY += point.y;
        }
        size_t count = points.size();
        return {centroidX / count, centroidY / count};
    }

    static void sortPoints(std::vector<MapPoint>& points) {
        MapPoint centroid = calculateCentroid(points);
        std::sort(points.begin(), points.end(), [centroid](const MapPoint& a, const MapPoint& b) {
            double angleA = std::atan2(a.y - centroid.y, a.x - centroid.x);
            double angleB = std::atan2(b.y - centroid.y, b.x - centroid.x);
            return angleA < angleB;
        });
    }

    static void fillPolygon(const std::vector<MapPoint>& vertices, 
        const glm::vec3& color, std::vector<unsigned char>& image, 
        unsigned int imageWidth, unsigned int imageHeight, unsigned int gridSize, 
        const std::shared_ptr<MapCenter>& center) {
        if (vertices.empty()) return;

        MapPoint centroid = calculateCentroid(vertices);
        size_t n = vertices.size();

        for (size_t i = 0; i < n; ++i) {
            std::vector<MapPoint> triangle = {centroid, vertices[i], vertices[(i + 1) % n]};
            fillTriangle(triangle, color, image, imageWidth, imageHeight, gridSize, center);
        }
    }
    static void fillTriangle(const std::vector<MapPoint>& triangle, 
        const glm::vec3& color, std::vector<unsigned char>& image, 
        unsigned int imageWidth, unsigned int imageHeight, unsigned int gridSize, 
        const std::shared_ptr<MapCenter>& center) {
        if (triangle.size() < 3) return; // Safety check

        auto minmaxX = std::minmax_element(triangle.begin(), triangle.end(), [](const MapPoint& a, const MapPoint& b) { return a.x < b.x; });
        auto minmaxY = std::minmax_element(triangle.begin(), triangle.end(), [](const MapPoint& a, const MapPoint& b) { return a.y < b.y; });

        int minX = std::max(0, static_cast<int>(floor(minmaxX.first->x)));
        int maxX = std::min(static_cast<int>(imageWidth - 1), static_cast<int>(ceil(minmaxX.second->x)));
        int minY = std::max(0, static_cast<int>(floor(minmaxY.first->y)));
        int maxY = std::min(static_cast<int>(imageHeight - 1), static_cast<int>(ceil(minmaxY.second->y)));

        auto edgeFunction = [](const MapPoint& a, const MapPoint& b, const MapPoint& c) {
            return (c.x - a.x) * (b.y - a.y) - (c.y - a.y) * (b.x - a.x);
        };

        const MapPoint& v0 = triangle[0];
        const MapPoint& v1 = triangle[1];
        const MapPoint& v2 = triangle[2];

        float area = edgeFunction(v0, v1, v2);
        bool hasPositiveArea = area > 0;

        for (int y = minY; y <= maxY; y++) {
            for (int x = minX; x <= maxX; x++) {
                MapPoint pt = {static_cast<double>(x), static_cast<double>(y)};
                float w0 = edgeFunction(v1, v2, pt);
                float w1 = edgeFunction(v2, v0, pt);
                float w2 = edgeFunction(v0, v1, pt);

                if ((w0 == 0 || w1 == 0 || w2 == 0) || (hasPositiveArea ? (w0 > 0 && w1 > 0 && w2 > 0) : (w0 < 0 && w1 < 0 && w2 < 0))) {
                    double gridX = pt.x * gridSize / imageWidth;
                    double gridY = pt.y * gridSize / imageHeight;
                    double elevation = center->getElevation(gridX, gridY);
                    
                    int index = (y * imageWidth + x) * 3;
                    image[index]     = static_cast<unsigned char>(color.r);
                    image[index + 1] = static_cast<unsigned char>(color.g);
                    image[index + 2] = static_cast<unsigned char>(color.b);
                }
            }
        }
    }
    static void drawRiver(const std::shared_ptr<MapEdge>& edge, std::vector<std::vector<Pixel>>& pixels, 
                      unsigned int imageWidth, unsigned int imageHeight, int riverSize, int gridSize, Config* config) {
        int x0 = static_cast<int>((edge->v0->point.x / gridSize) * imageWidth);
        int y0 = static_cast<int>((edge->v0->point.y / gridSize) * imageHeight);
        int x1 = static_cast<int>((edge->v1->point.x / gridSize) * imageWidth);
        int y1 = static_cast<int>((edge->v1->point.y / gridSize) * imageHeight);

        int dx = std::abs(x1 - x0), sx = x0 < x1 ? 1 : -1;
        int dy = -std::abs(y1 - y0), sy = y0 < y1 ? 1 : -1;
        int err = dx + dy, e2;
        
        glm::vec3 riverColor = config->riverColor;

        while (true) {
            for (int i = -riverSize; i <= riverSize; ++i) {
                for (int j = -riverSize; j <= riverSize; ++j) {
                    int nx = x0 + i;
                    int ny = y0 + j;
                    if (nx >= 0 && nx < imageWidth && ny >= 0 && ny < imageHeight) {
                        pixels[ny][nx].color = riverColor;
                        pixels[ny][nx].hasRiver = true;
                    }
                }
            }

            if (x0 == x1 && y0 == y1) break;
            e2 = 2 * err;
            if (e2 >= dy) { err += dy; x0 += sx; }
            if (e2 <= dx) { err += dx; y0 += sy; }
        }
    }

    static void drawRoad(const std::shared_ptr<MapEdge>& edge, std::vector<std::vector<Pixel>>& pixels, 
                        unsigned int imageWidth, unsigned int imageHeight, int roadSize, int gridSize, Config* config) {
        int x0 = static_cast<int>((edge->d0->point.x / gridSize) * imageWidth);
        int y0 = static_cast<int>((edge->d0->point.y / gridSize) * imageHeight);
        int x1 = static_cast<int>((edge->d1->point.x / gridSize) * imageWidth);
        int y1 = static_cast<int>((edge->d1->point.y / gridSize) * imageHeight);

        int dx = std::abs(x1 - x0), sx = x0 < x1 ? 1 : -1;
        int dy = -std::abs(y1 - y0), sy = y0 < y1 ? 1 : -1;
        int err = dx + dy, e2;
        
        glm::vec3 roadColor = config->roadColor;

        while (true) {
            for (int i = -roadSize; i <= roadSize; ++i) {
                for (int j = -roadSize; j <= roadSize; ++j) {
                    int nx = x0 + i;
                    int ny = y0 + j;
                    if (nx >= 0 && nx < imageWidth && ny >= 0 && ny < imageHeight) {
                        pixels[ny][nx].color = roadColor;
                        pixels[ny][nx].hasRoad = true;
                    }
                }
            }

            if (x0 == x1 && y0 == y1) break;
            e2 = 2 * err;
            if (e2 >= dy) { err += dy; x0 += sx; }
            if (e2 <= dx) { err += dx; y0 += sy; }
        }
    }

};

#endif // DRAWER_BIOMES_H



#ifndef DRAWER_BIOMES_H
#define DRAWER_BIOMES_H

#include <vector>
#include <random>
#include <iostream>
#include <glm/glm.hpp>
#include <mapgen/data.h>
#include <mapgen/config.h>
#include <mapgen/drawer.h>

struct Pixel {
    glm::vec3 color;
    bool hasRiver;
    bool hasRoad;

    Pixel() : color(255.0f, 255.0f, 255.0f), hasRiver(false), hasRoad(false) {} // white background
};

class MapRasterizer {
public:
    MapRasterizer() {}

    static glm::vec3 getColorFromConfig(const std::string& colorName, const Config& config) {
        if (colorName == "ocean") return config.ocean;
        if (colorName == "lake") return config.lake;
        if (colorName == "coast") return config.coast;
        if (colorName == "ground") return config.ground;
        if (colorName == "marsh") return config.marsh;
        if (colorName == "ice") return config.ice;
        if (colorName == "beach") return config.beach;
        if (colorName == "snow") return config.snow;
        if (colorName == "tundra") return config.tundra;
        if (colorName == "bare") return config.bare;
        if (colorName == "scorched") return config.scorched;
        if (colorName == "taiga") return config.taiga;
        if (colorName == "shrubland") return config.shrubland;
        if (colorName == "temperate_desert") return config.temperateDesert;
        if (colorName == "temperate_rain_forest") return config.temperateRainForest;
        if (colorName == "temperate_decidious_forest") return config.temperateDecidiousForest;
        if (colorName == "grassland") return config.grassland;
        if (colorName == "tropical_rain_forest") return config.tropicalRainForest;
        if (colorName == "tropical_seasonal_forest") return config.tropicalSeasonalForest;
        if (colorName == "subtropical_desert") return config.subtropicalDesert;

        return glm::vec3(0.0); 
    }

    static void rasterize(unsigned int gridSize, unsigned int rasterSize, 
                    unsigned int startX, unsigned int startY,
                    std::vector<std::shared_ptr<MapCenter>>& centers, Config* config,
                    std::vector<std::vector<Pixel>>& pixels) {
        int imageWidth = rasterSize;
        int imageHeight = rasterSize;

        std::mt19937 rng(std::random_device{}());
        std::uniform_real_distribution<> dist(0, 1);

        std::vector<glm::vec3> centerColors;
        for (const auto& center : centers) {
            glm::vec3 color = getColorFromConfig(center->biome, *config);
            centerColors.push_back(color);
        }

        int i = 0;
        for (auto& center : centers) {
            std::vector<MapPoint> vertices;
            for (auto& edge : center->borders) {
                for (auto& p : edge->noisyPoints0) {
                    double normalizedX = (p.x / gridSize) * rasterSize;
                    double normalizedY = (p.y / gridSize) * rasterSize;
                    vertices.push_back({normalizedX, normalizedY});
                }
            }

            sortPoints(vertices);
            fillPolygon(vertices, centerColors[i], pixels, imageWidth, imageHeight, gridSize, center);
            i++;
        }

        for (auto& center : centers) {
            for (auto& edge : center->borders) {
                if (edge->v1->ocean) {
                    continue;
                }
                if (edge->river > 0) {
                    drawRiver(edge, pixels, imageWidth, imageHeight, edge->river + config->riverFactor, gridSize, config);
                }
            }
        }

        for (auto& center : centers) {
            for (auto& edge : center->borders) {
                if (edge->road) {
                    drawRoad(edge, pixels, imageWidth, imageHeight, config->roadSize, gridSize, config);
                }
            }
        }
    }

private:
    static MapPoint calculateCentroid(const std::vector<MapPoint>& points) {
        double centroidX = 0, centroidY = 0;
        for (const auto& point : points) {
            centroidX += point.x;
            centroidY += point.y;
        }
        size_t count = points.size();
        return {centroidX / count, centroidY / count};
    }

    static void sortPoints(std::vector<MapPoint>& points) {
        MapPoint centroid = calculateCentroid(points);
        std::sort(points.begin(), points.end(), [centroid](const MapPoint& a, const MapPoint& b) {
            double angleA = std::atan2(a.y - centroid.y, a.x - centroid.x);
            double angleB = std::atan2(b.y - centroid.y, b.x - centroid.x);
            return angleA < angleB;
        });
    }

    static void fillPolygon(const std::vector<MapPoint>& vertices, 
                            const glm::vec3& color, std::vector<std::vector<Pixel>>& pixels, 
                            unsigned int imageWidth, unsigned int imageHeight, 
                            unsigned int gridSize, const std::shared_ptr<MapCenter>& center) {
        if (vertices.empty()) return;

        MapPoint centroid = calculateCentroid(vertices);
        size_t n = vertices.size();

        for (size_t i = 0; i < n; ++i) {
            std::vector<MapPoint> triangle = {centroid, vertices[i], vertices[(i + 1) % n]};
            fillTriangle(triangle, color, pixels, imageWidth, imageHeight, gridSize, center);
        }
    }

    static void fillTriangle(const std::vector<MapPoint>& triangle, 
                             const glm::vec3& color, std::vector<std::vector<Pixel>>& pixels, 
                             unsigned int imageWidth, unsigned int imageHeight, 
                             unsigned int gridSize, const std::shared_ptr<MapCenter>& center) {
        if (triangle.size() < 3) return; // Safety check

        auto minmaxX = std::minmax_element(triangle.begin(), triangle.end(), [](const MapPoint& a, const MapPoint& b) { return a.x < b.x; });
        auto minmaxY = std::minmax_element(triangle.begin(), triangle.end(), [](const MapPoint& a, const MapPoint& b) { return a.y < b.y; });

        int minX = std::max(0, static_cast<int>(std::floor(minmaxX.first->x)));
        int maxX = std::min(static_cast<int>(imageWidth - 1), static_cast<int>(std::ceil(minmaxX.second->x)));
        int minY = std::max(0, static_cast<int>(std::floor(minmaxY.first->y)));
        int maxY = std::min(static_cast<int>(imageHeight - 1), static_cast<int>(std::ceil(minmaxY.second->y)));

        auto edgeFunction = [](const MapPoint& a, const MapPoint& b, const MapPoint& c) {
            return (c.x - a.x) * (b.y - a.y) - (c.y - a.y) * (b.x - a.x);
        };

        const MapPoint& v0 = triangle[0];
        const MapPoint& v1 = triangle[1];
        const MapPoint& v2 = triangle[2];

        float area = edgeFunction(v0, v1, v2);
        bool hasPositiveArea = area > 0;

        for (int y = minY; y <= maxY; y++) {
            for (int x = minX; x <= maxX; x++) {
                MapPoint pt = {static_cast<double>(x), static_cast<double>(y)};
                float w0 = edgeFunction(v1, v2, pt);
                float w1 = edgeFunction(v2, v0, pt);
                float w2 = edgeFunction(v0, v1, pt);

                if ((w0 == 0 || w1 == 0 || w2 == 0) || (hasPositiveArea ? (w0 > 0 && w1 > 0 && w2 > 0) : (w0 < 0 && w1 < 0 && w2 < 0))) {
                    if (x >= 0 && x < imageWidth && y >= 0 && y < imageHeight) {
                        pixels[y][x].color = color;
                    }
                }
            }
        }
    }

    static void drawRiver(const std::shared_ptr<MapEdge>& edge, std::vector<std::vector<Pixel>>& pixels, unsigned int imageWidth, unsigned int imageHeight, int riverSize, int gridSize, Config*
