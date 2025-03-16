#ifndef CONFIG_H
#define CONFIG_H

#include <map>
#include <string>
#include <vector>
#include <glm/glm.hpp>

// Define a struct to hold the noise island parameters
struct NoiseConfig {
    int seed = 251;
    double frequency = 0.08;
    int type = 0;
    int fractalType = 1;
    int octaves = 5;
    double lacunarity = 2.0;
    double gain = 0.5;
    double weightedStrength = 0.45;
};

// Define a class to hold all the configuration data
class Config {
public:
    bool isLog = false;

    std::string filepath = "";

    // basic settings
    double wavelength = 20.0;
    int gridSize = 50;
    double jitter = 0.5;
    int imageSize = 1024;
    int seed = 251;
    double borderLength = 1.0;

    // noises
    NoiseConfig noiseIsland;

    double thresholdWater = 0.3;
    double thresholdWaterCount = 2;
    int riverCount = 25;
    int roadSize = 1;
    int riverFactor = 1;

    // colors
    glm::vec3 riverColor = glm::vec3(94, 182, 223);
    glm::vec3 roadColor = glm::vec3(0, 0, 0);

    glm::vec3 ocean = glm::vec3(94, 182, 223);
    glm::vec3 lake = glm::vec3(94, 182, 223);
    glm::vec3 coast = glm::vec3(223, 205, 87);
    glm::vec3 ground = glm::vec3(125, 185, 75);
    glm::vec3 marsh = glm::vec3(33, 94, 33);
    glm::vec3 ice = glm::vec3(210, 255, 252);
    glm::vec3 beach = glm::vec3(245, 222, 179);
    glm::vec3 snow = glm::vec3(255, 250, 250);
    glm::vec3 tundra = glm::vec3(169, 169, 169);
    glm::vec3 bare = glm::vec3(201, 180, 155);
    glm::vec3 scorched = glm::vec3(153, 130, 109);
    glm::vec3 taiga = glm::vec3(51, 102, 0);
    glm::vec3 shrubland = glm::vec3(128, 128, 0);
    glm::vec3 temperateDesert = glm::vec3(238, 214, 175);
    glm::vec3 temperateRainForest = glm::vec3(85, 107, 47);
    glm::vec3 temperateDecidiousForest = glm::vec3(34, 139, 34);
    glm::vec3 grassland = glm::vec3(124, 252, 0);
    glm::vec3 tropicalRainForest = glm::vec3(0, 100, 0);
    glm::vec3 tropicalSeasonalForest = glm::vec3(107, 142, 35);
    glm::vec3 subtropicalDesert = glm::vec3(250, 250, 210);

    // bools
    bool isWater = true;
    bool isCoast = true;
    bool isElevation = true;
    bool isRivers = true;
    bool isMoisture = true;
    bool isBiomes = true;
    bool isRoads = true;
    bool isTowns = true;
    bool isNoisyEdges = true;
};

#endif // GRAPH_DRAWER_H
