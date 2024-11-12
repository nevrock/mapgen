#include <iostream>
#include <fstream>
#include <string>
#include <map>

#include <mapgen/config.h>
#include <mapgen/generator.h>
//#include <mapgen/graph_generator.h>

#ifdef __cplusplus
extern "C" {
#endif
    __declspec(dllexport) void GenerateMap(
            const char* filepath,
            
            int seed = 251,
            int imageSize = 1024,

            double wavelength = 20.0,
            int gridSize = 50,
            double jitter = 0.5,

            double thresholdWater = 0.3,
            int riverCount = 25,
            int roadSize = 2,
            int riverFactor = 1
        ) {

        // Create a Config object
        Config config; 

        config.filepath = filepath;

        // Modify the config based on input arguments
        config.seed = seed;
        config.imageSize = imageSize;
        config.wavelength = wavelength;
        config.gridSize = gridSize;
        config.jitter = jitter;
        config.thresholdWater = thresholdWater;
        config.riverCount = riverCount;
        config.roadSize = roadSize;
        config.riverFactor = riverFactor;

        std::ofstream outputFile(std::string(filepath) + ".txt"); 

        if (outputFile.is_open()) {
            outputFile << "filepath: " << config.filepath << std::endl;
            outputFile << "seed: " << config.seed << std::endl;
            outputFile << "imageSize: " << config.imageSize << std::endl;
            outputFile << "wavelength: " << config.wavelength << std::endl;
            outputFile << "gridSize: " << config.gridSize << std::endl;
            outputFile << "jitter: " << config.jitter << std::endl;
            outputFile << "thresholdWater: " << config.thresholdWater << std::endl;
            outputFile << "riverCount: " << config.riverCount << std::endl;
            outputFile << "roadSize: " << config.roadSize << std::endl;
            outputFile << "riverFactor: " << config.riverFactor << std::endl;

            outputFile.close(); 
        } else {
            // Handle the error (e.g., print an error message)
            std::cerr << "Error opening file for writing!" << std::endl;
        }

        Generator generator(config);
        generator.loadFromConfig();
        generator.saveToPNG((std::string(filepath) + ".png").c_str());
    }
#ifdef __cplusplus
}
#endif