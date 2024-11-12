#include <iostream>
#include <fstream>
#include <string>
#include <map>

#include <mapgen/config.h>
#include <mapgen/generator.h>


void GenerateMap(
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

int main(int argc, char* argv[]) {

    std::cout << "Map generator from cli!" << std::endl;

    // Default values (same as your original code)
    std::string filepath = "C:\\Users\\trave\\Documents\\mapgen\\file_out"; 
    int seed = 251;
    int imageSize = 1024;
    double wavelength = 20.0;
    int gridSize = 50;
    double jitter = 0.5;
    double thresholdWater = 0.3; 
    int riverCount = 25;
    int roadSize = 2;
    int riverFactor = 1;

    // Process command line arguments
    for (int i = 1; i < argc; i += 2) { 
        std::string argName = argv[i];

        // Make sure there's a value after the argument name
        if (i + 1 < argc) {
            std::string argValue = argv[i + 1];

            if (argName == "filepath") {
                filepath = argValue;
            } else if (argName == "seed") {
                seed = std::stoi(argValue); 
            } else if (argName == "imageSize") {
                imageSize = std::stoi(argValue);
            } else if (argName == "wavelength") {
                wavelength = std::stod(argValue); 
            } else if (argName == "gridSize") {
                gridSize = std::stoi(argValue);
            } else if (argName == "jitter") {
                jitter = std::stod(argValue);
            } else if (argName == "thresholdWater") {
                thresholdWater = std::stod(argValue);
            } else if (argName == "riverCount") {
                riverCount = std::stoi(argValue);
            } else if (argName == "roadSize") {
                roadSize = std::stoi(argValue);
            } else if (argName == "riverFactor") {
                riverFactor = std::stoi(argValue);
            } else {
                std::cerr << "Unknown argument: " << argName << std::endl;
            }
        } else {
            std::cerr << "Missing value for argument: " << argName << std::endl;
        }
    }

    // Now call GenerateMap with the potentially modified values
    GenerateMap(filepath.c_str(), seed, imageSize, wavelength, gridSize, jitter, 
                thresholdWater, riverCount, roadSize, riverFactor);

    return 0;
}

