#ifndef BOUNDS_H
#define BOUNDS_H

#include <glm/glm.hpp>
#include <algorithm>

class Bounds {
public:
    glm::vec3 min;
    glm::vec3 max;

    // Constructors
    Bounds() : min(glm::vec3(std::numeric_limits<float>::max())), max(glm::vec3(std::numeric_limits<float>::lowest())) {}
    Bounds(const glm::vec3& point) : min(point), max(point) {}
    Bounds(const glm::vec3& min, const glm::vec3& max) : min(min), max(max) {}

    // Expand the bounds to include a point
    void expandToFit(const glm::vec3& point) {
        min = glm::min(min, point);
        max = glm::max(max, point);
    }

    // Expand the bounds to include another bounds
    void expandToFit(const Bounds& other) {
        min = glm::min(min, other.min);
        max = glm::max(max, other.max);
    }

    // Check if a point is inside the bounds
    bool contains(const glm::vec3& point) const {
        return (point.x >= min.x && point.x <= max.x &&
                point.y >= min.y && point.y <= max.y &&
                point.z >= min.z && point.z <= max.z);
    }

    // Check if another bounds intersects with this bounds
    bool intersects(const Bounds& other) const {
        return (min.x <= other.max.x && max.x >= other.min.x &&
                min.y <= other.max.y && max.y >= other.min.y &&
                min.z <= other.max.z && max.z >= other.min.z);
    }

    // Calculate the center of the bounds
    glm::vec3 center() const {
        return (min + max) * 0.5f;
    }

    // Calculate the size (extent) of the bounds
    glm::vec3 size() const {
        return max - min;
    }

    // Calculate the volume of the bounds
    float volume() const {
        glm::vec3 s = size();
        return s.x * s.y * s.z;
    }

    // Calculate the surface area of the bounds
    float surfaceArea() const {
        glm::vec3 s = size();
        return 2.0f * (s.x * s.y + s.x * s.z + s.y * s.z);
    }

    // Check if the bounds are valid (min <= max for all axes)
    bool isValid() const {
        return (min.x <= max.x && min.y <= max.y && min.z <= max.z);
    }

    // Reset bounds to an invalid state
    void reset() {
        min = glm::vec3(std::numeric_limits<float>::max());
        max = glm::vec3(std::numeric_limits<float>::lowest());
    }
};

#endif // BOUNDS_H
