#ifndef BOUNDING_BOX_H
#define BOUNDING_BOX_H

#include "point.h"
#include "../utils/math_util.h"

/*
    BoundingBox class represents an axis-aligned bounding box defined by its minimum and maximum coordinates.
    It provides methods to create a bounding box from a set of points, combine two bounding boxes, and check for intersection.
*/
class BoundingBox {
private:
    // Coordinates of the bounding box
    int64_t minX_, minY_, maxX_, maxY_;
public:
    // Constructors
    BoundingBox() : minX_(0), minY_(0), maxX_(0), maxY_(0) {}
    BoundingBox(int64_t minX, int64_t minY, int64_t maxX, int64_t maxY) : minX_(minX), minY_(minY), maxX_(maxX), maxY_(maxY) {}
    BoundingBox(const std::vector<std::shared_ptr<geom::Point>>& points) {
        if (points.empty()) {
            minX_ = minY_ = maxX_ = maxY_ = 0;
            return;
        }

        minX_ = maxX_ = points[0]->X();
        minY_ = maxY_ = points[0]->Y();

        for (const auto& p : points) {
            minX_ = MathUtil::min<int64_t>(minX_, p->X());
            minY_ = MathUtil::min<int64_t>(minY_, p->Y());
            maxX_ = MathUtil::max<int64_t>(maxX_, p->X());
            maxY_ = MathUtil::max<int64_t>(maxY_, p->Y());
        }
    }
    BoundingBox(const std::vector<Vec2>& points) {
        if (points.empty()) {
            minX_ = minY_ = maxX_ = maxY_ = 0;
            return;
        }

        minX_ = maxX_ = geom::toGrid(points[0].x);
        minY_ = maxY_ = geom::toGrid(points[0].y);

        for (const auto& p : points) {
            minX_ = MathUtil::min<int64_t>(minX_, geom::toGrid(p.x));
            minY_ = MathUtil::min<int64_t>(minY_, geom::toGrid(p.y));
            maxX_ = MathUtil::max<int64_t>(maxX_, geom::toGrid(p.x));
            maxY_ = MathUtil::max<int64_t>(maxY_, geom::toGrid(p.y));
        }
    }

    // Static method to create a bounding box that encompasses two given bounding boxes
    static BoundingBox createBox(BoundingBox* a, BoundingBox* b) {
        return BoundingBox(
            MathUtil::min(a->minX(), b->minX()),
            MathUtil::min(a->minY(), b->minY()),
            MathUtil::max(a->maxX(), b->maxX()),
            MathUtil::max(a->maxY(), b->maxY()));
    }

    std::shared_ptr<geom::Point> getMin() {
        return std::make_shared<geom::Point>(geom::fromGrid(minX_), geom::fromGrid(minY_));
    }

    std::shared_ptr<geom::Point> getMax() {
        return std::make_shared<geom::Point>(geom::fromGrid(maxX_), geom::fromGrid(maxY_));
    }

    // Getters for the bounding box coordinates
    int64_t minX() {
        return minX_;
    }
    int64_t minY() {
        return minY_;
    }
    int64_t maxX() {
        return maxX_;
    }
    int64_t maxY() {
        return maxY_;
    }

    // Method to check if this bounding box intersects with another bounding box
    bool intersects(const BoundingBox& other) const {
        return !(maxX_ < other.minX_ || other.maxX_ < minX_ ||
                 maxY_ < other.minY_ || other.maxY_ < minY_);
    }

    // Method to reset and expand the bounding box to include a new point
    void reset(const Vec2& p) {
        minX_ = geom::toGrid(p.x);
        minY_ = geom::toGrid(p.y);
        maxX_ = minX_;
        maxY_ = minY_;
    }

    // Expand the bbox to include a new point
    void expand(const Vec2& p) {
        int64_t x = geom::toGrid(p.x);
        int64_t y = geom::toGrid(p.y);
        if (x < minX_) minX_ = x;
        if (y < minY_) minY_ = y;
        if (x > maxX_) maxX_ = x;
        if (y > maxY_) maxY_ = y;
    }
};

#endif // BOUNDING_BOX_H