#ifndef POINT_H
#define POINT_H

#include <vector>
#include <string>
#include <cmath> 
#include <cstdint>
#include <memory>
#include <unordered_set>
#include "../utils/geom_const.h"

// Simple 2D vector struct for convenience in calculations
struct Vec2 {
    double x, y;
};

/*
    Point class represents a point in 2D space with both double precision coordinates (x, y) and int64_t grid coordinates (X, Y).
    The grid coordinates are used for precise geometric calculations while the double precision coordinates are used for plotting and other operations.
    The class provides methods for basic geometric transformations such as shifting, rotating, and mirroring points.
    It also includes a method to calculate the distance to another point and to check if two points are considered equal within a certain tolerance.
*/
namespace geom {

    // Forward declaration of Line class to avoid circular dependency
    class Line;

    // Constants for grid scaling and tolerance are defined in geom_const.h
    inline int64_t toGrid(double x) {
        return static_cast<int64_t>(x * GRID);
    }
    inline double fromGrid(int64_t X) {
        return static_cast<double>(X) / GRID;
    }

    /*
        Point class to represent a point in 2D space. 
    */
    class Point {
    private:
        // Double precision coordinates for plotting and general use
        double x_;
        double y_;

        // Integer grid coordinates for precise geometric calculations
        int64_t X_;
        int64_t Y_;

    public:
        // Set of parent lines that this point belongs to, useful for intersection logic
        std::unordered_set<std::shared_ptr<Line>> parentLines;

    public:
        // Constructors
        Point() : x_(0.0), y_(0.0) { updateCoors(); }
        Point(double x, double y) : x_(x), y_(y) { updateCoors(); }

        // Operator overloads for basic arithmetic operations
        Point operator*(double s) const {
            return Point(x() * s, y() * s);
        }

        Point operator+(const Point& other) const {
            return Point(x() + other.x(), y() + other.y());
        }

        // Getters and setters for both double and int64_t coordinates
        double x() const { return x_; }
        double y() const { return y_; }

        double X() const { return X_; }
        double Y() const { return Y_; }

        void x(double x) { 
            x_ = x; 
            X_ = toGrid(x);
        }

        void y(double y) {
            y_ = y; 
            Y_ = toGrid(y);
        }

        void X(int64_t X) {
            X_ = X;
            x_ = fromGrid(X);
            
        }

        void Y(int64_t Y) {
            Y_ = Y;
            y_ = fromGrid(Y);
        }

        // Geometric transformations
        void shift(double dx, double dy) { 
            x(x_ += dx);
            y(y_ += dy); 
        }

        void rotate(const Point& center, double angle) {
            double s = std::sin(angle);
            double c = std::cos(angle);

            double x0 = x_ - center.x();
            double y0 = y_ - center.y();

            double x1 = x0 * c - y0 * s;
            double y1 = x0 * s + y0 * c;

            x(x1 + center.x());
            y(y1 + center.y());
        }

        void mirror(const Point& a, const Point& b) {
            double dx = b.x() - a.x();
            double dy = b.y() - a.y();
            double dx2_dy2 = dx * dx + dy * dy;
            if (dx2_dy2 == 0) { 
                // Zero division
                return;
            }

            double t = ((x_ - a.x()) * dx + (y_ - a.y()) * dy) / dx2_dy2;

            // Projected point Q
            double qx = a.x() + t * dx;
            double qy = a.y() + t * dy;

            // Mirror P across Q
            x(2 * qx - x_);
            y(2 * qy - y_);
        }

        // Method to calculate distance to another point
        double distanceTo(const std::shared_ptr<Point> other) const {
            double dx = x_ - other->x();
            double dy = y_ - other->y();
            return std::sqrt(dx * dx + dy * dy);
        }

        // Method to check if two points are considered equal within a certain tolerance
        bool isEqual(std::shared_ptr<Point> p) {
            return fabs(X() - p->X()) <= SAME_POINT && fabs(Y() - p->Y()) <= SAME_POINT;
        }

    private:
        // Helper method to update grid coordinates based on double precision coordinates
        void updateCoors() {
            x(x());
            y(y());
        }
    };
}
#endif // POINT_H