#ifndef GEO_H
#define GEO_H

#include "line.h"

/*
    * This file defines the Path and Polygon classes, which are fundamental geometric structures.
    * The Path class represents a sequence of points/lines that form a closed shape, while the Polygon class represents a collection of paths (contours) that together define a polygonal shape. 
    * Both classes include methods for geometric transformations such as shifting, rotating, and mirroring, as well as methods for determining lines and calculating areas.
*/
namespace geom {
    class Path {
    private:
        // A vector of shared pointers to Point objects that define the vertices of the path.
        std::vector<std::shared_ptr<Point>> points_;
        /*
        Only necessary to determine this during boolean operation
        Idea is that we start from a line, which is the start and check if any of the 
        other polygons' lines intersects. 
        If yes, then collect points and write them in contours.
        Later create the new polygon definitions.
        */ 
        // A vector of shared pointers to Line objects that represent the edges of the path.
        std::vector<std::shared_ptr<Line>> lines_;

        // A boolean flag to indicate whether the path is a hole (negative area) or not, and a double to store the area of the path.
        bool isHole_;
        double area_;
    public:
        // Constructors
        Path(std::vector<std::shared_ptr<Point>> p) : points_(std::move(p)) {
            double a = signed_area();
            area_ = std::abs(a);
            isHole_ = a < 0;
        }

        Path(std::vector<std::shared_ptr<Line>> l) : lines_(std::move(l)) {
            determinePoints();
            double a = signed_area();
            area_ = std::abs(a);
            isHole_ = a < 0;
        }

        // Destructor
        ~Path() = default;

        // Getter and setter for points
        void points(std::vector<std::shared_ptr<Point>> p) { points_ = std::move(p); }

        std::vector<std::shared_ptr<Point>>& points() {
            return points_;
        }

        // Getter for hole status
        bool isHole() { return isHole_;  }

        // Geometrical operations
        void shift(double dx, double dy) { 
            for (const auto& p : points_) {
                p->shift(dx, dy);
            }
        }

        void rotate(const Point& center, double angle) {
            for (const auto& p : points_) {
                p->rotate(center, angle);
            }
        }

        void mirror(const Point& a, const Point& b) {
            for (const auto& p : points_) {
                p->mirror(a, b);
            }
        }

        void mirror(const Line& axis) {
            mirror(*axis.s(), *axis.e());
        }

        // Method to determine the sequence of points from the given lines, ensuring that they form a valid path.
        void determinePoints() {
            points_.clear();
            
            if (lines_.size() < 3) {
                return;
            }

            std::shared_ptr<Point> s = lines_[0]->s();
            std::shared_ptr<Point> e = lines_[0]->e();

            // Check which one is found in the next line
            std::shared_ptr<Point> sNext = lines_[1]->s();
            std::shared_ptr<Point> eNext = lines_[1]->e();

            if (s->isEqual(sNext)) {
                points_.push_back(e);
                points_.push_back(s);
                points_.push_back(eNext);

            } else if (s->isEqual(eNext)) {
                points_.push_back(e);
                points_.push_back(s);
                points_.push_back(sNext);

            } else if (e->isEqual(sNext)) {
                points_.push_back(s);
                points_.push_back(e);
                points_.push_back(eNext);

            } else if (e->isEqual(eNext)) {
                points_.push_back(s);
                points_.push_back(e);
                points_.push_back(sNext);

            } else {
                throw std::runtime_error("Invalid line sequence in path");
            }

            for (size_t i = 2; i < lines_.size(); ++i) {
                if (lines_[i]->s()->isEqual(e)) {
                    e = lines_[i]->e();
                    points_.push_back(e);
                } else if (lines_[i]->e()->isEqual(e)) {
                    e = lines_[i]->s();
                    points_.push_back(e);
                } else {
                    throw std::runtime_error("Invalid line sequence in path");
                }
            }
        }

        // Method to determine the sequence of lines from the given points, ensuring that they form a valid path.
        void determineLines() {
            if (!lines_.empty()) { 
                return;
            }
            for (size_t i = 0; i < points_.size(); ++i) {
                lines_.push_back(std::make_shared<Segment>(points_[i], points_[(i+1)%points_.size()]));
            }
        }

        // Getter for lines
        std::vector<std::shared_ptr<Line>>& getLines() {
            return lines_;
        }

        // Method to clear the lines vector. Might not need since we can't
        // generate the same path with different lines, but just in case.
        void deleteLines() {
            lines_.clear();
        }

    private:
        // Method to calculate the signed area of the path using the shoelace formula.
        // The sign of the area indicates whether the path is a hole (negative) or not (positive). 
        double signed_area() const {
            const size_t n = points_.size();
            if (n < 3) return 0.0; // Not a polygon

            double sum = 0.0;
            for (size_t i = 0; i < n; ++i) {
                const Point& p1 = *points_[i];
                const Point& p2 = *points_[(i + 1) % n];
                sum += p1.x() * p2.y() - p2.x() * p1.y();
            }
            return sum / 2.0;
        }
    };

    /*
        The Polygon class represents a polygonal shape defined by one or more contours (paths).
        Each contour is a closed sequence of points/lines that can represent either an outer boundary or a hole within the polygon. 
        The class includes methods for geometric transformations, determining lines from points, and calculating areas.
    */
    class Polygon {
    private:
        // A vector of shared pointers to Path objects that define the contours of the polygon. Each contour can represent either an outer boundary or a hole.
        std::vector<std::shared_ptr<Path>> contours_;
        // Not used yet, but can be useful for boolean operations to keep track of which components belong to which original polygons.
        std::vector<std::string> components_;
        // Not used yet, layer name can be useful for visualization or for keeping track of different types of polygons in a larger application.
        std::string layer_name_;

    public:
        // Constructors
        Polygon(std::vector<std::shared_ptr<Path>> p) : contours_(std::move(p)) {}
        Polygon(std::shared_ptr<Path> p) : contours_{std::move(p)} {}

        // Destructor
        ~Polygon() = default;

        // Getter and setter for contours
        void contours(std::vector<std::shared_ptr<Path>> p) { contours_ = std::move(p); }

        std::vector<std::shared_ptr<Path>>& contours() {
            return contours_;
        }

        // Geometrical operations
        void shift(double dx, double dy) {
            for (std::shared_ptr<Path> p : contours_) {
                p->shift(dx, dy);
            }
        }

        void rotate(const Point& center, double angle) {
            for (std::shared_ptr<Path> p : contours_) {
                p->rotate(center, angle);
            }
        }

        void mirror(const Point& a, const Point& b) {
            for (std::shared_ptr<Path> p : contours_) {
                p->mirror(a, b);
            }
        }

        void mirror(const Line& axis) {
            mirror(*axis.s().get(), *axis.e().get());
        }

        // Method to determine the sequence of lines for each path in the polygon, ensuring that they form valid contours.
        void determineLines() {
            for (std::shared_ptr<Path> p : contours_) {
                p->determineLines();
            }
        }

        // Method to collect all lines from all paths in the polygon into a single vector,
        // which can be useful for operations that need to consider the entire polygon as a collection of edges.
        std::vector<std::shared_ptr<Line>> getLines() {
            std::vector<std::shared_ptr<Line>> lines, tmp;

            for (std::shared_ptr<Path> p : contours_) {
                tmp = p->getLines();
                lines.reserve(lines.size() + tmp.size());
                lines.insert(lines.end(), tmp.begin(), tmp.end());
            }

            return lines;
        }

        // Method to print the points of each path in the polygon, which can be useful for debugging or visualization purposes.
        void print() {
            for (const auto& path : contours()) {
                std::cout << "Path: ";
                for (const auto& point : path->points()) {
                    std::cout << "P(" << point->x() << ", " << point->y() << ") -> ";
                }
                std::cout << "\b\b\b\b"  << "    "  << "\b\b\b\b"; 
                std::cout << std::endl;
            }
        }

        // Method to clear the lines from all paths in the polygon.
        void deleteLines() {
            for (std::shared_ptr<Path> p : contours_) {
                p->deleteLines();
            }
        }
    };
}
#endif // GEO_H