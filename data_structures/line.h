#ifndef LINE_H
#define LINE_H

#include "point.h"
#ifndef NOMINMAX
#define NOMINMAX
#endif
#undef Polygon
#include <stdexcept>
#include <algorithm>
#include <functional>
#include "../bezier_intersector/bezier_intersection_util.h"
#include "../utils/math_util.h"
#include "bounding_box.h"
#include <iostream>

/*
    Line class is an abstract base class representing a geometric line segment.
    It defines the common interface for all types of lines (e.g., segments, arcs, Bezier curves, B-splines).
    The class includes methods for getting the bounding box, plotting the line, and finding intersections with other lines.
    It also provides utility methods for sorting the endpoints and checking if a point lies on the segment.
*/
namespace geom {

    // Forward declarations of other line types to avoid circular dependencies
    class Segment;
    class ArcSegment;
    class Bezier;
    class BSpline;

    // Intersection methods for different combinations of line types
    std::vector<std::shared_ptr<Point>> intersectSegmentSegment(std::shared_ptr<Segment> seg1, std::shared_ptr<Segment> seg2);
    std::vector<std::shared_ptr<Point>> intersectSegmentArc(std::shared_ptr<Segment> seg, std::shared_ptr<ArcSegment> arc);
    std::vector<std::shared_ptr<Point>> intersectSegmentBezier(std::shared_ptr<Segment> seg, std::shared_ptr<Bezier> bez);
    std::vector<std::shared_ptr<Point>> intersectSegmentBSpline(std::shared_ptr<Segment> seg, std::shared_ptr<BSpline> bspl);

    std::vector<std::shared_ptr<Point>> intersectArcArc(std::shared_ptr<ArcSegment> arc1, std::shared_ptr<ArcSegment> arc2);
    std::vector<std::shared_ptr<Point>> intersectArcBezier(std::shared_ptr<ArcSegment> arc, std::shared_ptr<Bezier> bez);
    std::vector<std::shared_ptr<Point>> intersectArcBSpline(std::shared_ptr<ArcSegment> arc, std::shared_ptr<BSpline> bspl);

    std::vector<std::shared_ptr<Point>> intersectBezierBezier(std::shared_ptr<Bezier> bez1, std::shared_ptr<Bezier> bez2);
    std::vector<std::shared_ptr<Point>> intersectBezierBSpline(std::shared_ptr<Bezier> bez, std::shared_ptr<BSpline> bspl); 

    std::vector<std::shared_ptr<Point>> intersectBSplineBSpline(std::shared_ptr<BSpline> bspl1, std::shared_ptr<BSpline> bspl2);

    /*
        Abstract Line class definition.
        The intersect method is a virtual method that will be overridden by derived classes to implement specific 
        intersection logic based on the types of lines involved.
        It takes another line as input and returns a vector of points where the two lines intersect.
    */
    class Line : public std::enable_shared_from_this<Line> {
    protected:
        // Start and end points of the line segment
        std::shared_ptr<Point> s_;
        std::shared_ptr<Point> e_;

    public:
        // Virtual destructor to ensure proper cleanup of derived classes
        virtual ~Line() = default;
        // Pure virtual methods to be implemented by derived classes
        virtual BoundingBox getBoundingBox() const = 0;
        virtual std::vector<std::shared_ptr<Point>> intersect(std::shared_ptr<Line>& other) = 0;

        // Getters and setters for start and end points
        void s(Point& s) {
            s_ = std::make_shared<Point>(s);
        }
        void e(Point& e) {
            e_ = std::make_shared<Point>(e);
        }

        std::shared_ptr<Point> s() const { return s_; }
        std::shared_ptr<Point> e() const { return e_; }

        void s(const std::shared_ptr<Point>& p)  { s_ = p; }
        void e(const std::shared_ptr<Point>& p)  { e_ = p; }

        // Method to get the middle point of the line segment
        std::vector<std::shared_ptr<Point>> getMiddleLeftRight() {
            std::vector<std::shared_ptr<Point>> res(2);

            double dx = e()->x() - s()->x();
            double dy = e()->y() - s()->y();

            double dist = std::sqrt(dx*dx + dy*dy);

            // Normalize
            dx/=dist;
            dy/=dist;

            // Scale
            dx*=MIDDLE_POINT_SHIFT;
            dy*=MIDDLE_POINT_SHIFT;

            res[0] = getMiddle();
            res[1] = getMiddle();

            // Shift with normal vector
            // [-dy, dx] and [dy, -dx]
            res[0]->x(res[0]->x() - dy);
            res[0]->y(res[0]->y() + dx);

            res[1]->x(res[1]->x() + dy);
            res[1]->y(res[1]->y() - dx);

            return res;
        }

        // Method to sort the endpoints of the line segment based on their coordinates
        bool generalSort() {
            int64_t dx = fabs(s()->X() - e()->X());
            int64_t dy = fabs(s()->Y() - e()->Y());
            bool hor;

            if (dx >= dy) {
                // Line is more horizontal → sort by X
                hor = true;
                sortByX();
            } else {
                // Line is more vertical → sort by Y
                hor = false;
                sortByY();
            }

            return hor;
        }

        void sortByX() {
            if (e_->X() < s_->X()) {
                flip();
            }
        }

        void sortByY() {
            if (e_->Y() < s_->Y()) {
                flip();
            }
        }

        // Method to flip the start and end points of the line segment
        void flip() {
            std::swap(s_, e_);
        }

        // Getter methods for BoundingBox calculations
        int64_t getMinX() {
            return MathUtil::min(s_->X(), e_->X());
        }

        int64_t getMinY() {
            return MathUtil::min(s_->Y(), e_->Y());
        }

        int64_t getMaxX() {
            return MathUtil::max(s_->X(), e_->X());
        }

        int64_t getMaxY() {
            return MathUtil::max(s_->Y(), e_->Y());
        }

        // Method to check if a point lies on the line segment
        bool isPointOnSegment(const std::shared_ptr<Point>& p) {
            std::shared_ptr<Point> a = s();
            std::shared_ptr<Point> b = e();

            // Collinear check using cross product
            double cross = (b->y() - a->y()) * (p->x() - a->x()) - (b->x() - a->x()) * (p->y() - a->y());
            if (std::fabs(cross) > EPS) return false;

            // Check if p is within bounding box of a and b
            int64_t minX = getMinX(), maxX = getMaxX();
            int64_t minY = getMinY(), maxY = getMaxY();

            if (p->X() + EPS_INT64 < minX || p->X() - EPS_INT64 > maxX) return false;
            if (p->Y() + EPS_INT64 < minY || p->Y() - EPS_INT64 > maxY) return false;

            return true;
        }

        // Method to get the middle point of the line segment
        std::shared_ptr<Point> getMiddle() {
            double x = (s()->x() + e()->x()) / 2;
            double y = (s()->y() + e()->y()) / 2;

            std::shared_ptr<Point> p = std::make_shared<Point>(x,y);
            p->parentLines.insert(self());
            return p;
        }

    protected:
        // Helper method to get a shared pointer to the current instance
        std::shared_ptr<Line> self() {
            return shared_from_this();
        }
    };

    /*
        Segment is the simplest Line definition, containing only start and end point
    */
    class Segment : public Line {
    public:
        // Constructors
        Segment() {
            s_ = std::make_shared<Point>();
            e_ = std::make_shared<Point>();
        }
        Segment(std::shared_ptr<Point> s, std::shared_ptr<Point> e) {
            s_ = std::move(s);
            e_ = std::move(e);
        }

        // Method to get the bounding box of the line segment
        BoundingBox getBoundingBox() const override {
            return BoundingBox(
                MathUtil::min(s_->X(), e_->X()),
                MathUtil::min(s_->Y(), e_->Y()),
                MathUtil::max(s_->X(), e_->X()),
                MathUtil::max(s_->Y(), e_->Y())
            );
        }

        // Method to find the intersection points between this line segment and another line
        std::vector<std::shared_ptr<Point>> intersect(std::shared_ptr<Line>& other) override {
            if (std::shared_ptr<Segment> s = std::dynamic_pointer_cast<Segment>(other)) {
                return intersectSegmentSegment(self(), s);
            } else if (std::shared_ptr<ArcSegment> a = std::dynamic_pointer_cast<ArcSegment>(other)) {
                return intersectSegmentArc(self(), a);
            } else if (std::shared_ptr<Bezier> b = std::dynamic_pointer_cast<Bezier>(other)) {
                return intersectSegmentBezier(self(), b);
            } else if (std::shared_ptr<BSpline> bs = std::dynamic_pointer_cast<BSpline>(other)) {
                return intersectSegmentBSpline(self(), bs);
            } else {
                throw std::runtime_error("Unknown line type for intersection");
            }
        }

        // Method to convert this line segment to a Bezier curve (which is just a straight line in this case)
        std::shared_ptr<Bezier> toBezier() {
            return std::make_shared<Bezier>(std::vector<std::shared_ptr<Point>>{s(), e()});
        }

    private:
        // Helper method to get a shared pointer to the current instance
        std::shared_ptr<Segment> self() {
            return std::dynamic_pointer_cast<Segment>(shared_from_this());
        }
    };

    /*
        ArcSegment represents a circular arc defined by a center point, a start point, and an end point.
        It includes methods to calculate the angle to the center, check if a point lies on the arc, and 
        convert the arc to Bezier curves for rendering.
    */
    class ArcSegment : public Line {
    private:
        // Center point of the arc
        std::shared_ptr<Point> center_;

    public:
        // Constructors
        ArcSegment(std::shared_ptr<Point> center, std::shared_ptr<Point> s, std::shared_ptr<Point> e) {
            center_ = std::move(center);
            s_ = std::move(s);
            e_ = std::move(e);
        }

        // Getters and setters for the center point
        void center(Point& center) {
            center_ = std::make_shared<Point>(center);
        }

        std::shared_ptr<Point> center() const { 
            return center_; 
        }

        void center(const std::shared_ptr<Point>& center)  {
            center_ = center; 
        }

        // Method to calculate the angle from the center to a given point
        double angleToCenter(std::shared_ptr<Point> p) {
            return std::atan2(p->y() - center()->y(), p->x() - center()->x());
        }

        // Methods to get the start and end angles of the arc relative to the center
        double thetaStart() {
            return angleToCenter(s());
        }

        // Method to get the end angle of the arc relative to the center
        double thetaEnd() {
            return angleToCenter(e());
        }

        // Method to calculate the radius of the arc based on the distance from the center to the start or end point
        double radius() const {
            return std::hypot(s()->x() - center()->x(), s()->y() - center()->y());
        }

        // Helper method to normalize an angle to the range [0, 2*PI)
        double normalizeAngle(double a) {
            // Normalize angle to [0, 2*PI)
            while (a < 0) a += 2*M_PI;
            while (a >= 2*M_PI) a -= 2*M_PI;
            return a;
        }

        // Method to check if a point lies on the arc segment
        bool isPointOnArc(std::shared_ptr<Point> p) {
            // Check if point is on arc segment
            double aStart = normalizeAngle(thetaStart());
            double aEnd = normalizeAngle(thetaEnd());
            double aP = normalizeAngle(angleToCenter(p));

            // Check angle in arc range (counter-clockwise)
            if (aStart < aEnd) return aP >= aStart && aP <= aEnd;
            else return aP >= aStart || aP <= aEnd; // wrap around 2*PI
        }

        // Method to get the bounding box of the arc segment 
        // (currently a placeholder that uses the radius to create a bounding box around the center)
        BoundingBox getBoundingBox() const override {
            // TODO: real arc bounding box (angle-aware)
            // placeholder: bounding box of center + radius
            int64_t R = toGrid(MathUtil::max(
                center()->distanceTo(s_),
                center()->distanceTo(e_)
            ));

            return BoundingBox(
                center_->X() - R,
                center_->Y() - R,
                center_->X() + R,
                center_->Y() + R
            );
        }

        // Method to find the intersection points between this arc segment and another line
        std::vector<std::shared_ptr<Point>> intersect(std::shared_ptr<Line>& other) override {
            if (std::shared_ptr<Segment> s = std::dynamic_pointer_cast<Segment>(other)) {
                return intersectSegmentArc(s, self());
            } else if (std::shared_ptr<ArcSegment> a = std::dynamic_pointer_cast<ArcSegment>(other)) {
                return intersectArcArc(self(), a);
            } else if (std::shared_ptr<Bezier> b = std::dynamic_pointer_cast<Bezier>(other)) {
                return intersectArcBezier(self(), b);
            } else if (std::shared_ptr<BSpline> bs = std::dynamic_pointer_cast<BSpline>(other)) {
                return intersectArcBSpline(self(), bs);
            } else {
                throw std::runtime_error("Unknown line type for intersection");
            }
        }

        // Method to convert this arc segment to a series of Bezier curves for rendering
        std::vector<std::shared_ptr<Bezier>> toBezier() {
            std::vector<std::shared_ptr<Bezier>> beziers;
            // Convert arc to cubic bezier approximation
            double r = radius();
            double theta1 = thetaStart();
            double theta2 = thetaEnd();

            // Ensure theta2 > theta1 for proper segment length calculation
            if (theta2 < theta1) {
                theta2 += 2 * M_PI;
            }

            double deltaTheta = theta2 - theta1;

            // Number of segments (max 90° each)
            int segments = static_cast<int>(std::ceil(std::abs(deltaTheta) / (M_PI / 2.0)));
            double step = deltaTheta / segments;

            for (int i = 0; i < segments; ++i) {
                double a0 = theta1 + i * step;
                double a1 = a0 + step;

                // Calculate control points
                double k = 4.0 / 3.0 * std::tan(step / 4.0);

                std::shared_ptr<Point> p0 = std::make_shared<Point>(
                    center()->x() + radius() * std::cos(a0),
                    center()->y() + radius() * std::sin(a0)
                );

                std::shared_ptr<Point> p3 = std::make_shared<Point>(
                    center()->x() + radius() * std::cos(a1),
                    center()->y() + radius() * std::sin(a1)
                );

                std::shared_ptr<Point> p1 = std::make_shared<Point>(
                    center()->x() + radius() * std::cos(a0) - k * radius() * std::sin(a0),
                    center()->y() + radius() * std::sin(a0) + k * radius() * std::cos(a0)
                );

                std::shared_ptr<Point> p2 = std::make_shared<Point>(
                    center()->x() + radius() * std::cos(a1) + k * radius() * std::sin(a1),
                    center()->y() + radius() * std::sin(a1) - k * radius() * std::cos(a1)
                );

                beziers.push_back(std::make_shared<Bezier>(std::vector<std::shared_ptr<Point>>{p0, p1, p2, p3}));
            }
            
            return beziers;
        }

    private:
        // Helper method to get a shared pointer to the current instance
        std::shared_ptr<ArcSegment> self() {
            return std::dynamic_pointer_cast<ArcSegment>(shared_from_this());
        }
    };

    /*
        Bezier class represents a cubic Bezier curve defined by a set of control points.
        It includes methods to evaluate the curve at a given parameter t, calculate the bounding box, 
        and find intersections with other lines.
    */
    class Bezier : public Line {
    private:
        // Control points defining the Bezier curve (including start and end points)
        std::vector<std::shared_ptr<Point>> control_;
    public:
        // Constructors
        Bezier(std::vector<std::shared_ptr<Point>> control) {
            if (control.size() < 2) {
                throw std::invalid_argument("Bezier requires at least 2 control points");
            }

            control_ = std::move(control);
            s_ = control_.front();
            e_ = control_.back();
        }

        Bezier(std::vector<double> x, std::vector<double> y) : Bezier([&]() {
            std::vector<std::shared_ptr<Point>> control(x.size());
            for (size_t i = 0; i < x.size(); ++i) {
                control[i] = std::make_shared<Point>(x[i], y[i]);
            }
            return control; }()) {}

        // Getter for control points
        std::vector<std::shared_ptr<Point>> control() const { return control_; }

        // Method to get the bounding box of the Bezier curve by evaluating the curve at multiple points and finding the min/max coordinates
        BoundingBox getBoundingBox() const override {
            double minX = MathUtil::min(s()->X(), e()->X());
            double minY = MathUtil::min(s()->Y(), e()->Y());
            double maxX = MathUtil::max(s()->X(), e()->X());
            double maxY = MathUtil::max(s()->Y(), e()->Y());

            for (const auto& p : control_) {
                minX = MathUtil::min(minX, p->X());
                minY = MathUtil::min(minY, p->Y());
                maxX = MathUtil::max(maxX, p->X());
                maxY = MathUtil::max(maxY, p->Y());
            }

            return BoundingBox((int64_t) minX, (int64_t) minY, (int64_t) maxX, (int64_t) maxY);
        }

        // Method to find the intersection points between this Bezier curve and another line by
        // checking the type of the other line and calling the appropriate intersection function
        std::vector<std::shared_ptr<Point>> intersect(std::shared_ptr<Line>& other) override {
            if (std::shared_ptr<Segment> s = std::dynamic_pointer_cast<Segment>(other)) {
                return intersectSegmentBezier(s, self());
            } else if (std::shared_ptr<ArcSegment> a = std::dynamic_pointer_cast<ArcSegment>(other)) {
                return intersectArcBezier(a, self());
            } else if (std::shared_ptr<Bezier> b = std::dynamic_pointer_cast<Bezier>(other)) {
                return intersectBezierBezier(self(), b);
            } else if (std::shared_ptr<BSpline> bs = std::dynamic_pointer_cast<BSpline>(other)) {
                return intersectBezierBSpline(self(), bs);
            } else {
                throw std::runtime_error("Unknown line type for intersection");
            }
        }
    private: 
        // Method to evaluate the Bezier curve at a given parameter t using De Casteljau's algorithm
        std::shared_ptr<Point> evaluate(double t) const {
            std::vector<std::shared_ptr<Point>> temp;
            temp.reserve(control_.size());

            for (const auto& p : control_) {
                temp.push_back(std::make_shared<Point>(p->x(), p->y()));
            }

            for (size_t k = temp.size() - 1; k > 0; --k) {
                for (size_t i = 0; i < k; ++i) {
                    double x = (1 - t) * temp[i]->x() + t * temp[i + 1]->x();
                    double y = (1 - t) * temp[i]->y() + t * temp[i + 1]->y();
                    temp[i]->x(x);
                    temp[i]->y(y); 
                }
            }
            return temp[0];
        }

    private:
        // Helper method to get a shared pointer to the current instance
        std::shared_ptr<Bezier> self() {
            return std::dynamic_pointer_cast<Bezier>(shared_from_this());
        }
    };

    /*
        BSpline class represents a B-spline curve defined by a set of control points, a degree, and a knot vector.
        It includes methods to evaluate the curve at a given parameter t using De Boor's algorithm, calculate the bounding box, 
        and find intersections with other lines. It also has a method to convert the B-spline to a series of Bezier curves for rendering.
    */
    class BSpline : public Line {
    private:
        // Control points defining the B-spline curve
        std::vector<std::shared_ptr<Point>> control_;
        // Knot vector defining the parameterization of the B-spline curve
        std::vector<double> knots_;
        // Degree of the B-spline curve
        int degree_;
    public:
        // Constructors
        BSpline(std::vector<std::shared_ptr<Point>> control,
                int degree)
            : control_(std::move(control)),
            degree_(degree) {
            buildUniformClampedKnots();

            s_ = control_.front();
            e_ = control_.back();
        }

        BSpline(std::vector<std::shared_ptr<Point>> control,
                std::vector<double> knots,
                int degree)
            : control_(std::move(control)),
            knots_(std::move(knots)),
            degree_(degree) {
            validate();

            s_ = control_.front();
            e_ = control_.back();
        }

        // Method to evaluate the B-spline curve at a given parameter t using De Boor's algorithm
        std::shared_ptr<Point> evaluate(double t) const {
            int k = findKnotSpan(t);
            return deBoor(k, t);
        }

        // Method to get the bounding box of the B-spline curve by evaluating the curve at multiple points and finding the min/max coordinates
        BoundingBox getBoundingBox() const override {
            double minX = MathUtil::min(s()->X(), e()->X());
            double minY = MathUtil::min(s()->Y(), e()->Y());
            double maxX = MathUtil::max(s()->X(), e()->X());
            double maxY = MathUtil::max(s()->Y(), e()->Y());

            for (const auto& p : control_) {
                minX = MathUtil::min(minX, p->X());
                minY = MathUtil::min(minY, p->Y());
                maxX = MathUtil::max(maxX, p->X());
                maxY = MathUtil::max(maxY, p->Y());
            }

            return BoundingBox((int64_t) minX, (int64_t) minY, (int64_t) maxX, (int64_t) maxY);
        }

        // Method to find the intersection points between this B-spline curve and another line by
        std::vector<std::shared_ptr<Point>> intersect(std::shared_ptr<Line>& other) override {
            if (std::shared_ptr<Segment> s = std::dynamic_pointer_cast<Segment>(other)) {
                return intersectSegmentBSpline(s, self());
            } else if (std::shared_ptr<ArcSegment> a = std::dynamic_pointer_cast<ArcSegment>(other)) {
                return intersectArcBSpline(a, self());
            } else if (std::shared_ptr<Bezier> b = std::dynamic_pointer_cast<Bezier>(other)) {
                return intersectBezierBSpline(b, self());
            } else if (std::shared_ptr<BSpline> bs = std::dynamic_pointer_cast<BSpline>(other)) {
                return intersectBSplineBSpline(self(), bs);
            } else {
                throw std::runtime_error("Unknown line type for intersection");
            }
        }

        /*
        This implementation:
            Exact conversion (no loss)
            Resulting Bézier segments are C² continuous

        Limitations:
            Assumes uniform cubic B-spline
            For non-uniform knots, you must insert knots first
            For NURBS, weights must be handled (rational Bézier)
        */
        // Method to convert this B-spline curve to a series of Bezier curves for rendering by calculating
        // the appropriate control points for each segment
        std::vector<std::shared_ptr<Bezier>> toBezier(){
            std::vector<std::shared_ptr<Bezier>> beziers;

            if (control_.size() < 4)
                return beziers;

            for (size_t i = 0; i + 3 < control_.size(); ++i) {
                std::shared_ptr<Point> P0 = control_[i];
                std::shared_ptr<Point> P1 = control_[i + 1];
                std::shared_ptr<Point> P2 = control_[i + 2];
                std::shared_ptr<Point> P3 = control_[i + 3];

                std::vector<std::shared_ptr<Point>> control(4);
                control[0] = std::make_shared<Point>(
                    (P0->x() + 4.0 * P1->x() + P2->x()) / 6.0,
                    (P0->y() + 4.0 * P1->y() + P2->y()) / 6.0
                );
                control[1] = std::make_shared<Point>(
                    (4.0 * P1->x() + 2.0 * P2->x()) / 6.0,
                    (4.0 * P1->y() + 2.0 * P2->y()) / 6.0
                );
                control[2] = std::make_shared<Point>(
                    (2.0 * P1->x() + 4.0 * P2->x()) / 6.0,
                    (2.0 * P1->y() + 4.0 * P2->y()) / 6.0
                );
                control[3] = std::make_shared<Point>(
                    (P1->x() + 4.0 * P2->x() + P3->x()) / 6.0,
                    (P1->y() + 4.0 * P2->y() + P3->y()) / 6.0
                );

                beziers.push_back(std::make_shared<Bezier>(control));
            }

            return beziers;
        }

        std::vector<std::shared_ptr<Point>> control() const { return control_; }

    private:
        // Method to build a uniform clamped knot vector based on the number of control points and the degree of the B-spline
        void buildUniformClampedKnots() {
            int n = static_cast<int>(control_.size()) - 1;
            int m = n + degree_ + 1;

            knots_.resize(m + 1);

            for (int i = 0; i <= m; ++i) {
                if (i <= degree_)
                    knots_[i] = 0.0;
                else if (i >= m - degree_)
                    knots_[i] = 1.0;
                else
                    knots_[i] = double(i - degree_) / double(m - 2 * degree_);
            }
        }

        // Method to validate the B-spline parameters, ensuring that there are enough control points
        // for the specified degree and that the knot vector has the correct size
        void validate() const {
            if (control_.size() < static_cast<size_t>(degree_ + 1))
                throw std::runtime_error("Not enough control points for degree");

            if (knots_.size() != control_.size() + degree_ + 1)
                throw std::runtime_error("Invalid knot vector size");
        }

        // Method to find the knot span index for a given parameter t,
        // which is used in De Boor's algorithm to determine which control points influence the curve at that parameter value
        int findKnotSpan(double t) const {
            int n = static_cast<int>(control_.size()) - 1;

            if (t >= knots_[n + 1])
                return n;

            for (int i = degree_; i <= n; ++i) {
                if (t >= knots_[i] && t < knots_[i + 1])
                    return i;
            }

            return degree_;
        }

        // Method to perform De Boor's algorithm to evaluate the B-spline curve at a given parameter t,
        // which involves iteratively blending the control points based on the knot vector and the degree of the curve
        std::shared_ptr<Point> deBoor(int k, double t) const {
            std::vector<Point> d(degree_ + 1);

            for (int j = 0; j <= degree_; ++j)
                d[j] = *control_[k - degree_ + j];

            for (int r = 1; r <= degree_; ++r) {
                for (int j = degree_; j >= r; --j) {
                    int i = k - degree_ + j;
                    double alpha =
                        (t - knots_[i]) /
                        (knots_[i + degree_ - r + 1] - knots_[i]);

                    d[j] = d[j - 1] * (1.0 - alpha) + d[j] * alpha;
                }
            }

            return std::make_shared<Point>(d[degree_]);
        }

    private:
        // Helper method to get a shared pointer to the current instance
        std::shared_ptr<BSpline> self() {
            return std::dynamic_pointer_cast<BSpline>(shared_from_this());
        }
    };

    /*
        Implementation of the intersection logic for the sweep line algorithm. This is used in the boolean operations to find all intersection points between lines.
        The logic is based on the type of the line (Segment, ArcSegment, Bezier, BSpline) and the intersection functions defined above. 
        The intersection points are collected in a vector and returned. Additionally, if the lines are collinear or parallel, we check if the endpoints of one line lie on the other line and add those points to the intersection set as well.
    */

    // Intersection functions for different line types
    std::vector<std::shared_ptr<Point>> intersectSegmentSegment(std::shared_ptr<Segment> seg1, std::shared_ptr<Segment> seg2) {
        std::vector<std::shared_ptr<Point>> intersections;

        double x1 = seg1->s()->x(), y1 = seg1->s()->y();
        double x2 = seg1->e()->x(), y2 = seg1->e()->y();
        double x3 = seg2->s()->x(), y3 = seg2->s()->y();
        double x4 = seg2->e()->x(), y4 = seg2->e()->y();

        double dx1 = x2 - x1, dy1 = y2 - y1;
        double dx2 = x4 - x3, dy2 = y4 - y3;

        double denom = dx1 * dy2 - dy1 * dx2;
        auto insertIfOnSegment = [](const std::shared_ptr<Point>& p, const std::shared_ptr<Segment>& seg) {
            if (seg->isPointOnSegment(p)) p->parentLines.insert(seg);
        };

        if (fabs(denom) < EPS) {
            // Collinear or parallel case
            insertIfOnSegment(seg1->s(), seg2);
            insertIfOnSegment(seg1->e(), seg2);
            insertIfOnSegment(seg2->s(), seg1);
            insertIfOnSegment(seg2->e(), seg1);
            return intersections;
        }

        double t = ((x3 - x1) * dy2 - (y3 - y1) * dx2) / denom;
        double u = ((x3 - x1) * dy1 - (y3 - y1) * dx1) / denom;
        if (t < 0 || t > 1 || u < 0 || u > 1) return intersections;

        auto p = std::make_shared<Point>(x1 + t * dx1, y1 + t * dy1);
        p->parentLines.insert(seg1);
        p->parentLines.insert(seg2);

        intersections.push_back(p);
        return intersections;
    }

    std::vector<std::shared_ptr<Point>> intersectSegmentArc(std::shared_ptr<Segment> seg, std::shared_ptr<ArcSegment> arc) {
        std::vector<std::shared_ptr<Point>> intersections;

        // compute arc angles and radius
        double r = arc->radius();

        double theta_start = arc->thetaStart();
        double theta_end   = arc->thetaEnd();

        if (theta_start < 0) theta_start += 2*M_PI;
        if (theta_end   < 0) theta_end   += 2*M_PI;

        double dx = seg->e()->x() - seg->s()->x();
        double dy = seg->e()->y() - seg->s()->y();
        double fx = seg->s()->x() - arc->center()->x();
        double fy = seg->s()->y() - arc->center()->y();

        double A = dx*dx + dy*dy;
        double B = 2*(fx*dx + fy*dy);
        double C = fx*fx + fy*fy - r*r;

        double disc = B*B - 4*A*C;
        if (disc < 0) return intersections;

        double sqrtD = std::sqrt(disc);
        double ts[2] = { (-B + sqrtD)/(2*A), (-B - sqrtD)/(2*A) };

        for (double t : ts) {
            if (t < 0 || t > 1) continue;

            double ix = seg->s()->x() + t*dx;
            double iy = seg->s()->y() + t*dy;

            double angle = std::atan2(iy - arc->center()->y(), ix - arc->center()->x());
            if (angle < 0) angle += 2*M_PI;

            bool insideArc;
            if (theta_start <= theta_end) {
                insideArc = (angle >= theta_start && angle <= theta_end);
            } else {
                insideArc = (angle >= theta_start || angle <= theta_end);
            }

            if (insideArc) intersections.push_back(std::make_shared<Point>(ix, iy));
        }

        return intersections;
    }
    
    std::vector<std::shared_ptr<Point>> intersectSegmentBezier(std::shared_ptr<Segment> seg, std::shared_ptr<Bezier> bez) {
        return intersectBezierBezier(seg->toBezier(), bez);
    }
    
    std::vector<std::shared_ptr<Point>> intersectSegmentBSpline(std::shared_ptr<Segment> seg, std::shared_ptr<BSpline> bspl) {
        std::vector<std::shared_ptr<Bezier>> bezier = bspl->toBezier();
        std::vector<std::shared_ptr<Point>> intersections;

        for (const auto& b : bezier) {
            std::vector<std::shared_ptr<Point>> inters = intersectSegmentBezier(seg, b);
            intersections.insert(intersections.end(), inters.begin(), inters.end());
        }
        return intersections;
    }
    
    std::vector<std::shared_ptr<Point>> intersectArcArc(std::shared_ptr<ArcSegment> arc1, std::shared_ptr<ArcSegment> arc2) {
        std::vector<std::shared_ptr<Point>> intersections;

        double dx = arc2->center()->x() - arc1->center()->x();
        double dy = arc2->center()->y() - arc1->center()->y();
        double d = std::hypot(dx, dy);

        if (d < 1e-9) {
            if (std::abs(arc1->radius() - arc2->radius()) < 1e-9) {
                // Coincident arcs — return overlapping endpoints
                if (arc1->isPointOnArc(arc2->s())) intersections.push_back(arc2->s());
                if (arc1->isPointOnArc(arc2->e())) intersections.push_back(arc2->e());
                if (arc2->isPointOnArc(arc1->s())) intersections.push_back(arc1->s());
                if (arc2->isPointOnArc(arc1->e())) intersections.push_back(arc1->e());
            }
            return intersections;
        }

        // No intersection: circles too far apart or one inside another
        if (d > arc1->radius() + arc2->radius() || d < std::abs(arc1->radius() - arc2->radius())) return intersections;

        // Circle-circle intersection formulas
        double a = (arc1->radius()*arc1->radius() - arc2->radius()*arc2->radius() + d*d) / (2*d);
        double h = std::sqrt(arc1->radius()*arc1->radius() - a*a);

        double cx2 = arc1->center()->x() + a * dx / d;
        double cy2 = arc1->center()->y() + a * dy / d;

        // Two intersection points
        std::shared_ptr<Point> p1 = std::make_shared<Point>(cx2 + h * dy / d, cy2 - h * dx / d);
        std::shared_ptr<Point> p2 = std::make_shared<Point>(cx2 - h * dy / d, cy2 + h * dx / d);

        if (arc1->isPointOnArc(p1) && arc2->isPointOnArc(p1)) intersections.push_back(p1);
        if (arc1->isPointOnArc(p2) && arc2->isPointOnArc(p2)) intersections.push_back(p2);

        return intersections;
    }
    
    std::vector<std::shared_ptr<Point>> intersectArcBezier(std::shared_ptr<ArcSegment> arc, std::shared_ptr<Bezier> bez) {
        std::vector<std::shared_ptr<Bezier>> bezier = arc->toBezier();
        std::vector<std::shared_ptr<Point>> intersections;

        for (const auto& b : bezier) {
            std::vector<std::shared_ptr<Point>> inters = intersectBezierBezier(b, bez);
            intersections.insert(intersections.end(), inters.begin(), inters.end());
        }

        return intersections;
    }
    
    std::vector<std::shared_ptr<Point>> intersectArcBSpline(std::shared_ptr<ArcSegment> arc, std::shared_ptr<BSpline> bspl) {
        std::vector<std::shared_ptr<Bezier>> bezier1 = arc->toBezier();
        std::vector<std::shared_ptr<Bezier>> bezier2 = bspl->toBezier();
        std::vector<std::shared_ptr<Point>> intersections;

        for (const auto& b1 : bezier1) {
            for (const auto& b2 : bezier2) {
                std::vector<std::shared_ptr<Point>> inters = intersectBezierBezier(b1, b2);
                intersections.insert(intersections.end(), inters.begin(), inters.end());
            }
        }
        return intersections;
    }
    
    std::vector<std::shared_ptr<Point>> intersectBezierBSpline(std::shared_ptr<Bezier> bez, std::shared_ptr<BSpline> bspl) {
        std::vector<std::shared_ptr<Bezier>> bezier = bspl->toBezier();
        std::vector<std::shared_ptr<Point>> intersections;

        for (const auto& b : bezier) {
            std::vector<std::shared_ptr<Point>> inters = intersectBezierBezier(bez, b);
            intersections.insert(intersections.end(), inters.begin(), inters.end());
        }
        return intersections;
    }
    
    std::vector<std::shared_ptr<Point>> intersectBSplineBSpline(std::shared_ptr<BSpline> bspl1, std::shared_ptr<BSpline> bspl2){
        std::vector<std::shared_ptr<Bezier>> bezier1 = bspl1->toBezier();
        std::vector<std::shared_ptr<Bezier>> bezier2 = bspl2->toBezier();
        std::vector<std::shared_ptr<Point>> intersections;

        for (const auto& b1 : bezier1) {
            for (const auto& b2 : bezier2) {
                std::vector<std::shared_ptr<Point>> inters = intersectBezierBezier(b1, b2);
                intersections.insert(intersections.end(), inters.begin(), inters.end());
            }
        }
        return intersections;
    }

    std::vector<std::shared_ptr<Point>> intersectBezierBezier(std::shared_ptr<Bezier> bez1, std::shared_ptr<Bezier> bez2) {
        std::vector<std::shared_ptr<Point>> intersections;
        std::vector<BEZIER_INTERSECTION_UTIL::IntersectionCanditate> candidates;

        std::vector<Vec2> control1 = BEZIER_INTERSECTION_UTIL::toVec2(bez1->control());
        std::vector<Vec2> control2 = BEZIER_INTERSECTION_UTIL::toVec2(bez2->control());

        BEZIER_INTERSECTION_UTIL::findIntersectionCandidates(control1, control2, candidates);

        for (auto& c : candidates) {
            GaussNewtonSolver solver(control1, control2, c);
            auto pt = solver.solve();
            if (pt) {
                intersections.push_back(pt);
            }
        }

        return intersections;
    }

    std::vector<std::shared_ptr<Point>> intersectBezierBezier(std::vector<std::shared_ptr<Bezier>> bez1, std::shared_ptr<Bezier> bez2) {
        std::vector<std::shared_ptr<Point>> intersections;
        for (const auto& b1 : bez1) {
            std::vector<std::shared_ptr<Point>> inters = intersectBezierBezier(b1, bez2);
            intersections.insert(intersections.end(), inters.begin(), inters.end());
        }
        return intersections;
    }

    std::vector<std::shared_ptr<Point>> intersectBezierBezier(std::vector<std::shared_ptr<Bezier>> bez1, std::vector<std::shared_ptr<Bezier>> bez2) {
        std::vector<std::shared_ptr<Point>> intersections;
        for (const auto& b1 : bez1) {
            for (const auto& b2 : bez2) {
                std::vector<std::shared_ptr<Point>> inters = intersectBezierBezier(b1, b2);
                intersections.insert(intersections.end(), inters.begin(), inters.end());
            }
        }
        return intersections;
    }

}
#endif // LINE_H