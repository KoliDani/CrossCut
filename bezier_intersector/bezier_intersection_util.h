#ifndef BEZIER_INTERSECTION_UTIL_H
#define BEZIER_INTERSECTION_UTIL_H

#include "../data_structures/point.h"
#include <vector>
#include "../data_structures/bounding_box.h"
#include "projected_gauss_newton.h"

/*
    This file contains utility functions and data structures for finding intersection candidates between two Bezier curves.
    It includes methods for subdividing Bezier curves, calculating pseudo-curvature, extracting segments of Bezier curves, and checking bounding box intersections.
    The main purpose of these utilities is to support the Projected Gauss-Newton method for accurately finding intersection points between Bezier curves.
    Algorithm: https://www.mdpi.com/2227-7390/12/9/1344
*/
namespace BEZIER_INTERSECTION_UTIL {

    // Constants for the Bezier curve intersection algorithm, including maximum curvature for subdivision, maximum recursion depth, and epsilon for range comparison.
    // Original KAPPA_MAX = 1.0 + 5 * 1e-6;
    // Original MAX_DEPTH = 20;
    constexpr double KAPPA_MAX = 1.001;
    const int MAX_DEPTH = 10;
    const double RANGE_EPS = 1e-6;

    // Convert everything to Vec2 for easier handling in Gauss-Newton solver
    inline std::vector<Vec2> toVec2(const std::vector<std::shared_ptr<geom::Point>>& ctrl) {
        std::vector<Vec2> vecs;
        vecs.reserve(ctrl.size());
        for (const auto& p : ctrl) {
            vecs.push_back({p->x(), p->y()});
        }
        return vecs;
    }

    // Data structure to represent a segment of a Bezier curve, 
    // including its control points and the parameter range on the original curve that corresponds to this segment.
    struct BezierSegment {
        std::vector<Vec2> control;
        BoundingBox bbox;
        // The parameter range on the original curve that corresponds to this segment, 
        // useful for mapping back to the original curves after subdivision.
        double t0, t1;
    };

    // Data structure to represent a candidate for intersection between two Bezier curves, consisting of parameter ranges for each curve.
    struct RangeDividedBezier {
        BezierSegment left;
        BezierSegment right;
    };

    // Bezier subdivision using De Casteljau's algorithm, which takes a set of control points and a parameter t,
    // and returns the control points of the left and right subdivided curves.
    RangeDividedBezier subdivideBezier(BezierSegment in, double t) {

        const size_t n = in.control.size();

        BezierSegment left, right;
        left.control.resize(n);
        right.control.resize(n);

        std::vector<Vec2> work(in.control.begin(), in.control.end());

        left.control[0] = work[0];
        right.control[n - 1] = work[n - 1];

        // Initialize bounding boxes with the first points
        left.bbox.reset(left.control[0]);
        right.bbox.reset(right.control[n - 1]);

        // De Casteljau triangle
        for (size_t level = 1; level < n; ++level) {
            for (size_t i = 0; i < n - level; ++i) {
                work[i].x = (1.0 - t) * work[i].x + t * work[i + 1].x;
                work[i].y = (1.0 - t) * work[i].y + t * work[i + 1].y;
            }

            left.control[level] = work[0];
            right.control[n - level - 1] = work[n - level - 1];

            // Incrementally expand bounding boxes using the new points
            left.bbox.expand(work[0]);
            right.bbox.expand(work[n - level - 1]);
        }

        // Parameter adjustment for the subdivided segments, which is important for mapping back to the original curves after subdivision.
        double tm = in.t0 + (in.t1 - in.t0) * t;
        left.t0 = in.t0;
        left.t1 = tm;
        right.t0 = tm;
        right.t1 = in.t1;

        return {left, right};
    }

    // Pseudo-curvature calculation for a Bezier curve segment, which is used to determine if the curve is flat enough for subdivision.
    // Returns: true if the curve is flat enough (kappa <= KAPPA_MAX), false otherwise
    inline bool isFlatEnough(const std::vector<Vec2>& cp) {
        const size_t n = cp.size();
        if (n < 2) return true;

        const Vec2& a = cp.front();
        const Vec2& b = cp.back();

        double dx = b.x - a.x;
        double dy = b.y - a.y;
        double chordSq = dx*dx + dy*dy;

        if (chordSq < geom::EPS * geom::EPS)
            return false;

        double polySq = 0.0;

        for (size_t i = 0; i < n - 1; ++i)
        {
            double ddx = cp[i+1].x - cp[i].x;
            double ddy = cp[i+1].y - cp[i].y;
            polySq += ddx*ddx + ddy*ddy;
        }

        // equivalent to polyLen/chordLen <= KAPPA_MAX
        return polySq <= chordSq * (KAPPA_MAX * KAPPA_MAX);
    }

    // Find the intersection candidates by recursively subdividing the Bezier curves based on their pseudo-curvature and bounding box intersection.
    inline void findIntersectionCandidatesReq(const BezierSegment& cp1, const BezierSegment& cp2,
                                std::vector<IntersectionCanditate>& candidates, int depth) {
        
        if (depth > MAX_DEPTH || !cp1.bbox.intersects(cp2.bbox)) {
            return;
        }

        bool ka = isFlatEnough(cp1.control);
        bool kb = isFlatEnough(cp2.control);

        if (ka && kb) {
            ParamRange range1{cp1.t0, cp1.t1};
            ParamRange range2{cp2.t0, cp2.t1};
            candidates.push_back({range1, range2});
            return;
        }

        if (ka) {
            RangeDividedBezier b1 = subdivideBezier(cp1, 0.5); 

            findIntersectionCandidatesReq(b1.left, cp2, candidates, depth + 1);
            findIntersectionCandidatesReq(b1.right, cp2, candidates, depth + 1);
            return;
        } else {
            RangeDividedBezier b2 = subdivideBezier(cp2, 0.5); 

            findIntersectionCandidatesReq(cp1, b2.left, candidates, depth + 1);
            findIntersectionCandidatesReq(cp1, b2.right, candidates, depth + 1);
        }
    }

    // Data structure to represent a node in the stack for the iterative 
    // subdivision process, which includes the two Bezier segments being compared and the current depth of recursion.
    struct Node {
        BezierSegment a;
        BezierSegment b;
        int depth;
    };

    // Find the intersection candidates by subdividing the Bezier curves based on their pseudo-curvature and bounding box intersection.
    inline void findIntersectionCandidatesStack(const BezierSegment& b1, const BezierSegment& b2,
                                std::vector<IntersectionCanditate>& candidates) {
        
        std::vector<Node> stack;
        stack.push_back({b1, b2, 0});

        while (!stack.empty()) {
            Node node = stack.back();
            stack.pop_back();

            const BezierSegment& cp1 = node.a;
            const BezierSegment& cp2 = node.b;
            int depth = node.depth;

            if (depth > MAX_DEPTH || !cp1.bbox.intersects(cp2.bbox))
                continue;

            bool ka = isFlatEnough(cp1.control);
            bool kb = isFlatEnough(cp2.control);

            if (ka && kb) {
                ParamRange range1{cp1.t0, cp1.t1};
                ParamRange range2{cp2.t0, cp2.t1};
                candidates.push_back({range1, range2});
                continue;
            }

            // Subdivide only the curve that is less flat
            if (ka) {
                RangeDividedBezier subdiv = subdivideBezier(cp1, 0.5);
                stack.push_back({subdiv.left, cp2, depth + 1});
                stack.push_back({subdiv.right, cp2, depth + 1});
            } else {
                RangeDividedBezier subdiv = subdivideBezier(cp2, 0.5);
                stack.push_back({cp1, subdiv.left, depth + 1});
                stack.push_back({cp1, subdiv.right, depth + 1});
            }
        }
    }

    // Main function to find intersection candidates between two Bezier curves,
    // which initializes the process by creating the initial segments and calling the subdivision function.
    inline void findIntersectionCandidates(const std::vector<Vec2>& bez1, const std::vector<Vec2>& bez2,
                                std::vector<IntersectionCanditate>& candidates) {

        BoundingBox bb1(bez1);
        BoundingBox bb2(bez2);
        if (!bb1.intersects(bb2)) {
            return;
        }

        bool ka = isFlatEnough(bez1);
        bool kb = isFlatEnough(bez2);

        if (ka && kb) {
            candidates.push_back({{0.0, 1.0}, {0.0, 1.0}});
            return;
        }

        if (ka) {
            RangeDividedBezier b1 = subdivideBezier(BezierSegment{bez1, bb1, 0.0, 1.0}, 0.5); 
            BezierSegment b2{bez2, bb2, 0.0, 1.0};

            findIntersectionCandidatesStack(b1.left, b2, candidates);
            findIntersectionCandidatesStack(b1.right, b2, candidates);
        } else {
            BezierSegment b1{bez1, bb1, 0.0, 1.0};
            RangeDividedBezier b2 = subdivideBezier(BezierSegment{bez2, bb2, 0.0, 1.0}, 0.5); 

            findIntersectionCandidatesStack(b1, b2.left, candidates);
            findIntersectionCandidatesStack(b1, b2.right, candidates);
        }
    }
}

#endif /* BEZIER_INTERSECTION_UTIL_H */