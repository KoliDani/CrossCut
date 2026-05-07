#ifndef PROJECTED_GAUSS_NEWTON_H
#define PROJECTED_GAUSS_NEWTON_H

#include <cmath>
#include <algorithm>
#include <iostream>
#include "../data_structures/point.h"

/*
    This file implements the Projected Gauss-Newton method for finding the intersection of two Bezier curves.
    The algorithm iteratively refines parameter estimates for points on each curve to minimize the distance between them, 
    while ensuring that the parameters remain within valid ranges.
    The main class is GaussNewtonSolver, which takes control points of two Bezier curves and an initial candidate for intersection as input, 
    and provides a solve() method to compute the intersection point.
*/

// Constants for the Projected Gauss-Newton algorithm, including maximum iterations, tolerance for convergence, and step size.
namespace PROJ_GAUSS_NEWTON {
    // Original MAX_ITERATION = 20;
    const int MAX_ITERATION = 10;
    const double TOL_RES = 1e-5;
    const double TOL_STEP = 1e-5;
    const double TOL_SAME_POINT = 1e-1;
}

// Namespace for utility functions related to Bezier curve intersection, 
// including data structures for subdividing curves, calculating pseudo-curvature, and checking bounding box intersections.
namespace BEZIER_INTERSECTION_UTIL {
    struct ParamRange {
        double tMin;
        double tMax;
    };

    struct IntersectionCanditate {
        ParamRange a; // bezier 1 specific property
        ParamRange b; // bezier 2 specific property
    };
}

// Function to calculate the Euclidean norm of a 2D vector, used to measure the distance between points in the Gauss-Newton algorithm.
inline double norm(const Vec2& v) {
    return std::sqrt(v.x*v.x + v.y*v.y);
}

// Distance calculation between two 2D points represented as Vec2, used for measuring the distance between points in the algorithm.
inline double distanceTo(const Vec2& a, const Vec2& b) {
    double dx = a.x - b.x;
    double dy = a.y - b.y;
    return std::sqrt(dx*dx + dy*dy);
}

/*
    GaussNewtonSolver class implementation for 2D Bezier curve intersection. It uses the Projected Gauss-Newton method to iteratively find the 
    parameters on each curve that minimize the distance between them, effectively finding the intersection point.
*/
class GaussNewtonSolver {
private:
    // Control points of the two Bezier curves and the initial candidate for intersection
    std::vector<Vec2> controlA;
    std::vector<Vec2> controlB;
    BEZIER_INTERSECTION_UTIL::IntersectionCanditate candidate;

public:
    // Constructor to initialize the solver with the control points of the two curves and the initial candidate for intersection
    GaussNewtonSolver(const std::vector<Vec2>& ctrlA, 
        const std::vector<Vec2>& ctrlB, BEZIER_INTERSECTION_UTIL::IntersectionCanditate candidate)
        : controlA(ctrlA), controlB(ctrlB), candidate(candidate) {}

    // Main method to solve for the intersection point using the Projected Gauss-Newton algorithm.
    // It returns a shared pointer to the intersection point if found, or nullptr if the algorithm fails to converge or the result is not valid.
    std::shared_ptr<geom::Point> solve() {
        double xi = (candidate.a.tMin + candidate.a.tMax) / 2;
        double eta = (candidate.b.tMin + candidate.b.tMax) / 2;

        double xmin = candidate.a.tMin, xmax = candidate.a.tMax;
        double ymin = candidate.b.tMin, ymax = candidate.b.tMax;

        if (projectedGaussNewton(xi, eta, xmin, xmax, ymin, ymax)) {
            return getPoint(xi, eta); // CA(xi) or CB(eta) should be the same point;
        } else {
            return nullptr; // Failed to converge
        }
    }

private:
    // Methods to evaluate the Bezier curves and their derivatives at given parameter values,
    // which are used in the Gauss-Newton iterations to compute the function values and Jacobian matrix for the optimization step.
    Vec2 CA(double t) {
        return evaluate(controlA, t);
    }

    Vec2 CB(double t) {
        return evaluate(controlB, t);
    }

    // Methods to compute the derivatives of the Bezier curves at given parameter values, 
    // which are needed to construct the Jacobian matrix for the Gauss-Newton step.
    Vec2 dCA(double t) {
        return derivate(controlA, t);
    }

    Vec2 dCB(double t) {
        return derivate(controlB, t);
    }

    // Nonlinear system
    Vec2 F(double xi, double eta) {
        Vec2 a = CA(xi);
        Vec2 b = CB(eta);
        return {a.x - b.x, a.y - b.y};
    }

    // Solve 2x2 system for Gauss-Newton step
    bool solveStep(double xi, double eta, double& dxi, double& deta) {
        Vec2 da = dCA(xi);
        Vec2 db = dCB(eta);

        double J11 = da.x, J12 = -db.x;
        double J21 = da.y, J22 = -db.y;

        double det = J11*J22 - J12*J21;
        if (std::abs(det) < 1e-12) return false; // singular

        Vec2 f = F(xi, eta);
        dxi  = ( J22*f.x - J12*f.y) / det;
        deta = (-J21*f.x + J11*f.y) / det;
        return true;
    }

    // Project parameters back to valid range
    void project(double& xi, double& eta, double xmin, double xmax, double ymin, double ymax) {
        xi  = std::max(xmin, std::min(xi, xmax));
        eta = std::max(ymin, std::min(eta, ymax));
    }

    // Projected Gauss-Newton solver implementation that iteratively refines the parameter estimates for the intersection point.
    bool projectedGaussNewton(
        double& xi, double& eta,
        double xmin, double xmax,
        double ymin, double ymax) {

        for (int k = 0; k < PROJ_GAUSS_NEWTON::MAX_ITERATION; ++k) {
            Vec2 f = F(xi, eta);
            double r = norm(f);
            if (r < PROJ_GAUSS_NEWTON::TOL_RES) return true;

            double dxi, deta;
            if (!solveStep(xi, eta, dxi, deta)) return false;

            double newXi  = xi - dxi;
            double newEta = eta - deta;

            project(newXi, newEta, xmin, xmax, ymin, ymax);

            double step = std::sqrt((newXi-xi)*(newXi-xi) + (newEta-eta)*(newEta-eta));
            xi = newXi; eta = newEta;

            if(step < PROJ_GAUSS_NEWTON::TOL_STEP) return true;
        }
        return false;
    }

    // Static method to evaluate a Bezier curve at a given parameter t using De Casteljau's algorithm,
    // which is a numerically stable way to compute points on a Bezier curve based on its control points.
    static Vec2 evaluate(const std::vector<Vec2>& control, double t) {
        std::vector<Vec2> temp;
        temp.reserve(control.size());
        for (const auto& p : control) {
            temp.push_back(Vec2{p.x, p.y});
        }
        for (size_t k = temp.size() - 1; k > 0; --k) {
            for (size_t i = 0; i < k; ++i) {
                double x = (1 - t) * temp[i].x + t * temp[i + 1].x;
                double y = (1 - t) * temp[i].y + t * temp[i + 1].y;
                temp[i] = Vec2{x, y};
            }
        }
        return temp[0];
    }

    // Static method to compute the derivative of a Bezier curve at a given parameter t, which is needed for the Gauss-Newton iterations.
    static Vec2 derivate(const std::vector<Vec2>& control, double t) {
        int n = control.size() - 1;
        if(n == 0) return Vec2{0.0, 0.0};

        std::vector<Vec2> dP(n);
        for(int i = 0; i < n; ++i){
            dP[i] = Vec2{
                n * (control[i+1].x - control[i].x),
                n * (control[i+1].y - control[i].y)
            };
        }

        return evaluate(dP, t);
    }

    // This function can be used to return the intersection point after convergence
    // by evaluating either CA(xi) and CB(eta). It checks if the two points are close enough to be considered the 
    // same point within a certain tolerance, and returns the average of the two points for better accuracy.
    // If the points are not within the tolerance, it returns nullptr to indicate that the result is not valid.
    std::shared_ptr<geom::Point> getPoint(double xi, double eta) {
        // This function can be used to return the intersection point after convergence
        // by evaluating either CA(xi) and CB(eta)

        Vec2 ptA = evaluate(controlA, xi);
        Vec2 ptB = evaluate(controlB, eta);

        if (distanceTo(ptA, ptB) <= PROJ_GAUSS_NEWTON::TOL_SAME_POINT) {
            // Average of the two points for better accuracy
            return std::make_shared<geom::Point>((ptA.x + ptB.x) / 2, (ptA.y + ptB.y) / 2); 
        } else {
            // This happens, when the algorithm converges to a point where CA and CB are close but not within the tolerance.
            return nullptr; 
        }
    }
};

#endif // PROJECTED_GAUSS_NEWTON_H