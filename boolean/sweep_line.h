#ifndef SWEEP_LINE_H
#define SWEEP_LINE_H

#include "../data_structures/geo.h"
#include <set>
#include <vector>
#include <memory>
#include <algorithm>
#include <cmath>
#include <unordered_map>
#include <unordered_set>

/*
Original code is modified in order to return the std::vector<std::shared_ptr<Point>>, which is the all
intersection point with the given line set
Later the original line points is also append to this vector, then after unique operation the final point 
set will be defined. Since intersection point will contain, which edges defined the point, we only need to
create a look-up data structure Original Line - Point, in this way, we can "split" the original Lines.
Algorithm: https://www.geeksforgeeks.org/dsa/given-a-set-of-line-segments-find-if-any-two-segments-intersect/
*/

// An event for the sweep line algorithm
struct Event {
    std::shared_ptr<geom::Point> point;
    bool isLeft;
    int index;
};

// Returns the predecessor iterator in set s.
std::set<Event>::iterator pred(std::set<Event> &s, std::set<Event>::iterator it) {
    auto tmp = it;
    return tmp == s.begin() ? s.end() : --tmp;
}

// Returns the successor iterator in set s.
std::set<Event>::iterator succ(std::set<Event> &s, std::set<Event>::iterator it) {
    auto tmp = it;
    return ++tmp;
}

// Returns the number of intersections found between segments.
std::vector<std::shared_ptr<geom::Point>> sweepLineIntersection(std::vector<std::shared_ptr<geom::Line>> &lines,
     std::unordered_map<std::shared_ptr<geom::Line>, std::shared_ptr<geom::Polygon>>& doNotCheck,
     std::unordered_map<std::shared_ptr<geom::Polygon>, std::unordered_set<std::shared_ptr<geom::Line>>>& polygonLineMap) {
    // Fill the parentLine property of the points
    for (size_t i = 0; i < lines.size(); i++) {
        lines[i]->s()->parentLines.insert(lines[i]);
        lines[i]->e()->parentLines.insert(lines[i]);
    }

    // To avoid duplicate intersection reports.
    std::vector<std::shared_ptr<geom::Point>> intersections; 
    std::vector<Event> events;

    // Build event queue 
    int i = 0;
    for (auto& l : lines) {
        // If the type of line is not segment, 
        // then instead of start and end place the bounding box min and max
        if (auto seg = std::dynamic_pointer_cast<geom::Segment>(l)) {
            events.push_back({l->s(), true, i});
            events.push_back({l->e(), false, i++});
        } else {
            BoundingBox box = l->getBoundingBox();
            events.push_back({box.getMin(), true, i});
            events.push_back({box.getMax(), false, i++});
        }
    }

    // Sort event queue
    std::sort(events.begin(), events.end(),
    [](const Event& a, const Event& b) {
        if (a.point->X() != b.point->X())
            return a.point->X() < b.point->X();
        return a.point->Y() < b.point->Y();
    });

    std::unordered_set<int> active;

    // Process all events.
    for (const Event& curr : events) {
        int index = curr.index;
        if (curr.isLeft) {
            for (int other : active) {
                if (polygonLineMap[doNotCheck[lines[other]]]
                        .find(lines[index]) !=
                    polygonLineMap[doNotCheck[lines[other]]].end())
                    continue;

                std::vector<std::shared_ptr<geom::Point>> pts = lines[index]->intersect(lines[other]);
                for (const auto& p : pts) {
                    intersections.push_back(p);
                }
            }
            active.insert(index);
        }
        else {
            active.erase(index);
        }
    }

    return intersections;
}

#endif // SWEEP_LINE_H
