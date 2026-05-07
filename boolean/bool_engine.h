#ifndef BOOL_ENGINE_H
#define BOOL_ENGINE_H

#include "../data_structures/geo.h"
#include "./sweep_line.h"
#include "./bool_options.h"
#include <optional>
#include <unordered_map>
#include <unordered_set>
#include <deque>
#include <algorithm>
#include <memory>
// #include <omp.h>

/*
This file will contain the main functionalities for the boolean operations
*/

struct InLine {
    int R;
    int L;
    std::shared_ptr<geom::Line> line;
};

class InPolygonChecker {
private:
    std::vector<std::shared_ptr<geom::Line>> lines;
    std::vector<std::shared_ptr<geom::Polygon>> polygons;
    std::vector<std::vector<bool>> pointInPolygons;

    std::vector<std::vector<InLine>> lineInPolygonMap;

public:
    InPolygonChecker(std::vector<std::shared_ptr<geom::Line>>& lines, std::vector<std::shared_ptr<geom::Polygon>>& polygons) :
          lines(lines), polygons(polygons) {}

    void execute() {
        // Create the middle points on left and right side on edge
        std::vector<std::shared_ptr<geom::Point>> checkPoints;
        for (const auto& line : lines) {
            std::vector<std::shared_ptr<geom::Point>> tmp = line->getMiddleLeftRight();
            checkPoints.push_back(tmp[0]);
            checkPoints.push_back(tmp[1]);
        }
        
        pointInPolygons.resize(checkPoints.size());
        for (size_t i = 0; i < checkPoints.size(); ++i) {
            pointInPolygons[i].resize(polygons.size());
        }

        // #pragma omp parallel for collapse(2)
        // #pragma omp parallel for schedule(dynamic, 100)
        for (size_t i = 0; i < checkPoints.size(); ++i) {
            for (size_t j = 0; j < polygons.size(); ++j) {
                pointInPolygons[i][j] = isInPoly(checkPoints[i], polygons[j]);
            }
        }

        for (size_t i = 0; i < checkPoints.size()/2; ++i) {
            std::vector<InLine> lineCollector;
            for (size_t j = 0; j < polygons.size(); ++j) {
                InLine tmp;
                tmp.line = lines[i];

                tmp.L = pointInPolygons[2*i][j];
                tmp.R = pointInPolygons[(2*i) + 1][j];

                lineCollector.push_back(tmp);
            }
            lineInPolygonMap.push_back(lineCollector);
        }
    } 

    std::vector<std::vector<InLine>> getLineInPolygonMap() {
        return lineInPolygonMap;
    }

private:
    bool isInPoly(std::shared_ptr<geom::Point> p, std::shared_ptr<geom::Polygon> poly) {
        bool inside = false;

        const auto& lines = poly->getLines();
        for (auto& l : lines) {
            // Check if point is on edge
            if (l->isPointOnSegment(p)) return true;

            // Ray-casting algorithm
            if ((l->s()->y() > p->y()) != (l->e()->y() > p->y())) {
                double xIntersect = (l->e()->x() - l->s()->x()) * (p->y() - l->s()->y()) / (l->e()->y() - l->s()->y()) + l->s()->x();
                if (p->x() < xIntersect) {
                    inside = !inside;
                }
            }
        }
        return inside;
    }
};

class BooleanEngine {
private:
    // Input properties
    std::vector<std::shared_ptr<geom::Polygon>> base;
    std::vector<std::shared_ptr<geom::Polygon>> clip;
    BooleanOptions type;

    // Output property
    std::vector<std::shared_ptr<geom::Polygon>> solu;

private: 
    // Inner private data structures
    // Input properties
    std::vector<std::shared_ptr<geom::Point>> allPoints;
    std::vector<std::shared_ptr<geom::Line>> allLines;

    // Line intersection properties
    // Used to quick filter, which lines has not to be checked
    std::unordered_map<std::shared_ptr<geom::Polygon>, std::unordered_set<std::shared_ptr<geom::Line>>> polygonLineMap;
    std::unordered_map<std::shared_ptr<geom::Line>, std::shared_ptr<geom::Polygon>> linePolygonMap;

    // Result of the sweep line method
    std::vector<std::shared_ptr<geom::Point>> intersectionPoints;
    
    // Line split properties
    std::vector<std::shared_ptr<geom::Line>> splitedLines;

    // RayTracing properties
    std::vector<std::vector<InLine>> lineInPolygonMap;

public:
    BooleanEngine(const std::shared_ptr<geom::Polygon>& p1, const std::shared_ptr<geom::Polygon>& p2, BooleanOptions type)
        : base{p1}, clip{p2}, type(type) {}

    BooleanEngine(const std::vector<std::shared_ptr<geom::Polygon>>& p1, const std::shared_ptr<geom::Polygon>& p2, BooleanOptions type)
        : base(p1), clip{p2}, type(type) {}

    BooleanEngine(const std::shared_ptr<geom::Polygon>& p1, const std::vector<std::shared_ptr<geom::Polygon>>& p2, BooleanOptions type)
        : base{p1}, clip(p2), type(type) {}

    BooleanEngine(const std::vector<std::shared_ptr<geom::Polygon>>& p1,
                const std::vector<std::shared_ptr<geom::Polygon>>& p2,
                BooleanOptions type)
        : base(p1), clip(p2), type(type) {}

    void execute() {
        // Get all the lines and store the data 
        collectLines();

        // Sort the lines
        // sortLines();

        // Sweep line to determine intersection points 
        intersectionPoints = sweepLineIntersection(allLines, linePolygonMap, polygonLineMap);
        intersectionPoints.insert(intersectionPoints.end(), allPoints.begin(), allPoints.end());
        uniqueIntersectionPoints();
        
        // Create the splited lines 
        splitLines();

        // Check which line is inside which polygon
        fillRayTraceResults();

        // Switch logic based on type
        switch (type) {
            case (BooleanOptions::Union):
                doUnion();
                break;
            case (BooleanOptions::Intersection):
                doIntersection();
                break;
            case (BooleanOptions::Difference):
                doDifference();
                break;
            case (BooleanOptions::Xor):
                doXor();
                break;
        }

        deleteLines();
    }

    std::vector<std::shared_ptr<geom::Polygon>> getSolu() {
        return solu;
    }

    ~BooleanEngine() = default;

private:
    void collectLines() {
        std::unordered_set<std::shared_ptr<geom::Line>> lineSet;

        auto processPolygons = [&](const std::vector<std::shared_ptr<geom::Polygon>>& polys){
            for (auto poly : polys) {
                poly->determineLines();
                std::vector<std::shared_ptr<geom::Line>> tmp = poly->getLines();
                std::unordered_set<std::shared_ptr<geom::Line>> lines(tmp.begin(), tmp.end());
                polygonLineMap[poly] = std::move(lines);
                lineSet.insert(polygonLineMap[poly].begin(), polygonLineMap[poly].end());
                for (std::shared_ptr<geom::Line> line : polygonLineMap[poly])
                    linePolygonMap[line] = poly;
            }
        };
        processPolygons(base);
        processPolygons(clip);

        allLines.insert(allLines.end(), lineSet.begin(), lineSet.end());

        getAllPoints();
    }

    void deleteLines() {
        auto processPolygons = [&](const std::vector<std::shared_ptr<geom::Polygon>>& polys){
            for (auto poly : polys) {
                poly->deleteLines();
            }
        };
        processPolygons(base);
        processPolygons(clip);
    }
    
    void getAllPoints() {
        std::unordered_set<std::shared_ptr<geom::Point>> uniquePoints; 
        for (size_t i = 0; i < allLines.size(); ++i) {
            uniquePoints.insert(allLines[i]->s());
            uniquePoints.insert(allLines[i]->e());
        }

        allPoints.clear();
        allPoints.insert(allPoints.end(), uniquePoints.begin(), uniquePoints.end());
    }

    void uniqueIntersectionPoints() {
        std::sort(intersectionPoints.begin(), intersectionPoints.end(),
          [](const std::shared_ptr<geom::Point>& a, const std::shared_ptr<geom::Point>& b) {
              if (a->X() != b->X()) return a->X() < b->X();
              return a->Y() < b->Y();
          });

        // Loop over the points and check if the next point is within the same
        // tolerance
        std::vector<std::shared_ptr<geom::Point>> unique;
        for (size_t i = 0; i < intersectionPoints.size(); ++i) {
            auto& p = intersectionPoints[i];
            if (!p) continue;

            for (size_t j = i+1; j < intersectionPoints.size(); ++j) {
                auto& q = intersectionPoints[j];
                if (!q) continue;

                if (q->X() > p->X() + geom::SAME_POINT) {
                    break;
                }

                if (p->isEqual(q)) {
                    p->parentLines.insert(q->parentLines.begin(), q->parentLines.end());
                    intersectionPoints[j] = nullptr;
                }
            }

            unique.push_back(p);
        }
        
        allPoints.clear();
        allPoints.insert(allPoints.end(), unique.begin(), unique.end());
    }

    void uniqueSplitedLines() {
        // Merge the lines with the same ends
        std::sort(splitedLines.begin(), splitedLines.end(),
            [](const std::shared_ptr<geom::Line>& a,
            const std::shared_ptr<geom::Line>& b) {

                if (a->s()->x() != b->s()->x())
                    return a->s()->x() < b->s()->x();
                if (a->s()->y() != b->s()->y())
                    return a->s()->y() < b->s()->y();
                if (a->e()->x() != b->e()->x())
                    return a->e()->x() < b->e()->x();
                return a->e()->y() < b->e()->y();
            }
        );

        std::vector<std::shared_ptr<geom::Line>> keepSplited;
        keepSplited.reserve(splitedLines.size());

        auto sameLine = [](const std::shared_ptr<geom::Line>& a,
                           const std::shared_ptr<geom::Line>& b) {
            return a->s() == b->s() && a->e() == b->e();
        };

        for (const auto& line : splitedLines) {
            if (keepSplited.empty() ||
                !sameLine(keepSplited.back(), line)) {
                keepSplited.push_back(line);
            }
        }

        splitedLines = std::move(keepSplited);
    }

    void updateParentLines() {
        for (size_t i = 0; i < splitedLines.size(); ++i) {
            splitedLines[i]->s()->parentLines.clear();
            splitedLines[i]->e()->parentLines.clear();
        }

        for (size_t i = 0; i < splitedLines.size(); ++i) {
            splitedLines[i]->s()->parentLines.insert(splitedLines[i]);
            splitedLines[i]->e()->parentLines.insert(splitedLines[i]);
        }
    }

    void splitLines() {
        std::unordered_map<std::shared_ptr<geom::Line>, std::unordered_set<std::shared_ptr<geom::Point>>> linePointMap;

        for (size_t i = 0; i < allLines.size(); ++i) {
            std::unordered_set<std::shared_ptr<geom::Point>> tmp;
            linePointMap[allLines[i]] = tmp;
        } 

        for (size_t i = 0; i < allPoints.size(); ++i) {
            std::shared_ptr<geom::Point> current = allPoints[i];
            for (auto& line : current->parentLines) {
                linePointMap[line].insert(current);
            }
        }

        for (size_t i = 0; i < allLines.size(); ++i) {
            std::vector<std::shared_ptr<geom::Point>> tmpPoints;
            tmpPoints.insert(tmpPoints.end(), linePointMap[allLines[i]].begin(), linePointMap[allLines[i]].end());
            
            if (allLines[i]->generalSort()) {
                std::sort(tmpPoints.begin(), tmpPoints.end(),
                    [](const std::shared_ptr<geom::Point>& a, const std::shared_ptr<geom::Point>& b) {
                        return a->x() < b->x();
                    });
            } else {
                std::sort(tmpPoints.begin(), tmpPoints.end(),
                    [](const std::shared_ptr<geom::Point>& a, const std::shared_ptr<geom::Point>& b) {
                        return a->y() < b->y();
                    });
            }

            // Create the lines and fill the map
            for (size_t j = 0; j < tmpPoints.size() - 1; ++j) {
                std::shared_ptr<geom::Line> split = std::make_shared<geom::Segment>(tmpPoints[j], tmpPoints[(j+1)%tmpPoints.size()]);
                splitedLines.push_back(split);
            }
        } 

        uniqueSplitedLines();
        updateParentLines();
    }

    void fillRayTraceResults() {
        std::vector<std::shared_ptr<geom::Polygon>> allPoly;
        allPoly.insert(allPoly.end(), base.begin(), base.end());
        allPoly.insert(allPoly.end(), clip.begin(), clip.end());

        std::unique_ptr<InPolygonChecker> inPoly = 
            std::make_unique<InPolygonChecker>(splitedLines, allPoly);
        inPoly->execute();
        lineInPolygonMap = inPoly->getLineInPolygonMap();
    }

    std::vector<std::vector<std::shared_ptr<geom::Point>>> buildPaths(std::unordered_set<std::shared_ptr<geom::Line>>& lines) {
        // Since multiple contour can be created store the result in vector
        std::vector<std::vector<std::shared_ptr<geom::Point>>> contours;

        int32_t count = 0;
        auto it = lines.begin();
        std::shared_ptr<geom::Line> start = *it;

        std::vector<std::shared_ptr<geom::Point>> tmp;
        tmp.push_back(start->s());
        tmp.push_back(start->e());

        int32_t pathSize = 2;

        lines.erase(it);
        contours.push_back(tmp);

        // Exclude logic to avoid false cases:
        int prevCount = 0;
        int prevPathSize = 0;
        while (!lines.empty()) {
            if (count == prevCount && pathSize == prevPathSize) {
                // Dead end the contour is not closed, drop it and 
                // start again
                prevCount = 0;
                prevPathSize = 0;

                contours[count].clear();

                auto it = lines.begin();
                std::shared_ptr<geom::Line> startContour = *it;
                contours[count].push_back(startContour->s());
                contours[count].push_back(startContour->e());

                lines.erase(it);
                pathSize = 2;
            }

            prevCount = count;
            prevPathSize = pathSize;
            std::shared_ptr<geom::Point> check = contours[count][pathSize-1];
            for (const auto& l : check->parentLines) {
                if (lines.find(l) != lines.end()) {
                    if (l->s() == check) {
                        contours[count].push_back(l->e());
                    }
                    if (l->e() == check) {
                        contours[count].push_back(l->s());
                    }

                    lines.erase(l);
                    pathSize++;

                    if (lines.empty()) {
                        break;
                    }

                    if (contours[count][0] == contours[count][pathSize-1]) {
                        // Loop closed, start a new contour
                        std::vector<std::shared_ptr<geom::Point>> newContour;
                        auto it = lines.begin();
                        std::shared_ptr<geom::Line> startContour = *it;
                        newContour.push_back(startContour->s());
                        newContour.push_back(startContour->e());

                        lines.erase(it);
                        contours.push_back(newContour);
                        pathSize = 2;
                        count++;
                    }
                }
            }
            
        }

        return contours;
    }

    void doUnion() {
        std::unordered_set<std::shared_ptr<geom::Line>> neededLines;
        for (size_t i = 0; i < lineInPolygonMap.size(); ++i) {
            int L = 0, R = 0;
            for (size_t j = 0; j < lineInPolygonMap[i].size(); ++j) {
                InLine tmp = lineInPolygonMap[i][j];
                L += tmp.L;
                R += tmp.R;
            }

            if (L != R && (L == 0 || R == 0)) {
                neededLines.insert(lineInPolygonMap[i][0].line);
            }
        }

        if (neededLines.empty()) {
            return;
        }

        std::vector<std::vector<std::shared_ptr<geom::Point>>> contours = buildPaths(neededLines);

        // First approach, just write the contours as solu
        // TODO: create a proper logic to rebuild the hierarchy of the contours and collect material and perforations
        for (const auto& data : contours) {
            solu.push_back(std::make_shared<geom::Polygon>(std::make_shared<geom::Path>(data)));
        }
    }

    void doIntersection() {
        int64_t nrPoly = base.size() + clip.size();
        std::unordered_set<std::shared_ptr<geom::Line>> neededLines;
        for (size_t i = 0; i < lineInPolygonMap.size(); ++i) {
            int L = 0, R = 0;
            for (size_t j = 0; j < lineInPolygonMap[i].size(); ++j) {
                InLine tmp = lineInPolygonMap[i][j];
                L += tmp.L;
                R += tmp.R;
            }

            if ((L > 0 && R > 0 && L != R) || (L == nrPoly || R == nrPoly)) {
                neededLines.insert(lineInPolygonMap[i][0].line);
            }
            if (L > 0 && R > 0 && L == R && L == nrPoly) {
                // Add, since the edge is in both poly's contour
                neededLines.insert(lineInPolygonMap[i][0].line);
            }
        }

        if (neededLines.empty()) {
            return;
        }

        std::vector<std::vector<std::shared_ptr<geom::Point>>> contours = buildPaths(neededLines);

        // First approach, just write the contours as solu
        // TODO: create a proper logic to rebuild the hierarchy of the contours and collect material and perforations
        for (const auto& data : contours) {
            solu.push_back(std::make_shared<geom::Polygon>(std::make_shared<geom::Path>(data)));
        }
    }

    void doDifference() {
        // Keep segments inside A, but not in B
        std::unordered_set<std::shared_ptr<geom::Line>> neededLines;

        // Collect the lines, which are contour base
        for (size_t i = 0; i < lineInPolygonMap.size(); ++i) {
            for (size_t j = 0; j < base.size(); ++j) {
                InLine tmp = lineInPolygonMap[i][j];
                if (tmp.L != tmp.R && (tmp.L == 0 || tmp.R == 0)) {
                    neededLines.insert(lineInPolygonMap[i][0].line);
                }
            }
        }
        // Exclude lines, which are inside contour clip
        for (size_t i = 0; i < lineInPolygonMap.size(); ++i) {
            for (size_t j = base.size(); j < clip.size() + base.size(); ++j) {
                InLine tmp = lineInPolygonMap[i][j];
                if (tmp.L == tmp.R && tmp.L != 0) {
                    neededLines.erase(lineInPolygonMap[i][0].line);
                }
                if (tmp.L != tmp.R && (tmp.L == 0 || tmp.R == 0)) {
                    neededLines.erase(lineInPolygonMap[i][0].line);
                }
            }
        }
        // Collect the lines, which are contour of clip and inside base or contour of base
        for (size_t i = 0; i < lineInPolygonMap.size(); ++i) {
            int LBase = 0, RBase = 0, LClip = 0, RClip = 0;
            for (size_t j = 0; j < lineInPolygonMap[i].size(); ++j) {
                InLine tmp = lineInPolygonMap[i][j];
                if (j < base.size()) {
                    LBase += tmp.L;
                    RBase += tmp.R;
                } else {
                    LClip += tmp.L; 
                    RClip += tmp.R;
                }

                if (LClip != RClip && (LClip == 0 || RClip == 0)) {
                    // Those lines, which are on contour of clip polygons
                    if (LBase == RBase && LBase != 0) {
                        // Those lines from set above, which inside base polygons
                        neededLines.insert(lineInPolygonMap[i][0].line);
                    }
                    if (LBase != RBase && (LBase == 0 || RBase == 0)) {
                        neededLines.insert(lineInPolygonMap[i][0].line);
                    }
                }
            }
        }
        
        if (neededLines.empty()) {
            return;
        }

        std::vector<std::vector<std::shared_ptr<geom::Point>>> contours = buildPaths(neededLines);

        // First approach, just write the contours as solu
        // TODO: create a proper logic to rebuild the hierarchy of the contours and collect material and perforations
        for (const auto& data : contours) {
            solu.push_back(std::make_shared<geom::Polygon>(std::make_shared<geom::Path>(data)));
        }
    }

    void doDifferenceFliped() {
        // Keep segments inside A, but not in B
        std::unordered_set<std::shared_ptr<geom::Line>> neededLines;

        // Collect the lines, which are contour clip
        for (size_t i = 0; i < lineInPolygonMap.size(); ++i) {
            for (size_t j = base.size(); j < base.size() + clip.size(); ++j) {
                InLine tmp = lineInPolygonMap[i][j];
                if (tmp.L != tmp.R && (tmp.L == 0 || tmp.R == 0)) {
                    neededLines.insert(lineInPolygonMap[i][0].line);
                }
            }
        }
        // Exclude lines, which are inside contour clip
        for (size_t i = 0; i < lineInPolygonMap.size(); ++i) {
            for (size_t j = 0; j < base.size(); ++j) {
                InLine tmp = lineInPolygonMap[i][j];
                if (tmp.L == tmp.R && tmp.L != 0) {
                    neededLines.erase(lineInPolygonMap[i][0].line);
                }
                if (tmp.L != tmp.R && (tmp.L == 0 || tmp.R == 0)) {
                    neededLines.erase(lineInPolygonMap[i][0].line);
                }
            }
        }
        // Collect the lines, which are contour of clip and inside base
        for (size_t i = 0; i < lineInPolygonMap.size(); ++i) {
            int LBase = 0, RBase = 0, LClip = 0, RClip = 0;
            for (size_t j = 0; j < lineInPolygonMap[i].size(); ++j) {
                InLine tmp = lineInPolygonMap[i][j];
                if (j < base.size()) {
                    LBase += tmp.L;
                    RBase += tmp.R;
                } else {
                    LClip += tmp.L; 
                    RClip += tmp.R;
                }

                if (LBase != RBase && (LBase == 0 || RBase == 0)) {
                    // Those lines, which are on contour of clip polygons
                    if (LClip == RClip && LClip != 0) {
                        // Those lines from set above, which inside base polygons
                        neededLines.insert(lineInPolygonMap[i][0].line);
                    }
                    if (LClip != RClip && (LClip == 0 || RClip == 0)) {
                        neededLines.insert(lineInPolygonMap[i][0].line);
                    }
                }
            }
        }
        
        if (neededLines.empty()) {
            return;
        }

        std::vector<std::vector<std::shared_ptr<geom::Point>>> contours = buildPaths(neededLines);

        // First approach, just write the contours as solu
        // TODO: create a proper logic to rebuild the hierarchy of the contours and collect material and perforations
        for (const auto& data : contours) {
            solu.push_back(std::make_shared<geom::Polygon>(std::make_shared<geom::Path>(data)));
        }
    }

    void doXor() {
        // Keep segments inside exactly one polygon.
        doDifference();
        doDifferenceFliped();
    }
};

class PolygonBuilder {
public:
    PolygonBuilder() {}
};

#endif // BOOL_ENGINE_H