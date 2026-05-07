#include "data_structures/geo.h"

int main() {
    double x[] = {0, 10, 10, 0, 0};
    double y[] = {0, 0, 10, 10, 0};

    // Create points
    std::vector<std::shared_ptr<geom::Point>> points;
    for (int i = 0; i < 5; ++i) {
        points.push_back(std::make_shared<geom::Point>(x[i], y[i]));
    }

    // Fill points to Path
    std::shared_ptr<geom::Path> path = std::make_shared<geom::Path>(points);

    // Write out the points
    for (std::shared_ptr<geom::Point> p : path->points()) {
        std::cout << "x: " << p->x() << ", y: " << p->y() << std::endl;
    }

    return 0;
}