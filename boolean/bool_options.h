#ifndef BOOL_OPTIONS_H
#define BOOL_OPTIONS_H

/*
    * This file defines the BooleanOptions enum, which represents different boolean operations
    * that can be performed on sets or collections. The options include Union, Intersection,
    * Difference, and Xor.
    *
*/
enum BooleanOptions {
    Union,        // 0
    Intersection, // 1
    Difference,   // 2
    Xor           // 3
};

#endif // BOOL_OPTIONS_H