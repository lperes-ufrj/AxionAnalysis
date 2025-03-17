#include <cmath>
#include <iostream>
#include <cstdlib>


#include "src/STClibrary.h"

int test() {
    // Example input values
    double latitude = 44.3;  // Observer's latitude
    double longitude = 103.70; // Observer's longitude
    double l = 0.0;        // Galactic longitude
    double b = 0.0;         // Galactic latitude
    double jd = 2459580.5;   // Julian Date
    double refjd = 2459580.0; // Reference Julian Date
    int eqn = 1;             // Example equation type

    double az, el;
    
    // Call the gal2azel function
    gal2azel(latitude, longitude, l, b, jd, &az, &el, refjd, eqn);

    // Print the results
    std::cout << "Azimuth: " << az << std::endl;
    std::cout << "Elevation: " << el << std::endl;

    return 0;
}