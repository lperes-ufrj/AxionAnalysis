#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <sstream>
#include "TH2F.h"
#include "TFile.h"

// Function to convert spherical (polar) coordinates (theta, phi) to angles
void unit_vector_to_angles_polar(double theta, double phi, double &theta_xz, double &theta_yz) {
    theta_xz = atan2(sin(theta) * cos(phi), cos(theta));
    theta_yz = asin(sin(theta) * sin(phi));
}

void roothistfile() {
    // Open the CSV file
    std::ifstream file("DetectorCoordSys_directions_over_10_years.csv"); // Ensure this file exists and contains theta, phi data
    if (!file.is_open()) {
        std::cerr << "Error: Could not open the CSV file!" << std::endl;
        return;
    }

    // Create a ROOT file to store the histogram
    TFile *root_file_cart = new TFile("angle_distributions.root", "RECREATE");

    // Create a 2D histogram
    const char* hist_name = "theta_xz_yz";
    TH2F *hist = new TH2F(hist_name, "Theta XZ vs YZ Distribution", 200, -3.14, 3.14, 200, -3.14, 3.14);

    std::vector<double> ang_xz;
    std::vector<double> ang_yz;

    std::string line;
    int line_count = 0;

    // Skip the header line if needed
    std::getline(file, line);

    while (std::getline(file, line) && line_count < 100) { // Read only first 100 rows
        std::stringstream ss(line);
        std::string value;
        std::vector<double> row;

        // Read comma-separated values
        while (std::getline(ss, value, ',')) {
            row.push_back(std::stod(value)); // Convert string to double
        }

        if (row.size() < 6) { // Ensure at least 6 columns exist (since theta is in col 5, phi in col 4)
            std::cerr << "Error: Row with insufficient data." << std::endl;
            continue;
        }

        double phi = row[4];   // Assuming phi is in column 5 (index 4)
        double theta = row[5]; // Assuming theta is in column 6 (index 5)

        // Compute angles
        double theta_xz, theta_yz;
        unit_vector_to_angles_polar(theta, phi, theta_xz, theta_yz);

        hist->Fill(theta_xz, theta_yz);
        ang_xz.push_back(theta_xz);
        ang_yz.push_back(theta_yz);

        line_count++;
    }

    // Write histogram to the ROOT file
    hist->Write();
    root_file_cart->Close();

    std::cout << "Histogram '" << hist_name << "' saved in 'angle_distributions.root'" << std::endl;

    file.close();
}
