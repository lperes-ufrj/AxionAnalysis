#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <sstream>
#include "TH2F.h"
#include "TFile.h"

// Function to convert unit vector to angles
void unit_vector_to_angles_cartesian(double x, double y, double z, double &theta_xz, double &theta_yz) {
    theta_xz = atan2(x, z);
    theta_yz = atan2(y, z);
}

int main() {
    // Open the CSV file
    std::ifstream file("DetectorCoordSys_directions_over_10_years.csv"); // Replace with your actual CSV file
    if (!file.is_open()) {
        std::cerr << "Error: Could not open the CSV file!" << std::endl;
        return 1;
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

    // Skip the header line if there is one
    std::getline(file, line); 

    while (std::getline(file, line) && line_count < 100) { // Read only first 100 rows
        std::stringstream ss(line);
        std::string value;
        std::vector<double> row;

        // Read comma-separated values
        while (std::getline(ss, value, ',')) {
            row.push_back(std::stod(value)); // Convert string to double
        }

        if (row.size() < 4) {
            std::cerr << "Error: Row with insufficient data." << std::endl;
            continue;
        }

        double x = row[1]; // Assuming x is in column 2 (index 1)
        double y = row[2]; // Assuming y is in column 3 (index 2)
        double z = row[3]; // Assuming z is in column 4 (index 3)

        double theta_xz, theta_yz;
        unit_vector_to_angles_cartesian(x, y, z, theta_xz, theta_yz);

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
    return 0;
}