#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <vector>

#include "easyio.h"

namespace utils {
    namespace easyio {

        // Writes by reference
        void readCsv(std::istream& in, std::vector<std::vector<double>>& data) {
            std::string line;
            
            while (getline(in, line)) {
                data.push_back({}); // Add a new row
            
                std::istringstream lineIn(line);
                int val;
                
                while (lineIn >> val) {
                    data.back().push_back(val);
            
                    char c;
                    lineIn >> c;
                    
                    if (lineIn && c != ',') {
                        std::cout << "Expected a comma - aborting." << std::endl;
                        std::cout << c << std::endl;
                        exit(0);
                    }
                }
            }
        }

        void outputToTxt(std::string filename, const std::vector<double>& data) {
            std::ofstream out(filename);
            if (!out.is_open()) {
                std::cerr << "Failed to open file for writing." << std::endl;
            }

            for (int i = 0; i < data.size(); ++i) {
                out << data[i] << std::endl;
            }
        }
    }
}