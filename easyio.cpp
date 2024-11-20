#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <vector>

#include "easyio.h"

namespace utils {
    namespace easyio {

        void readCsv(const std::string& fileName, std::vector<std::vector<double>>& data) {
            std::ifstream file(fileName);
            if (!file.is_open()) {
                std::cerr << "Error: Could not open file " << fileName << std::endl;
                exit(1);
            }

            std::string line;
            while (getline(file, line)) {
                data.push_back({}); // Add a new row
            
                std::istringstream lineIn(line);
                double val;
                
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

            file.close();
        }

        void outputToTxt(const std::string filename, const std::vector<double>& data) {
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