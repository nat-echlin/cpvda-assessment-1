#ifndef EASYIO_H
#define EASYIO_H

#include <vector>
#include <iostream>

namespace utils {
    namespace easyio {
        void readCsv(std::istream& in, std::vector<std::vector<double>>& data);
        void outputToTxt(std::string filename, const std::vector<double>& data);
    }
}

#endif
