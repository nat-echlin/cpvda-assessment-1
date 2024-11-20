#ifndef CSVIO_H
#define CSVIO_H

#include <vector>
#include <iostream>

namespace utils {
    namespace csv {
        void readCsv(std::istream& in, std::vector<std::vector<double>>& data);
        void outputToTxt(std::string filename, const std::vector<double>& data);
    }
}

#endif
