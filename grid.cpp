#include <vector>
#include <iostream>
#include <utility>
#include <omp.h>

#include "easyio.h"
#include "rng.h"

using namespace std;

typedef vector<vector<vector<pair<double, double>>>> Grid;


Grid buildGrid(vector<double> x, vector<double> y) {
    size_t n = x.size();
    
    // for calculating furthests we need that the furthest point will be in one of the corners almost surely
    // so pts per cell needs to be ~100
    // do some probability stuff on this to determine how many points this needs for us to be sure of it
    int ptsPerCell = 100;
    int cellCount = n / ptsPerCell;
    int groupsPerAxis = sqrt(cellCount);

    double cellSize = 1 / sqrt(cellCount); 
    
    // initialise how big the grid should be
    Grid grid(groupsPerAxis, vector<vector<pair<double, double>>>(groupsPerAxis, vector<pair<double, double>>()));

    for (size_t i = 0; i < n; ++i) {
        // use emplace to build the pair inside the vector, a bit faster
        int xInd = x[i] / cellSize;
        int yInd = y[i] / cellSize;
        
        grid[xInd][yInd].emplace_back(x[i], y[i]);        
    }

    return grid;
}

// Assumes that all cells in the grid have at least one point
void calculateDistancesWithGrid(Grid grid) {
    int cellCountPerD = grid.size();

    if (cellCountPerD != grid[0].size()) {
        throw invalid_argument("x and y must have the same number of elements");
    }
    
    // calculate nearests
    vector<double> nearests = {};
    
    // iter through all grids
    for (int i = 0; i < cellCountPerD; ++i) {
        for (int j = 0; j < cellCountPerD; ++j) {
            
            // iter through all points
            double runningNearest = 1.5; // arbitrary number > sqrt(2)
            
            for (int k = 0; k < grid[i][j].size(); ++i) {
                
                // now again iter through relevant grids s.t. we check all neighbour grids
                for (int l = max(0, i-1); l <= min(cellCountPerD, i+1); ++l) {
                    for (int m = max(0, j-1); m <= min(cellCountPerD, j+1); ++m) {
                        
                        // now iter through points in that grid
                        for (int n = 0; n < grid[l][m].size(); ++n) {
                            // compare point n to point k
                            double xdiff = abs(grid[i][j][k].first - grid[l][m][n].first);
                            double ydiff = abs(grid[i][j][k].second - grid[l][m][n].second);
                        
                            double dist = sqrt(pow(ydiff, 2) + pow(xdiff, 2));

                            if (dist < runningNearest) {
                                runningNearest = dist;
                            }
                        }
                    }
                }
            }
            
            nearests.push_back(runningNearest);
        }
    }    
    
    // print out summary here
    // and write extended results to txt
}   
