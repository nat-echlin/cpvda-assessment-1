#include <vector>
#include <iostream>
#include <fstream>
#include <chrono>
#include <utility>
#include <omp.h>

#include "easyio.h"
#include "rng.h"

// small program, no danger of namespace pollution. Makes code more readable, good to do.
using namespace std; 

// taken from https://stackoverflow.com/questions/22387586/measuring-execution-time-of-a-function-in-c
// changed slightly to give result in seconds rather than nanoseconds
typedef chrono::high_resolution_clock::time_point TimeVar;
#define duration(a) chrono::duration<double>(a).count()
#define timeNow() chrono::high_resolution_clock::now()

template<typename F, typename... Args>
void funcTime(F func, Args&&... args){
    TimeVar t1=timeNow();
    func(std::forward<Args>(args)...);
    cout << "time taken: " << duration(timeNow()-t1) << endl;
}
// end stackoverflow reference

// Returns 100k randomly initialised points, from 0.0 to 1.0 in the x & y direction. 
// In the format of {{x_1, x_2, ..., x_n},  {y_1, y_2, ..., x_n}}
// int n : the number of points to create
vector<vector<double>> randomPoints(int n) {
    double start = 0;
    double end = 1;
    
    vector<double> x(n, start);
    vector<double> y(n, start);

    for (int i = 0; i < n; ++i) {
        // gen a random x and y coord
        
        x[i] = utils::rand::randDouble(start, end);
        y[i] = utils::rand::randDouble(start, end);
    }

    return vector<vector<double>> {x, y};
}

// Calculates the shortest and longest distance, given the x and y distances, between two points.
// Places the shortest distance in arr[0], the longest distance in arr[1]. 
// 
// arr[0] must have size >= 2.
// Uses wraparound geometry.
void wraparoundGeo(double xdiff, double ydiff, double arr[]) {
    double closestDist = sqrt(pow(min(ydiff, 1 - ydiff), 2) + pow(min(xdiff, 1 - xdiff), 2));
    double furthestDist = sqrt(pow(max(ydiff, 1 - ydiff), 2) + pow(max(xdiff, 1 - xdiff), 2));

    arr[0] = closestDist;
    arr[1] = furthestDist;
}

// Calculates the shortest and longest distance, given the x and y distances, between two points.
// Places the shortest distance in arr[0], the longest distance in arr[1]. 
// 
// arr[0] must have size >= 2.
// Uses basic geometry.
void basicGeo(double xdiff, double ydiff, double arr[]) {
    double dist = sqrt(pow(ydiff, 2) + pow(xdiff, 2));
    
    arr[0] = dist;
    arr[1] = dist;
}

// For each given point (x_i, y_i), calculate the distance to the nearest and furthest point.
// For calcFunc, pass the function that describes the geometry to use. Basic and wraparound are
// both given in this file.
// 
// Nearest points will be written to nearest.txt, furthests to furthest.txt. 
// Prints to stdout the average nearest and furthest distances.
// Calculated in serial.
void calcNearestAndFurthestDistances_Serial(vector<double> x, vector<double> y, void (*calcFunc)(double, double, double[])) {
    int n = x.size();

    if (n != y.size()) {
        throw invalid_argument("x and y must have the same number of elements");
    }
    
    // setup
    vector<double> nearests(n, 0);
    vector<double> furthests(n, 0);
    
    double avgNearest = 0.0;
    double avgFurthest = 0.0;

    // contains {nearest dist, furthest dist}. No need to reinitialise ever, just reuse it.
    double res[2];
    
    // iter through all points
    for (int i = 0; i < n; ++i) {
        // calc for each point
        double runningNearest = 1.5; // arbitrary number > sqrt(2)
        double runningFurthest = -0.1; // arbitrary number < 0
        
        for (int j = 0; j < n; ++j) {
            // calculate euclidean distances in x & y
            double ydiff = abs(y[j] - y[i]);
            double xdiff = abs(x[j] - x[i]);

            // pass to calcFunc to compute nearest and furthest distances, basic dependency injection
            calcFunc(xdiff, ydiff, res);
            
            if (res[0] < runningNearest && i != j) {
                runningNearest = res[0];
            }
            if (res[1] > runningFurthest && i != j) {
                runningFurthest = res[1];
            }
        }

        nearests[i] = runningNearest;
        furthests[i] = runningFurthest;
    }

    double nearestMean = 0;
    double furthestMean = 0;
    for (int i = 0; i < n; ++i) {
        nearestMean += nearests[i];
        furthestMean += furthests[i];
    }
    avgNearest = nearestMean / n;
    avgFurthest = furthestMean / n;

    // write to plaintext file
    utils::easyio::outputToTxt("basicNearests.txt", nearests);
    utils::easyio::outputToTxt("basicFurthests.txt", furthests);

    // output averages
    cout << "mean nearest distance: " << avgNearest << endl
         << "mean furthest distance: " << avgFurthest << endl
         << "\n"; // to make it easier to read multiple outputs at once
}

// For each given point (x_i, y_i), calculate the distance to the nearest and furthest point.
// For calcFunc, pass the function that describes the geometry to use. Basic and wraparound are
// both given in this file.
// 
// Nearest points will be written to nearest.txt, furthests to furthest.txt. 
// Prints to stdout the average nearest and furthest distances.
// Calculated using parallelisation from OpenMP.
void calcNearestAndFurthestDistances_Parallel(vector<double> x, vector<double> y, void (*calcFunc)(double, double, double[])) {
    int n = x.size();

    if (n != y.size()) {
        throw invalid_argument("x and y must have the same number of elements");
    }
    
    // setup
    vector<double> nearests(n, 0);
    vector<double> furthests(n, 0);
    
    double avgNearest = 0.0;
    double avgFurthest = 0.0;

    // contains {nearest dist, furthest dist}. No need to reinitialise ever, just reuse it.
    double res[2];
    
    // iter through all points
    #pragma omp parallel for private(res) // make sure that each thread has its own threadbound`res`
    for (int i = 0; i < n; ++i) {
        // calc for each point
        double runningNearest = 1.5; // arbitrary number > sqrt(2)
        double runningFurthest = -0.1; // arbitrary number < 0
        
        for (int j = 0; j < n; ++j) {
            // calculate euclidean distances in x & y
            double ydiff = abs(y[j] - y[i]);
            double xdiff = abs(x[j] - x[i]);

            // pass to calcFunc to compute nearest and furthest distances, basic dependency injection
            calcFunc(xdiff, ydiff, res);
            
            if (res[0] < runningNearest && i != j) {
                runningNearest = res[0];
            }
            if (res[1] > runningFurthest && i != j) {
                runningFurthest = res[1];
            }
        }

        nearests[i] = runningNearest;
        furthests[i] = runningFurthest;
    }

    double nearestMean = 0;
    double furthestMean = 0;
    #pragma omp parallel for reduction(+:nearestSum, furthestSum) // prevent race conditions 
    for (int i = 0; i < n; ++i) {
        nearestMean += nearests[i];
        furthestMean += furthests[i];
    }
    avgNearest = nearestMean / n;
    avgFurthest = furthestMean / n;

    // write to plaintext file
    utils::easyio::outputToTxt("basicNearests.txt", nearests);
    utils::easyio::outputToTxt("basicFurthests.txt", furthests);

    // output averages
    cout << "mean nearest distance: " << avgNearest << endl
         << "mean furthest distance: " << avgFurthest << endl
         << "\n"; // to make it easier to read multiple outputs at once
}



int main () {
    int n = 1000;

    // vector<vector<double>> points_csvBas;
    // utils::easyio::readCsv("100000 locations.csv", points_csvBas);

    // for (int i = 0; i < 5; ++i) {
    //     cout << points_csvBas[0][0] << " " << points_csvBas[i][1] << endl;
    // }

    // examples

    // n randomly-initialised points, using basic geometry.
        // vector<vector<double>> points_randBas = randomPoints(n);
        
        // funcTime(
        //     calcNearestAndFurthestDistances_Serial, points_randBas[0], points_randBas[1], basicGeo
        // );

    // n randomly-initialised points, using wraparound geometry.
        // vector<vector<double>> points_randWrap = randomPoints(n);
        
        // funcTime(
        //     calcNearestAndFurthestDistances_Serial, points_randWrap[0], points_randWrap[1], wraparoundGeo
        // );

    // n csv-initialised points, using basic geometry.
        vector<vector<double>> points_csvBas;
        utils::easyio::readCsv("100000 locations.csv", points_csvBas);

        funcTime(
            calcNearestAndFurthestDistances_Serial, points_csvBas[0], points_csvBas[1], basicGeo
        );

    // n csv-initialised points, using wraparound geometry.
        vector<vector<double>> points_csvWrap;
        utils::easyio::readCsv("100000 locations.csv", points_csvWrap);

        funcTime(
            calcNearestAndFurthestDistances_Serial, points_csvWrap[0], points_csvWrap[1], wraparoundGeo
        );

    return 0;
}