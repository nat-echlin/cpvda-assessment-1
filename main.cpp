#include <vector>
#include <iostream>
#include <fstream>
#include <chrono>
#include <utility>

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

// For each given point (x_i, y_i), calculate the distance to the nearest and furthest point.
// For calcFunc, pass the function that describes the geometry to use. Basic and wraparound are
// both given in this file.
// 
// Nearest points will be written to nearest.txt, furthests to furthest.txt. 
// Prints to stdout the average nearest and furthest distances.
void calcNearestAndFurthestDistances(vector<double> x, vector<double> y, void (*calcFunc)(double, double, double[])) {
    int n = x.size();
    
    vector<double> nearests(n, 0);
    vector<double> furthests(n, 0);
    
    double avgNearest = 0.0;
    double avgFurthest = 0.0;

    // contains {nearest dist, furthest dist}. No need to reinitialise ever, just reuse it.
    double res[2];
    
    // iter through all points
    for (int i = 0; i < n; ++i) {
        // calc for each point
        double runningNearest = 10.0;
        double runningFurthest = -10.0;

        
        for (int j = 0; j < n; ++j) {
            // using pythagoras to calculate the difference
            double ydiff = abs(y[j] - y[i]);
            double xdiff = abs(x[j] - x[i]);

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

        // recalc means
        // avgNearest = (avgNearest * i + runningNearest) / (double) (i + 1);
        // avgFurthest = (avgFurthest * i + runningFurthest) / (double) (i + 1);
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
         << "mean furthest distance: " << avgFurthest << endl;
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

int main () {
    int n = 100000;

    // 100k randomly initialised points using basic geometry
        // vector<vector<double>> points = randomPoints(n);
        
        // funcTime(
        //     calcNearestAndFurthestDistances, points[0], points[1], basicGeo
        // );

    // 100k randomly initialised points using wraparound geometry
        vector<vector<double>> points = randomPoints(n);
        
        funcTime(
            calcNearestAndFurthestDistances, points[0], points[1], wraparoundGeo
        );

    return 0;
}