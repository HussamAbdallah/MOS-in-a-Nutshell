#pragma once
#pragma once


#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string.h>


using namespace std;

namespace MyApp
{
    void func();

    class Interpolator {
    public:
        // Linear interpolation function
        double linearInterpolate(double x0, double y0, double x1, double y1, double x);

        // Function to find the interval in which x lies and interpolate
        double interpolate(const vector<double>& xVals, const vector<double>& yVals, double x);
    };

    class LUT
    {
    public:
        string filename;
        vector<vector<double>> data;

        LUT operator/(const LUT& other);

        LUT operator/(double val);

        LUT operator*(const LUT& other);

        LUT()
        {

        }


        LUT(string filename)
        {

            data = Parse_data(filename);

        }

    private:
        vector<vector<double>> Parse_data(const string& filename);


    };
#pragma once


}
