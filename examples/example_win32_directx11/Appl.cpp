
#include "Appl.h"
#include "imgui.h"
#include "implot.h"
#include <algorithm>
#include <string.h>
/*
Length is sweeped (0.2um : 0.2um : 5um) // 25 rows
Vgs is sweeped (0.3 : 5m : 1)      // 141 column
VDS = VDD/3

*/

#define IM_ARRAYSIZE(_ARR)          ((int)(sizeof(_ARR) / sizeof(*(_ARR))))     // Size of a static C-style array. Don't use on pointers!



// operator overloading 

namespace MyApp
{
    double Interpolator::linearInterpolate(double x0, double y0, double x1, double y1, double x)
    {
        if (x1 == x0) {
            throw invalid_argument("x0 and x1 cannot be the same value");
        }
        return y0 + (y1 - y0) * ((x - x0) / (x1 - x0));
    }

    double Interpolator::interpolate(const vector<double>& xVals, const vector<double>& yVals, double x)
    {
        if (xVals.size() != yVals.size() || xVals.size() < 2) {
            throw invalid_argument("Vectors must be of the same size and contain at least two elements");
        }

        // Handle cases where x is out of range by clamping x to the nearest bounds
        if (x < xVals.front()) {
            return linearInterpolate(xVals[0], yVals[0], xVals[1], yVals[1], x);
        }
        if (x > xVals.back()) {
            return linearInterpolate(xVals[xVals.size() - 2], yVals[yVals.size() - 2], xVals.back(), yVals.back(), x);
        }

        // Find the interval [x0, x1] where x0 <= x <= x1
        for (size_t i = 0; i < xVals.size() - 1; ++i) {
            if (x >= xVals[i] && x <= xVals[i + 1]) {
                return linearInterpolate(xVals[i], yVals[i], xVals[i + 1], yVals[i + 1], x);
            }
        }

        // This should not happen due to the clamping logic above, but added as a safeguard
        throw out_of_range("x is out of the range of xVals");
    }

    LUT LUT::operator/(const LUT& other)
    {
        LUT Result;
        Result.data.resize(data.size(), vector<double>(data[0].size()));

        for (int i = 0; i < data.size(); i++)
        {
            for (int j = 0; j < data[0].size(); j++)
            {
                Result.data[i][j] = data[i][j] / other.data[i][j];

            }

        }
        return Result;
    }

    LUT LUT::operator/(double val)
    {
        LUT Result;
        Result.data.resize(data.size(), vector<double>(data[0].size()));

        for (int i = 0; i < data.size(); i++)
        {
            for (int j = 0; j < data[0].size(); j++)
            {
                Result.data[i][j] = data[i][j] / val;

            }

        }
        return Result;
    }

    LUT LUT::operator*(const LUT& other)
    {
        LUT Result;
        Result.data.resize(data.size(), vector<double>(data[0].size()));

        for (int i = 0; i < data.size(); i++)
        {
            for (int j = 0; j < data[0].size(); j++)
            {
                Result.data[i][j] = data[i][j] * other.data[i][j];

            }

        }
        return Result;
    }


    vector<vector<double>> LUT::Parse_data(const string& filename)
    {
        ifstream File;
        File.open(filename);
        if (!File.is_open()) {
            cerr << "Error opening file: " << filename << std::endl;
            return {};
        }
        string line;
        // vector<vector<double>> data;
        int row_size = 0;


        while (getline(File, line))    // Reading the file until eof is reached (getline returns false)
        {
            if (line[0] == 'M')
            {
                continue;
            }
            stringstream ss(line);     // initialize the stream with the string line (same as ss << line)
            vector<double> row;
            string value;

            while (getline(ss, value, ','))   // read each value seperated by a comma in the line 
            {
                row.push_back(stod(value));

            }
            data.push_back(row);
            row_size = row.size();

        }
        File.close(); // Close the file


        vector<double>* ptr;         // address of pointer to array
        vector<vector<double>> ID;   // rows represent lengths, columns represents VGS

        for (int L = 1; L < (row_size); L = L + 2)
        {
            vector<double> ID_L;
            for (const auto& row : data) {
                if (L < row.size()) {
                    ID_L.push_back(row[L]);
                }
            }
            ID.push_back(ID_L);
        }
        return ID;
    }


    /////////////////////// INTERPOLATION

    // Function to find the three closest points
    void findClosestPoints(const std::vector<double>& x, double x_target, int& idx0, int& idx1, int& idx2) {
        int n = x.size();
        std::vector<std::pair<double, int>> distances;

        for (int i = 0; i < n; ++i) {
            distances.push_back({ std::abs(x[i] - x_target), i });
        }

        // Sort distances based on the first element (the distance)
        std::sort(distances.begin(), distances.end());

        idx0 = distances[0].second;
        idx1 = distances[1].second;
        idx2 = distances[2].second;

        // Ensure idx0 < idx1 < idx2
        if (idx0 > idx1) std::swap(idx0, idx1);
        if (idx1 > idx2) std::swap(idx1, idx2);
        if (idx0 > idx1) std::swap(idx0, idx1);
    }

    // Function to perform second-order (quadratic) interpolation
    double quadraticInterpolation(const std::vector<double>& x, const std::vector<double>& y, double x_target) {
        if (x.size() < 3 || y.size() < 3) {
            throw std::invalid_argument("Need at least three points for second-order interpolation.");
        }

        // Find the three closest points to x_target
        int idx0, idx1, idx2;
        findClosestPoints(x, x_target, idx0, idx1, idx2);

        // Debugging: Print selected points
       // std::cout << "Using points: (" << x[idx0] << ", " << y[idx0] << "), ("
         //   << x[idx1] << ", " << y[idx1] << "), (" << x[idx2] << ", " << y[idx2] << ")\n";

        // Using points x[idx0], x[idx1], x[idx2] for interpolation
        double x0 = x[idx0], x1 = x[idx1], x2 = x[idx2];
        double y0 = y[idx0], y1 = y[idx1], y2 = y[idx2];

        // Check for distinct x values to avoid division by zero
        if ((x0 == x1) || (x1 == x2) || (x0 == x2)) {
            throw std::invalid_argument("Interpolation points must have distinct x values.");
        }

        double denom0 = (x0 - x1) * (x0 - x2);
        double denom1 = (x1 - x0) * (x1 - x2);
        double denom2 = (x2 - x0) * (x2 - x1);

        // Calculate Lagrange basis polynomials
        double L0 = ((x_target - x1) * (x_target - x2)) / denom0;
        double L1 = ((x_target - x0) * (x_target - x2)) / denom1;
        double L2 = ((x_target - x0) * (x_target - x1)) / denom2;

        // Calculate the interpolated value
        double y_target = y0 * L0 + y1 * L1 + y2 * L2;

        return y_target;
    }

    // Function to perform quadratic interpolation for given points
    double quadraticInterpolationp(double x0, double y0, double x1, double y1, double x2, double y2, double x) {
        // Calculate the differences
        double denom0 = (x0 - x1) * (x0 - x2);
        double denom1 = (x1 - x0) * (x1 - x2);
        double denom2 = (x2 - x0) * (x2 - x1);

        if (denom0 == 0 || denom1 == 0 || denom2 == 0) {
            throw std::invalid_argument("Interpolation points must have distinct x values.");
        }

        // Calculate Lagrange basis polynomials
        double L0 = ((x - x1) * (x - x2)) / denom0;
        double L1 = ((x - x0) * (x - x2)) / denom1;
        double L2 = ((x - x0) * (x - x1)) / denom2;

        // Print intermediate values for debugging

        // Calculate the interpolated value
        double y = y0 * L0 + y1 * L1 + y2 * L2;

        return y;
    }



    void func()
    {
      
        static bool opt_fullscreen = true;
        static bool opt_padding = false;
        static ImGuiDockNodeFlags dockspace_flags = ImGuiDockNodeFlags_None;

        // We are using the ImGuiWindowFlags_NoDocking flag to make the parent window not dockable into,
        // because it would be confusing to have two docking targets within each others.
        ImGuiWindowFlags window_flags = ImGuiWindowFlags_MenuBar | ImGuiWindowFlags_NoDocking;
        if (opt_fullscreen)
        {
            const ImGuiViewport* viewport = ImGui::GetMainViewport();
            ImGui::SetNextWindowPos(viewport->WorkPos);
            ImGui::SetNextWindowSize(viewport->WorkSize);
            ImGui::SetNextWindowViewport(viewport->ID);
            ImGui::PushStyleVar(ImGuiStyleVar_WindowRounding, 0.0f);
            ImGui::PushStyleVar(ImGuiStyleVar_WindowBorderSize, 0.0f);
            window_flags |= ImGuiWindowFlags_NoTitleBar | ImGuiWindowFlags_NoCollapse | ImGuiWindowFlags_NoResize | ImGuiWindowFlags_NoMove;
            window_flags |= ImGuiWindowFlags_NoBringToFrontOnFocus | ImGuiWindowFlags_NoNavFocus;
        }
        else
        {
            dockspace_flags &= ~ImGuiDockNodeFlags_PassthruCentralNode;
        }

        // When using ImGuiDockNodeFlags_PassthruCentralNode, DockSpace() will render our background
        // and handle the pass-thru hole, so we ask Begin() to not render a background.
        if (dockspace_flags & ImGuiDockNodeFlags_PassthruCentralNode)
            window_flags |= ImGuiWindowFlags_NoBackground;

        // Important: note that we proceed even if Begin() returns false (aka window is collapsed).
        // This is because we want to keep our DockSpace() active. If a DockSpace() is inactive,
        // all active windows docked into it will lose their parent and become undocked.
        // We cannot preserve the docking relationship between an active window and an inactive docking, otherwise
        // any change of dockspace/settings would lead to windows being stuck in limbo and never being visible.
        if (!opt_padding)
            ImGui::PushStyleVar(ImGuiStyleVar_WindowPadding, ImVec2(0.0f, 0.0f));
        ImGui::Begin("DockSpace Demo", nullptr, window_flags);
        if (!opt_padding)
            ImGui::PopStyleVar();

        if (opt_fullscreen)
            ImGui::PopStyleVar(2);

        // Submit the DockSpace
        ImGuiIO& io = ImGui::GetIO();
        if (io.ConfigFlags & ImGuiConfigFlags_DockingEnable)
        {
            ImGuiID dockspace_id = ImGui::GetID("MyDockSpace");
            ImGui::DockSpace(dockspace_id, ImVec2(0.0f, 0.0f), dockspace_flags);
        }


        if (ImGui::BeginMenuBar())
        {
            if (ImGui::BeginMenu("Options"))
            {
                // Disabling fullscreen would allow the window to be moved to the front of other windows,
                // which we can't undo at the moment without finer window depth/z control.
                ImGui::MenuItem("Fullscreen", NULL, &opt_fullscreen);
                ImGui::MenuItem("Padding", NULL, &opt_padding);
                ImGui::Separator();

                if (ImGui::MenuItem("Flag: NoDockingOverCentralNode", "", (dockspace_flags & ImGuiDockNodeFlags_NoDockingOverCentralNode) != 0)) { dockspace_flags ^= ImGuiDockNodeFlags_NoDockingOverCentralNode; }
                if (ImGui::MenuItem("Flag: NoDockingSplit", "", (dockspace_flags & ImGuiDockNodeFlags_NoDockingSplit) != 0)) { dockspace_flags ^= ImGuiDockNodeFlags_NoDockingSplit; }
                if (ImGui::MenuItem("Flag: NoUndocking", "", (dockspace_flags & ImGuiDockNodeFlags_NoUndocking) != 0)) { dockspace_flags ^= ImGuiDockNodeFlags_NoUndocking; }
                if (ImGui::MenuItem("Flag: NoResize", "", (dockspace_flags & ImGuiDockNodeFlags_NoResize) != 0)) { dockspace_flags ^= ImGuiDockNodeFlags_NoResize; }
                if (ImGui::MenuItem("Flag: AutoHideTabBar", "", (dockspace_flags & ImGuiDockNodeFlags_AutoHideTabBar) != 0)) { dockspace_flags ^= ImGuiDockNodeFlags_AutoHideTabBar; }
                if (ImGui::MenuItem("Flag: PassthruCentralNode", "", (dockspace_flags & ImGuiDockNodeFlags_PassthruCentralNode) != 0, opt_fullscreen)) { dockspace_flags ^= ImGuiDockNodeFlags_PassthruCentralNode; }
                ImGui::Separator();


                ImGui::EndMenu();
            }


            ImGui::EndMenuBar();
        }


        // READ LUTS
        LUT ID_n("Nmos_ID.csv");
        LUT GM_n("Nmos_GM.csv");
        LUT GDS_n("NMOS_GDS.csv");
        LUT CGG_n("NMOS_CGG.csv");
        LUT GMOVERID_n("NMOS_GMOVERID.csv");
        LUT VDSAT_n("NMOS_VDSAT.csv");
        double W = 10e-6;

        LUT ID_p("Pmos_ID.csv");
        LUT GM_p("Pmos_GM.csv");
        LUT GDS_p("PMOS_GDS.csv");
        LUT CGG_p("PMOS_CGG.csv");
        LUT GMOVERID_p("PMOS_GMOVERID.csv");
        LUT VDSAT_p("PMOS_VDSAT.csv");

        LUT VGS("VGS.csv");



        // MORE LUTS
        LUT GM_GDS_n = GM_n / GDS_n;
        LUT ID_W_n = ID_n / W;
        LUT GM_GDS_p = GM_p / GDS_p;
        ID_p = ID_p / -1;
        LUT ID_W_p = ID_p / W;
        VDSAT_p = VDSAT_p / -1;

        // row0 ID - row 1 GM - row 2 GDS - row 3 CGG and so on 
        LUT LUTs[9][2] = { {GMOVERID_n,GMOVERID_p},{VGS,VGS},{VDSAT_n,VDSAT_p}, {ID_n,ID_p} ,{GM_n,GM_p} , {GDS_n,GDS_p} ,{CGG_n,CGG_p}
                        , {GM_GDS_n , GM_GDS_p} , {ID_W_n,ID_W_p} , };



        // Choose the Device Type 
        const char* device[] = { "NMOS", "PMOS" };
        static int device_current = 0;
        ImGui::Combo("Device Type", &device_current, device, IM_ARRAYSIZE(device));

        // Choose X-axis
        const char* x_axis[] = { "GM/ID", "VGS" , "VDSAT" };
        static int x_current = 0;
        ImGui::Combo("X-axis", &x_current, x_axis, IM_ARRAYSIZE(x_axis));

        // Choose Y-axis
        const char* y_axis[] = { "GM/ID", "VGS" , "VDSAT","ID", "GM" , "GDS" , "CGG","GM/GDS", "ID/W" };
        static int y_current = 0;
        ImGui::Combo("Y-axis", &y_current, y_axis, IM_ARRAYSIZE(y_axis));

        double x_s[141], y_s[141];
        Interpolator intrp;


        static int clicked = 0;
        if (ImGui::Button("PLOT!"))
        {
            clicked++;
        }

        if (clicked & 1)
        {

            double achieved_target;
            ImGui::Begin("Working Area");

            static float slider_L = 0.2;
            float slider_Lmax = 5.0f;
            float slider_Lmin = 0.2f;

            // try to get the y-vector for a in-between L
            vector <double> y_for_L(141, 0.0);

            vector <vector<double>> temp_y(141, vector<double>(25, 0.0));
            vector <vector<double>> temp_xx(141, vector<double>(25, 0.0));

            vector <double> temp_x(25, 0.0);
            for (int i = 0; i < 25; i++)
            {
                temp_x[i] = (0.2 * 0.000001) + (i * 0.2 * 0.000001);
            }

            for (int i = 0; i < 141; i++)
            {
                for (int j = 0; j < 25; j++)      // j is the length index  
                {
                    temp_y[i][j] = LUTs[y_current][device_current].data[j][i];  // temp_y[i] is the y-val for all lengths at certain x-val
                    temp_xx[i][j] = LUTs[x_current][device_current].data[j][i];
                }
            }

            vector <double> x_for_L(141);
            for (int i = 0; i < 141; i++)
            {

                try {
                    // All is good but again gm/gds are giving off values  
                //    y_for_L[i] = quadraticInterpolationp(x_0 * 0.000001, LUTs[y_current][device_current].data[length_index][i]
                  //      , x_1 * 0.000001, LUTs[y_current][device_current].data[length_index + 1][i]
                    //    , x_2 * 0.000001, LUTs[y_current][device_current].data[length_index + 2][i], slider_L * 0.000001);
                  //  y_for_L[i] = quadraticInterpolation(temp_x, temp_y[i], slider_L * 0.000001);

                    y_for_L[i] = quadraticInterpolation(temp_x, temp_y[i], slider_L * 0.000001);
                    x_for_L[i] = quadraticInterpolation(temp_x, temp_xx[i], slider_L * 0.000001);

                }
                catch (const std::exception& e)
                {
                    std::cerr << "Error: " << e.what() << endl;
                }

            }

            // max_element returns an iterator to the first occurence of the max
            ImGui::SliderFloat("Length", &slider_L, slider_Lmin, slider_Lmax, "%.7e");

            int length_index = floor((slider_L - 0.2f) / 0.2f);

            auto slider_xmax = max_element(x_for_L.begin(), x_for_L.end());
            auto slider_xmin = min_element(x_for_L.begin(), x_for_L.end());
            // max_element returns an iterator to the first occurence of the max
            static float slider_x = *slider_xmin;
            ImGui::SliderFloat(x_axis[x_current], &slider_x, *slider_xmin, *slider_xmax, "%.7e");


            auto slider_ymax = max_element(y_for_L.begin(), y_for_L.end());
            auto slider_ymin = min_element(y_for_L.begin(), y_for_L.end());
            // max_element returns an iterator to the first occurence of the max
            static float slider_y = *slider_ymin;
            ImGui::SliderFloat(y_axis[y_current], &slider_y, *slider_ymin, *slider_ymax, "%.7e");



            try {
                // achieved_target = quadraticInterpolation(LUTs[x_current][device_current].data[length_index],
                 //    LUTs[y_current][device_current].data[length_index], slider_x);

                achieved_target = quadraticInterpolation(x_for_L,y_for_L, slider_x);
            }
            catch (const std::exception& e)
            {
                std::cerr << "Error: " << e.what() << endl;
            }


            slider_y = achieved_target;

            string len = to_string(slider_L);
            char* length_c = new char[len.size() + 1];
            strcpy(length_c, len.c_str());

            char print[20];
            strcat(print, "L= ");
            strcat(print, length_c);
            strcat(print, "um");


            string slid_X = to_string(slider_x);
            char* slid_x_c = new char[slid_X.size() + 1];
            strcpy(slid_x_c, slid_X.c_str());

            string slid_Y = to_string(slider_y);
            char* slid_y_c = new char[slid_Y.size() + 1];
            strcpy(slid_y_c, slid_Y.c_str());

            char slidx[50];
            strcat(slidx, slid_x_c);
            strcat(slidx, ",");
            strcat(slidx, slid_y_c);


            for (int i = 0; i < 141; i++)
            {
                // x_s[i] = LUTs[x_current][device_current].data[length_index][i];
                 //y_s[i] = LUTs[y_current][device_current].data[length_index][i];
                x_s[i] = x_for_L[i];
                y_s[i] = y_for_L[i];

            }
            if (ImPlot::BeginPlot("Plot", x_axis[x_current], y_axis[y_current])) {
                ImPlot::PlotLine( print , x_s, y_s, 141);
                ImVec4 color = ImVec4(0.9f, 0.1f, 0.1f, 1.0f); // Red color
                ImVec2 offset = ImVec2(10, 10); // Offset the text slightly
                ImPlot::Annotation(slider_x, slider_y,color, offset,true,slidx);
                ImPlot::EndPlot();

            }
            ImGui::End();




        }

        ImGui::End();
        //ImGui::ShowDemoWindow();




    }    // func 
}      // main
