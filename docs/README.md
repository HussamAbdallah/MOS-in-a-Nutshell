MOS in a NutShell
=====
"Mos in a NutShell" works as a Lookup table for different mosfet parameters!
The program uses data extracted from Spectre simulator for NMOS and PMOS, where length is sweeped as following (L= 0.2um : 0.2um : 5um)
VDS = VDD/3 
neglecting VSB,
W = 10um
The extracted parameters are saved as comma seperated values "csv", to operate on the values a parsing mechanism is needed. 
The "LUT" class contains a private parsing method "Parse_data(const string& filename)" the filename is passed as a parameter 
it returns a 2D vector where each row represent the values at a certain Length.

To have a MOS paramter at any length (L') given parameter values at (L= 0.2um : 0.2um : 5um) at a certain x a quadratic interpolation is applied at the required length (L').
This is done for all x-values (To get a 2-D graph of the paramter at certain L).
Then, a second quadratic interpolation is applied to get the parameter value at a certain x!

![GM/GDS vs GM/ID](https://i.postimg.cc/50BrTT50/image.png)

One Feature of the program is that it gives the valid range (Min-Max) of each paramter that can be controlled by a slider so that the designer can 
pick the sutiable paramters values to achieve requirements.
![](https://i.postimg.cc/FzWwVqXB/Animation.gif)

The Program uses "Dear ImGui" [ImGui](https://github.com/ocornut/imgui) for the graphical user interface.
