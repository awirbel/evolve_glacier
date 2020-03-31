// mesh resolution 
cl1 = 20;

// mesh corners
north = 1000;
south = 0;
west = 1000;
east = 0;

// surface points domain
Point(1) = {west, south, 0, cl1};
Point(2) = {east, south, 0, cl1};
Point(3) = {east, north, 0, cl1};
Point(4) = {west, north, 0, cl1};
// make the lines
Line(1) = {1, 2};
Line(2) = {4, 3};
Line(3) = {1, 4};
Line(4) = {2, 3};
//+
Line Loop(1) = {3, 2, -4, -1};
//+
Plane Surface(1) = {1};

// now build the extrude the surface
Extrude {0, 0, -1000} {
  Surface{1}; Layers{1};
}
