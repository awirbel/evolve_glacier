// grid resolution for outer and inner mesh
cl1 = 20;
cl2 = 20;

// mesh corners that are based on the geotiff edges
north = 1400.0;
south = 600.0;
west = 600.0;
east = 1400.0;

// here I move one half dx/dy inside the domain to avoid the edges vs the center coordinates
Point(1) = {west, south, 0, cl1};
Point(2) = {east, south, 0, cl1};
Point(3) = {east, north, 0, cl1};
Point(4) = {west, north, 0, cl1};
// make the transfinite lines
Line(1) = {1, 2};
Line(2) = {4, 3};
Line(3) = {1, 4};
Line(4) = {2, 3};

Line Loop(6) = {3, 2, -4, -1};

Plane Surface(8) = {6};
//+

