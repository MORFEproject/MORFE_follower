
// global dimensions

L=30; // length of beam
H=2; // thick. of beams

// mesh parameters
 
lc1=2;

// point coordinates

Point(1) = {0,0,0,lc1};
Point(2) = {L,0,0,lc1};
Point(3) = {L,H,0,lc1};
Point(4) = {0,H,0,lc1};

Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,1};
Line Loop(5) = {1,2,3,4};

// surface

Plane Surface(1) = {5};

// physical entities

Physical Line(1) = {4}; // left clamp
Physical Line(2) = {1};   // force lower side
Physical Line(3) = {2};   // force tip
Physical Surface(4) = {1}; 

// creates second order mesh and saves

Mesh.MshFileVersion=2;
Mesh.ElementOrder=2;
Mesh 2;
Save 'cantilever.msh';
