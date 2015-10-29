%% Test the coordinates system transform betweem Cartesian and Spherical coordinate
% The spherical coordinate display the tangential direction
% Spherical Coordinates
Theta=pi/2;
Phi=pi/4;
R=1;
% Cartesian Coordinates
X=R*sin(Theta)*cos(Phi);
Y=R*sin(Theta)*sin(Phi);
Z=R*cos(Theta);
% Base Vector in Cartesian Coordinates
Ele_Vector_CC=[X,0,0;0,Y,0;0,0,Z];
% Base Vector in Spherical Coordinates
CC_to_SC=...
    [sin(Theta)*cos(Phi),sin(Theta)*sin(Phi),cos(Theta);
    cos(Theta)*cos(Phi),cos(Theta)*sin(Phi),-sin(Theta);
    -sin(Phi),cos(Phi),0;];
Ele_Vector_SC=CC_to_SC*Ele_Vector_CC;
% Test Coordinates Transformation (Spherical Coordinates to Cartesian Coordinates)
SC_to_CC=...
    [sin(Theta)*cos(Phi),cos(Theta)*cos(Phi),-sin(Phi);
    sin(Theta)*sin(Phi),cos(Theta)*sin(Phi),cos(Phi);
    cos(Theta),-sin(Phi),0;];

Test_Ele_Vector_CC=SC_to_CC*Ele_Vector_SC;
% Test Coordinates Transformation (Cartesian Coordinates to Spherical Coordinates)
Test_Ele_Vector_SC=CC_to_SC*Test_Ele_Vector_CC;

