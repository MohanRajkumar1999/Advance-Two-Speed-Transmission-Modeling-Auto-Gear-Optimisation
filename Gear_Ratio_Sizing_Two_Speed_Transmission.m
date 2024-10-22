%Importing the drive cycle
Drive = readmatrix('UDDS.csv');

%Vehicle
%Vehicle mass
Mv = 1650.665; %kg
%Mass of drive
Md = 80;
%Gross vehicle mass
GVM = Mv+Md;
%Gravitational constant
g = 9.81;
%Gross vehicle weight
GVW = GVM*g;

%Vehicle Resistive
%Coefficient of rolling resistance
Crf = 0.015;
%Frontal area 
Fa = 3.8056;
%Air denisty
rho = 1.225;
%Drag coefficient
Cd = 0.28;

%Radius of the wheel
Rw = 0.2032;

%Two Speed Transmission
G1max = 11.075;
G1min = 6.06;
G2max = 7.96;
G2min = 5.41;

G1_range = [G1min,(G1min+G1max)/2,G1max];
G2_range = [G2min,(G2min+G2max)/2,G2max];

%Velocity at which gear ratio shifting occurs 
V_shift = 25; %Kmph
Teff = 0.89;

%%Motor
%%Motor Efficiency
data = readmatrix("Nissan Leaf Motor Efficiency.xls");
MMeff_Torque = data(:,17);
MMeff_Speed = data(:,1);
MMeff = data(:,2:16);

