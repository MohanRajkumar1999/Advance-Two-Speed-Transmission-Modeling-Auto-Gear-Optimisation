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

%Velocity at which gear ratio shifting occurs 
V_shift = 25; %Kmph
Teff = 0.89;

%%Motor
%%Motor Efficiency
data = readmatrix("Nissan Leaf Motor Efficiency.xls");
MMeff_Torque = data(:,17);
MMeff_Speed = data(:,1);
MMeff = data(:,2:16);

%Motor Characteristic Graph
MCG = readmatrix("Nissan Leaf MC.csv");
M_CG_S = MCG(:,1); % Motor Speed
M_CG_PT = MCG(:,2); % Motor Peak Torque
M_CG_CT = MCG(:,3); % Motor Continuous Torque
%Motor efficiency with less point
MLD = readmatrix("Motor Efficiency Data Less Points.csv");
MLD_MS = MLD(:,1);
MLD_MT = MLD(:,2);
MLD_ME = MLD(:,3);
%Coverting number to text
MLD_ME_text = string(MLD_ME*100);

%Optimization of two speed transmission
G_1 = linspace(G1min,G1max,3);
G_2 = linspace(G2min,G2max,3);
[G1_range,G2_range] = meshgrid(G_1,G_2);

%Motor Peak Values
MPT = 280; %Nm % Motor Peak Torque
MPS = 10390; %Rpm % Motor Peak Speed

%Empty variable to store data
Motor_eff = [];
M_Energy = [];
Motor_avg_eff = [];
Motor_Torque= [];
Motor_Speed = [];
M_Energy_mat = [];

%Calculating the no of interations

iteration = 1;
iteration_sel = 1;
irange = [];
k = 1;
%%Selected gear ratio
GS1 = 0;
GS2 = 0;

%Defining the motor energy
Motor_E_sel = inf;
Motor_avg_eff_sel = 0;

% %Defing the for loop
for i = 1:length(G_1)
    for j = 1:length(G_2)
        %Selecting the gear ratio values for the test
        G1 = G1_range(i,j);
        G2 = G2_range(i,j);
        
        %Simulating the simulink two speed transmission model 
        out = sim('Gear_Ratio_Sizing_Two_Speed.slx');
        
        %Storing the test data output
        Motor_avg_eff(i,j) = mean(out.MMo_eff.Data);
        M_Energy(iteration) = out.ME.Data(end);
        M_Energy_mat(i,j) = out.ME.Data(end);
        Motor_Torque(:,iteration)= out.M_Torque.Data;
        Motor_Speed(:,iteration)= out.Wheel_Speed.Data;
        
        %%Motor Operating Region
        
        figure('Name', 'Motor Charateristic graph')
        plot(M_CG_S,M_CG_PT,'r')
        text(1000,280,'Motor Peak Torque')
        hold on
        plot(M_CG_S,M_CG_CT,'g')
        text(1000,155,'Motor Continuous Torque')
        %%Motor Operating region
        scatter(out.M_Speed.data,out.M_Torque_GR1.Data,'b')
        scatter(out.M_Speed.data,out.M_Torque_GR2.Data,'r')
        legend('Gear Ratio 1','Gear Ratio 2')
        text(MLD_MS,MLD_MT,MLD_ME_text)
        %Motor Efficiency point
        scatter(MLD_MS,MLD_MT,MLD_ME,'r','LineWidth',5);
        xlabel('Motor Speed (rpm)')
        ylabel ('Motor Torque (Nm)')
        
        %Saving the figure as a image
        filename = ("Motor operating points" + " Gear ratio 1= " + string(G1)+" Gear ratio 2= "+ string(G2)+ ".png");
        saveas(gcf(),filename)
        figure('Name', 'Motor Efficiency Map')
        %%Motor Efficiency map
        surf(MMeff_Speed,MMeff_Torque,MMeff)
        hold on
        %Motor efficiency points wrt model
        scatter3(out.M_Speed.data,out.M_Torque.data,out.MMo_eff.data)
        xlabel('Motor Speed (rpm)')
        ylabel ('Motor Torque (Nm)')
        zlabel ('Motor Efficiency (%)')
        filename = ("Motor Efficiency Map" + " Gear ratio 1= " + string(G1)+" Gear ratio 2= "+ string(G2)+ ".png");
        saveas(gcf(),filename)
        
        %Verifying if the operating points are within the motor
        %Characteristics graph
        if (max(out.M_Torque.data) < MPT && max(out.M_Speed.data) < MPS)
            %identify the gear ratio which has the lowest energy
            %Cosumption the highest efficiency
            if (M_Energy_mat(i,j) < Motor_E_sel && Motor_avg_eff(i,j) > Motor_avg_eff_sel)
                %%Selection of energy consumed 
                Motor_E_sel = M_Energy_mat(i,j);
                Motor_avg_eff_sel = Motor_avg_eff(i,j);
                GS_1 = G1;
                GS_2 = G2;
                row = i;
                col = j;
                Motor_Torque_sel_G1 =  out.M_Torque_GR1.Data;
                Motor_Torque_sel_G2 =  out.M_Torque_GR2.Data;
                Motor_Speed_sel =  out.M_Speed.data;
                Motor_Eff_val_sel = out.MMo_eff.Data;
                %Storing the energy consumption values wrt each iteration
                irange(iteration) = iteration;
                io = iteration;
                Motor_Energy_e(iteration) = M_Energy_mat(i,j);
                Motor_Eff_e(iteration) = Motor_avg_eff(iteration)
            end
        end
        iteration = iteration + 1;
    end
     
end
iteration = iteration - 1;
%Display the selected gear ratio value
disp('Total number of iteration')
disp(iteration)
disp('Optimal gear ratio 1')
disp(GS_1)
disp('Optimal gear ratio 2')
disp(GS_2)
disp('Motor Energy Consumed (wh)')
disp(Motor_E_sel)
disp('Motor average efficiency(%)')
disp(Motor_avg_eff_sel)
disp('Total number of iteration')
disp(iteration)
%Motor Energy consumption based on gear ratio value
figure('Name','Motor Energy based on different gear ratio')
surf(G1_range,G2_range,M_Energy_mat)
xlabel('Gear ratio 1')
ylabel('Gear ratio 2')
zlabel('Motor Energy Consumed (Wh)')
title('Motor energy based on different gear ratio value')
filename = ("Motor energy based on different gear ratio value" + ".png");
saveas(gcf(),filename)

%Optimal motor energy consumed based on the algorithm 
figure('Name','Optimal motor energy consumed based on the algorithm')
scatter(io,Motor_E_sel,'r')
xlabel('Coverging interation');
ylabel('Optimal motor energy consumption (wh)')
title('Coverging interation to obtain the optimal energy consumption')

%Optimal motor average efficiency based on the algorithm 
figure('Name','Optimal motor average efficiency based on the algorithm ')
scatter(io,Motor_avg_eff_sel,'r')
xlabel('Coverging interation');
ylabel('Optimal motor average efficiency (%)')
title('Coverging interation to obtain the optimal motor efficiency')

%Conergence of the solution for each and every iteration
figure('Name','Motor Energy for each convergence of the solution')
plot(irange,Motor_Energy_e,'-o')
hold on
scatter(io,Motor_E_sel,'r')
text(io,Motor_E_sel+1,'Optimal motor energy consumption')
xlim([-1 max(irange)+10])
ylim([-100 max(Motor_Energy_e)+200])
xlabel('Interation at which solution converged');
ylabel('Optimal motor average efficiency (%)')
title('Motor energy based on the algorithm convergence')

%Conergence of the solution for each and every iteration
figure('Name','Motor average Effiency for each convergence of the solution')
plot(irange,Motor_Eff_e,'-o')
hold on
scatter(io,Motor_avg_eff_sel+1,'r')
text(io,Motor_avg_eff_sel,'Optimal motor average efficiency consumption')
xlim([-1 max(irange)+5])
ylim([-1 max(Motor_Eff_e)+10])
xlabel('Interation at which solution converged');
ylabel('Optimal motor average efficiency (%)')
title('Motor average energy based on the algorithm convergence')

%Plotting
%Motor Characteristic graph
figure('Name', 'Optimal Motor Charateristic graph')
plot(M_CG_S,M_CG_PT,'r')
text(1000,280,'Motor Peak Torque')
hold on
plot(M_CG_S,M_CG_CT,'g')
text(1000,155,'Motor Continuous Torque')
%%Motor Operating region
scatter(Motor_Speed_sel,Motor_Torque_sel_G1,'b')
scatter(Motor_Speed_sel,Motor_Torque_sel_G2,'r')
legend('Gear Ratio 1','Gear Ratio 2')
text(MLD_MS,MLD_MT,MLD_ME_text)
%Motor Efficiency point
scatter(MLD_MS,MLD_MT,MLD_ME,'r','LineWidth',5);
xlabel('Motor Speed (rpm)')
ylabel ('Motor Torque (Nm)')


figure('Name', 'Optimal Motor Efficiency Map')
%%Motor Efficiency map
surf(MMeff_Speed,MMeff_Torque,MMeff)
hold on
%Motor efficiency points wrt model
scatter3(Motor_Speed_sel,Motor_Torque_sel_G1+Motor_Torque_sel_G2,Motor_Eff_val_sel)
xlabel('Motor Speed (rpm)')
ylabel ('Motor Torque (Nm)')
zlabel ('Motor Efficiency (%)')
