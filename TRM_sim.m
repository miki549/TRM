clear all
close all

%physical parameters of continuous time model where the state variable is vehice density [vehicle / m]:
V_max = 5; %maximum speed [m/s]
l = 5; %average vehicle length including minimum space between vehicles [m]
rho_max = 1/l; %maximum vehicle density [vehicle / m]
omega = V_max/rho_max; %[m^2 / (vehicle * s)]
delta_x=200; %cell size [m]
k = omega/delta_x; %transition rate between cells [m / (vehicle * s)]
cm = delta_x/l; %maximum capacity of cells [vehicle]

%reaction (transition) rate of models where the state variable is the number of vehicles (spaces) in each cell:
%this is simply obtained by transforming the state variables in the kinetic ODEs from densities to vehicle numbers
%needed to compute directly in vehicle / space units (otherwise not important)
k_stoch = k / delta_x;

Nc = 6; %number of compartments
l=5; %vehicle length
cc = (delta_x / l) * ones(Nc,1); %capacities of compartments

K=zeros(Nc,Nc);
kij=k;
K(1,2) = kij;
K(2,3) = kij;
K(3,4) = kij;
K(4,5) = kij;
K(5,6) = kij;
K(6,1) = kij;

K(3,6) = kij;
K(1,5) = kij;
K(4,2) = kij;

%K=rand(Nc,Nc)

K=K/delta_x; %state variables are molecule (vehicle) numbers (not concentrations)

cap=40*ones(Nc,1); %capacities of compartments

n0 = [38; 38; 30; 2; 2; 2];
y0 = [n0];

tfinal=400
tspan = [0 tfinal]; %time is in seconds

options = odeset('RelTol',1e-8,'AbsTol',1e-8);

[t_det, y_det]=ode45(@(t,y) TRM_red_ode_mod(t,y,K,cap), tspan, y0, options);

figure

plot(t_det, y_det, 'LineWidth', 2)

newcolors = [0.83 0.14 0.14; ...
       	     1.00 0.54 0.00; ...
             0.2 0.1 1; ...
             0.25 0.80 0.54; ...
             0.47 0.25 0.80; ...
             0.5 0 0.5]
       
colororder(newcolors)

legend('n_1','n_2','n_3','n_4','n_5','n_6')
title('TRM simulation - Traffic states')
ylabel('states')
xlabel('time')
colororder(newcolors)
grid





