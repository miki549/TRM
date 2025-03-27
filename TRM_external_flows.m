function out = TRM_external_flows(t,N)
%function for defining external in and outflows for the networked TRM model
%used by both 'TRM_red_ode_sim' and 'TRM_obsv'
%INPUTS:
%t: time
%N: number of compartments
%
%OUTPUT:
%out first column: inflows, second column: outflows
%out(i,1): external inflow rate to compartment i
%out(i,2): outflow rate from compartment i

flow=zeros(N,2);

%inflow rates
flow(1,1) = 1*(sin(t/5)+1) + 2*(sin(pi*t/6)+1);
%flow(2,1) = 2*(cos(sqrt(2)*t/8)+1) + (sin(t/3)+1);

%outflow rates
flow(2,2) = 0.4;


out=flow;




