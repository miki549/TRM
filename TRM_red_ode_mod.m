function dy = TRM_red_ode_mod(t, y, K, cap)
%ODE model of generalized ribosome flow models with mass action kinetics for simulation
%reduced state space: we only consider particles (vehicles) but no space variables
%INPUTS
%K: interconnection matrix, where K(i,j) is the transition rate from compartment i to compartment j
%K(i,j)=0: no transition from comp. i to j
%cap: capacity vector of the compartments

%state vector: y = [n(1) ... n(m)]


m=size(K,1); %number of compartments
%c=40*ones(m,1); %capacities of compartments
n=y(1:m); %vehicles (particles)


%assemble monomial matrix
NS=zeros(m,m);
for i=1:m
    for j=1:m
        NS(i,j)=n(i)*(cap(j)-n(j));
    end
end

tmp_dy=zeros(m, 1); %for storing the derivatives

%constructing the automous part of the equations
for i=1:m
    NS_row = NS(i,:);
    NS_col = NS(:,i);
    K_row = K(i,:);
    K_col = K(:,i);
    outflow = K_row * NS_row';
    inflow = K_col' * NS_col;
    tmp = inflow - outflow;
    tmp_dy(i) = tmp; %derivative of particles in the i-th compartment 
end


%add external inflows and outflows using the rates defined
ext_flows = TRM_external_flows(t,m);
ext_inflow = ext_flows(:,1).*(cap-y);
ext_outflow = ext_flows(:,2).*y;

tmp_dy = tmp_dy + ext_inflow - ext_outflow;

dy = tmp_dy;



