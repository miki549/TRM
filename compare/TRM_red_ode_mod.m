function dy = TRM_red_ode_mod(t, y, K, cap, k_stoch_func)
%ODE model of generalized ribosome flow models with mass action kinetics for simulation
%reduced state space: we only consider particles (vehicles) but no space variables
%INPUTS
%K: base interconnection matrix (constant)
%cap: capacity vector of the compartments
%reaction_matrix: reaction pairs
%k_stoch_func: time-dependent rate function

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

%constructing the automous part of the equations with time-dependent rates
for i=1:m
    outflow = 0;
    inflow = 0;
    
    % Calculate outflows from compartment i
    for j=1:m
        if K(i,j) > 0
            % Use time-dependent rate function instead of constant K
            rate = k_stoch_func(t, i, j);
            outflow = outflow + rate * n(i) * (cap(j) - n(j));
        end
    end
    
    % Calculate inflows to compartment i
    for j=1:m
        if K(j,i) > 0
            % Use time-dependent rate function instead of constant K
            rate = k_stoch_func(t, j, i);
            inflow = inflow + rate * n(j) * (cap(i) - n(i));
        end
    end
    
    tmp_dy(i) = inflow - outflow; %derivative of particles in the i-th compartment 
end

%add external inflows and outflows using the rates defined
ext_flows = TRM_external_flows(t,m);
ext_inflow = ext_flows(:,1).*(cap-y);
ext_outflow = ext_flows(:,2).*y;

tmp_dy = tmp_dy + ext_inflow - ext_outflow;

dy = tmp_dy;



