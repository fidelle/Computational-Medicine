function y_diff = caffeine_cost(k0, t_exp, y_exp, tspan,tspan2, y0, p, weights)

% AMOXICILLIN_COST_2C y_diff = amoxicillin_cost_2c(k0, data, tspan, y0, p)
%   Sets up the cost function for the optimization of the absorption and
%   elimination rate constants from the provided patient data. Calls on
%   amoxicllin_main_2c for the ODE solver; used by lsqnonlin in the driver
%   file
%
%   INPUT:
%       k0 = initial guesses for the parameter(s) to be optimized
%       t_exp = time points (in hours) corresponding to patient data
%       y_exp = patient concentration data, input as a matrix where rows -
%       different time points and columns = different dosing methods
%       tspan = vector of times to use for the beginning and end of the
%       system, given in hours
%       y0 = initial values for the model, input as a matrix where the rows
%       are the equations and the columns are the different dosing methods
%       p = structure containing parameters for the model, input as a
%       structure
%   OUTPUT:
%       ydiff = matrix of the difference between the simulation and the
%       experimental data at each time point; rows = different time points
%       and columns = different dosing methods; lsqnonlin will minimize the
%       sum of squares of this matrix

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Model Solver



% Initialize output matrix
Y = [];

% For each patient 
for i = 1:5
    p.ka=k0(1,i);
    p.kc=k0(2,i);
    p.v=k0(3,i);
%    
% %    % p.ka=k0(1,i);
%     p.kc=k0(1,i);
%     p.v=k0(2,i);
% 
    % Run the ODE solver
    [T1i, Y1i, ~, ~] = caffeine_main(tspan,y0(:,i),p);
    
    y0new=Y1i(end,:);
    y0new(end,3)=y0new(end,3)+92;
    
    
    [T2i, Y2i, ~, ~] = caffeine_main(tspan2,y0new,p);
    
    Ti=[T1i;T2i];
    Yi=[Y1i;Y2i];
    
    % Select the output from the central compartment
    Yi = Yi(:,1);
    
    % Store output
    Y = cat(2,Y,Yi);
    
end

% Find difference between model output and patient data
% Initialize y_diff matrix
y_diff = zeros(length(t_exp),5);

% For each experimental time point
for k = 1:length(t_exp)
    
    % Difference between the output and experimental time
    t_diff = abs(Ti - t_exp(k));
    
    % Finding the value of T that minimizes the difference
    [~, t_index] = min(t_diff);
    
    % Find the difference between the experimental concentrations and model
    % output at that time point
    y_diff(k,:) = Y(t_index, :) - y_exp(k,:);
    
    
end















