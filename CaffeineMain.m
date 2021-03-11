function [T,Y,Balance,AUC] = caffeine_main(tspan,y0,p)

% AMOXICILLIN_MAIN_2C [T,Y,Balance,AUC] = amoxicillin_main_2c(tspan,y0,p)
%   Solves the system of ODEs given in the amoxicillin_eqns_2c file based
%   on the parameters given. Inputs are the length of time, initial
%   conditions, and parameters; outputs are the time, concentrations, mass
%   balance, and area under the curve from the ODE solver. Test code for 
%   Systems Pharmacology and Personalized Medicine
%
%   INPUT:
%       tspan = vector of times to use for the beginning and end of the
%       system, given in hours
%       y0 = initial values for the model, input as a vector where the
%       positions correspond to the equations
%       p = structure containing parameters for the model, input as a
%       structure
%   OUTPUT:
%       T = vector of times corresponding to the concentrations in Y, given
%       in seconds
%       Y = matrix of concentration values output from the ODEs; each
%       column corresponds to an equation in order, each row is a time
%       point
%       Balance = vector containing the mass balance for the ODE solution,
%       each row is a time point
%       AUC = the area under the time versus concentration curve in the 
%       central compartment

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Model Solver

% Options for ODE Solver
options = odeset('MaxStep',5e-2, 'AbsTol', 1e-5,'RelTol', 1e-5,...
    'InitialStep', 1e-2);

% Run ODE Solver
[T,Y] = ode45(@Caffeine_eqns, tspan, y0, options, p);

% Mass Balance
TotalD(:,1) = Y(:,1)*p.v;           % Drug in blood stream (mg)
TotalD(:,2) = Y(:,3);           % Drug to be absorbed (mg)
TotalD(:,3) = Y(:,2);           % Drug removed from system (mg)

DrugIn = p.q*T + p.dose;      % Total drug input to system (mg)
DrugOut = TotalD(:,3);              % Total drug removed from system (mg)
Balance = DrugIn - DrugOut - TotalD(:,1) - TotalD(:,2);

if any(Balance > 1e-6)
    fprintf('Mass is not balanced: %g\n', Balance)
end

% % AUC
AUC = trapz(T,Y(:,1));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




















