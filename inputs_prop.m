%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This file contains the properties of the evaporating compounds          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Mixture of sub-cooled succinic acid adipic acid

% Molar gas constant (J/mol/K)
R = 8.314472;

% Number of species
nspec = length(cstar);

% Molar masses of the compounds [kg/mol]
MW(1:length(cstar))=0.2;

% Volatilities of the surrogate compounds at T_ref (K). The volatilities can be
% given in either vapor pressures or saturation concentrations. 
T_ref = 298.15;
% % Vapor pressure in Pascals
% pstar = []; %
% Saturation concentration in ug/m3
%cstar = pstar.*MW.*1e9./R./T_ref 
%cstar = [0.01 0.1 1.0 10.0];
% Volatilities of the surrogate compounds at T_ref [Pa]
pstar = cstar*R.*T_ref./MW./1e9;

% Diffusion coefficients of the surrogate compounds in air at T_ref [m2/s]
% and the temperature-dependent factor
Dn(1:length(cstar)) =1.0e-5; 
mu(1:length(cstar)) =1.75;

% Densities of the surrogate compounds [kg/m3] 
%rho(1:length(cstar)) =1500; % 

% Surface tensions of the surrogate compounds [N/m]
sigma(1:length(cstar))=0.05; %

% Vaporization enthalpies of the surrogate compounds at 298 K [J/mol]
%dHvap = 1e3.*[189.3134  169.5311  150.8227  133.1875  116.6248  101.1333 ]; % 
%dHvap = 1e3.*(85 - 11.*log10(cstar));
%dHvap=[160000 160000 160000 160000];

% Mass and thermal accommodation coefficients of the surrogate species
%alpha_m = [0.8 0.8 0.8 0.8];
alpha_t(1:length(cstar)) =1.0;

% Initial mass fractions of the species in the aerosol
%X_i =[0.5 0.0 0.5 0.0];