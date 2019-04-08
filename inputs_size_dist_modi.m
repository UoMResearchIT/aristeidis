%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This file contains the temperature and the                               % 
%properties of the initial size distribution                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initial temperature in K
T_i = 25.0 + 273.15;
%T_i=xlsread('MFR_data.xls','A2:A2')+273.15;
% % Thermodenuder temperatures (K)
% if ntrials > 1
% for i = 1:ntrials
% T_f(i) = Tf_range(1) + (i-1).*diff(Tf_range)./(ntrials-1);
% end
% else
% T_f = mean(Tf_range);
% end

% Monodisperse case total aerosol mass concentration (kg/m3)
%c_aer_tot_i = 13.0./1e9;
% Total aerosol number (1/m3)
%n_tot_i = 4e5.*1e6;

% Peak size (m)
%dp_peak_i = 200e-9;
rp_peak_i = dp_peak_i./2;
% Here note that the density is assumed to be 1500 kg/m3
%n_tot_i = c_aer_tot_i./(4.*1500.*rp_peak_i.^3.*pi./3);
n_tot_i = c_aer_tot_i./(4.*rho(1).*rp_peak_i.^3.*pi./3);
% Size bins
% Number of bins
nbins = 1;

% Size range (nm)
%dprange = [200e-9];

% Bins (nm)
dp_int = logspace(log10(dprange(1)),log10(dprange(end)),nbins);

% Masses in each bin (kg/m3)
c_aer_dist_int = c_aer_tot_i;
% Number in each bin (1/m3)
n_dist_int = n_tot_i;

% Total mass (kg/m3)
c_aer_int_i = nansum(c_aer_dist_int);
n_int_i = nansum(n_dist_int);

% Collecting all to one vector 
dp_i = [dp_int dp_peak_i]';
rp_i = dp_i./2;

c_aer_dist_i = [c_aer_dist_int c_aer_tot_i]';
n_dist_i = [n_dist_int n_tot_i]';








