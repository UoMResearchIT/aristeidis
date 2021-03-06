%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is the main program to model the dynamic evaporation in the TD     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% close all
% clear all

% Loading the inputs:
% Residence time, TD length and constants
inputs_TD
[resultsX, qq] = cstar_combinations(cstar, step);

% Properties of the evaporating compounds
inputs_prop
% Size distribution
inputs_size_dist_modi

% Preparation for parfor loop
qi = 0;
len_index = qq * length(dH) * length(alp);
[qw_index, qk_index, qz_index] = deal(zeros(1, len_index));
qi_index = 1:len_index;

for qw=1:qq
    for qk=1:length(dH)
         for qz=1:length(alp)
             qi = qi+1;
             qw_index(qi) = qw;
             qk_index(qi) = qk;
             qz_index(qi) = qz;
         end
    end
end

n_cstar = length(cstar);

parfor qi = qi_index
    qw = qw_index(qi);
    qk = qk_index(qi);
    qz = qz_index(qi);
    
    X_i = resultsX(qw, :);
        
    %Initial mole fractions of the species in the aerosol
    n_i_apu = X_i./MW;
    n_i_tot_apu = sum(n_i_apu);
    Xm_i = n_i_apu./n_i_tot_apu;
    % Initial density of the aerosol
    rhol_i = sum(X_i.*rho); % Mass-weighted average
    % Initial surface tension of the aerosol
    sigmal_i = sum(Xm_i.*sigma); % Mole-weighted average
    
    dHvap = repmat(dH(qk), 1, n_cstar);
    % Saturation pressures at initial temperature
    psat_i = pstar.*exp(dHvap.*(1./T_ref - 1./T_i)./R);
    
    alpha_m = repmat(alp(qz), 1, n_cstar);

    %  Kelvin effect corresponding to the initial composition
    Ke_i = exp(2.0 .* MW .* sigmal_i ./(R * T_i .* rho .* rp_i));
    % Mass fraction based equilibrium pressures
    peq_i = X_i .* psat_i .* Ke_i;
   
    % Initial partial pressures of the species, assuming aerosol initially in
    % equilibrium with the aerosol corresponding to the peak size
    pv_i = peq_i(end, :);
    
    % Initial particle mass
    mp_i = rhol_i.*4.0*pi.*rp_i.^3.0./3.0;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Calculation for each TD temperature                                     %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    time = linspace(0,t_res_heat,10000);
    dt = mean(diff(time));

    % Preallocate arrays
    output = zeros(length(time), nbins+3, nspec);
    Pc_k_end = zeros(ntrials, length(1:nbins+1), size(output, 3));
    Gc_k_end = zeros(ntrials, length(nbins+2:nbins+3), size(output, 3));
    [~, ncols, ~] = size(Pc_k_end);
    mp_end = zeros(ntrials, ncols);
    ncols_end = size(mp_end, 2);
    X_end = zeros(1, ncols_end, nspec);
    apu = zeros(nspec, 1);
    rhol_end = zeros(ntrials, nbins+1);
    [rp_end, dp_end, c_aer_dist_end] = deal(zeros(ntrials, ncols_end));
    [c_aer_tot_end, mfr_dist, mfr_mono] = deal(zeros(1, ncols_end));

    for k = 1:ntrials
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Initializing the time-dependent variables at the TD temp for calculation %
        %of evaporation in the heating section                                    %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % Time
        time0 = 0;
        % Space
        x_coord0 = l_heat.*time0./t_res_heat;

        % if t_res_heat <17.0 && t_res_heat >13.0
        % % Temperature profile for 16 s residence time in the beginning of the
        % % reactor
        %  para0 = [2.4509 -8.0918 5.0610 0.1405];
        %  T_nd0 = para0(1).*x_coord0.^3 + para0(2).*x_coord0.^2 + para0(3).*x_coord0 + para0(4);
        %  T_TD = T_nd0.*(T_f(k) - T_i) + T_i;
        % else
        % % Flat temperature profile
        T_TD = T_f(k);
        assert(max(size(T_TD)) == 1)
        % end

        % Equilibrium pressures at the TD temperature
        psat = pstar.*exp(dHvap.*(1./T_ref - 1./T_TD)./R);
        csat = MW.*psat./R./T_TD;

        [Ke, peq, pv_a] = deal(zeros(nbins + 1, nspec));
        for i = 1:nspec
            % Kelvin effect corresponding to the initial composition
            Ke(1:nbins+1,i) = exp(2.0.*MW(i).*sigmal_i./R./T_TD./rho(i)./rp_i);
            % Mole fraction based equilibrium pressures
            %peq(1:nbins+1,i) = Xm_i(i).*psat(i).*Ke(1:nbins+1,i);
            % Mass fraction based equilibrium pressures
            peq(1:nbins+1,i) = X_i(i).*psat(i).*Ke(1:nbins+1,i);
            % Partial pressure at the particle surface
            pv_a(1:nbins+1,i) = peq(1:nbins+1,i);
        end

        % Total particle mass concentration (kg/m3)
        c_aer_tot = c_aer_tot_i.*T_i./T_TD;
        % Particle mass concentrations (kg/m3)
        c_aer_dist = c_aer_dist_i.*T_i./T_TD;
        % Particle number concentration (1/m3)
        n_tot = n_tot_i.*T_i./T_TD;
        % Particle number concentrations in each bin (1/m3)
        n_dist = n_dist_i.*T_i./T_TD;

        % Diffusion coefficient (m2/s)
        D = Dn.*(T_TD./T_ref).^mu;

        % Particle mass (kg) and size (m)
        mp = mp_i;
        rp = rp_i;

        % Vapor mass concentrations (kg/m3)
        cgas = pv_i.*MW./R./T_TD;

        % Total concentrations of the species
        ctot = cgas + X_i.*c_aer_dist(end);

        % Mean velocity of the gas molecules:
        c_ave = sqrt(8.*R.*T_TD./MW./pi);
        % Mean free path of the gas molecules:
        lambda = 3.*D./c_ave;

        [Pc, Kn, beta] = deal(zeros(nbins + 1, nspec));
        for i = 1:nspec
            % Masses of each species in the individual particles (kg)
            Pc(1:nbins+1,i) = X_i(i).*mp;
            % Knudsen number:
            Kn(1:nbins+1,i) = lambda(i)./rp_i;
            % Fuchs and Sutugin transition regime correction
            beta(1:nbins+1,i) = (1.0 + Kn(1:nbins+1,i))./(1.0 + (4./alpha_m(i)./3 + 0.377).*Kn(1:nbins+1,i) + 4.*Kn(1:nbins+1,i).^2./alpha_m(i)./3);
        end

        % Gas phase concentrations kg/m3 for the size distribution case and for the
        % monodisperse case

        Gc = [cgas; cgas];

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %     Calculating the time-dependent evaporation                          %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        input = [Pc; Gc];                

        % Converting the input matrix to column vector for the odesolver
        input0 = reshape(input, [], 1);

        % Solving the mass fluxes
        options = odeset('RelTol',1E-4,'AbsTol',1E-13);
        [tout, output0] = ode45(@fluxes_size_dist, time, input0,  options, dt,nbins,nspec,l_heat,t_res_heat,...
            T_f(k),T_i,n_tot_i,n_dist_i,pstar,dHvap,T_ref,MW,sigma,rho,Dn,mu,p,alpha_m,alpha_t);

        % % Converting back to our format
        output = reshape(output0, size(output, 1), size(output, 2), size(output, 3));
        
        % Temporary evolution of particle and gas phase masses (kg) and (kg/m3)
        Pc_t_k = output(:, 1:nbins+1, :);
        Gc_t_k = output(:, nbins+2:nbins+3, :);

        findex1 = find(Pc_t_k <= 0.0);
        Pc_t_k(findex1) = 0.0;
        findex2 = find(Gc_t_k <= 0.0);
        Gc_t_k(findex2) = 0.0;

        % Particle and gas phase masses at the end of the heating section (kg) and
        % (kg/m3)
        Pc_k_end(k, :, :) = Pc_t_k(end, :, :);
        Gc_k_end(k, :, :) = Gc_t_k(end, :, :);

        % Particle masses in the end (kg)
        mp_end(k, :) = sum(Pc_k_end(k, :, :), 3);

        findex3 = find(mp_end <=0.0);
        mp_end(findex3) = 0.0;

        % Particle composition in the end
        for i = 1:nspec
            X_end(k, :, i) = Pc_k_end(k, :, i)./mp_end(k, :);
        end

        findex4 = find(X_end <= 0.0 | isnan(X_end));
        X_end(findex4) = 0.0;

        % Particle density (kg/m3)
        for j = 1:nbins+1
            for i = 1:nspec
                apu(i) = X_end(k,j,i).*rho(i);
            end
            rhol_end(k,j) = sum(apu);
        end

        findex5 = find(rhol_end <= 0.0);
        rhol_end(findex5) = 1000;
        % Particle radii and diameters
        rp_end = (3 * mp_end ./ (4 * pi .* rhol_end)).^(1/3);
        dp_end = rp_end * 2;

        % Final concentrations corrected back to T_i
        % The number concentrations in each bin are the same as initially
        c_aer_dist_end(k, :) = n_dist_i(:)'.*mp_end(k, :);
        % Total mass
        c_aer_tot_end(k) = sum(c_aer_dist_end(k,1:nbins));

        % Mass fraction remaining
        mfr_dist(k) = c_aer_tot_end(k)./c_aer_int_i;
        mfr_mono(k) = c_aer_dist_end(k,nbins+1)./n_tot_i./mp_i(end);
    end

    % different combinations of properties examined
    X_a(qi,:)=X_i;
    enthalpia(qi)=dHvap(1);
    alpha(qi)=alpha_m(1);
    
    for qd=1:ntrials
        MFR(qi,qd)=mfr_dist(qd);
    end
    % calculating the error
    error(qi)=sqrt(sum((experimental-mfr_dist).^2))/ntrials;
end
% all the combinations of properties with the corresponding error
results = [X_a, enthalpia', alpha', error', MFR];
