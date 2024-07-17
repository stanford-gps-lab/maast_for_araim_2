function [vpl, hpl, sig_acc, emt, subsets, pfault_inst, pfault_exp, p_not_monitored, alloc] = mhss_raim_baseline_v5(G, sigpr2_int, sigpr2_acc, nom_bias_int, nom_bias_acc, ...
    p_sat, p_const,tau_sat, tau_const, opt_flag, prm, rho_j, subsets, pfault_inst,pfault_exp, p_not_monitored)

%*************************************************************************
%*     Copyright c 2024 The board of trustees of the Leland Stanford     *
%*                      Junior University. All rights reserved.          *
%*     This script file may be distributed and used freely, provided     *
%*     this copyright notice is always kept with it.                     *
%*                                                                       *
%*     Questions and comments should be directed to Juan Blanch at:      *
%*     blanch@stanford.edu                                               *
%*************************************************************************



%Created 14 August 2012 by Juan Blanch
%Modified 17 March 2015 by Juan Blanch (refinement of exclusion list)
%Modified 24 August 2017 by Juan Blanch (allocation to exclusion, modified
%Q function, and subset selection algorithm)
%Updated April 2024 by Juan Blanch
%MHSS_RAIM_BASELINE_V5 implements the vpl, hpl, sigma accuracy, and emt
%computation described in the ARAIM ADD v4.2 as well as speed execution improvements
% (based on Blanch, J., Lee, Y., Walter, T., Enge, P., Pervan, B., Belabbas, B., Spletter, A., Rippl, M.,
%"Advanced RAIM User Algorithm Description: Integrity Support Message Processing,
%Fault Detection, Exclusion, and Protection Level Calculation," Proceedings of the 25th 
%International Technical Meeting of The Satellite Division of the Institute
%of Navigation (ION GNSS 2012), Nashville, September 2012.)

%G (Nsat X (3+Nconst)) geometry matrix in ENU. G(i,3+j)=1 if satellite i belongs to constellation j and zero otherwise
%sig2pr_int (Nsat X 1) nominal variance of the pseudorange error for integrity
%sig2pr_acc (Nsat X 1) nominal variance of the pseudorange error for accuracy
%nom_bias_int (Nsat X 1) nominal bias of the pseudorange error for integrity
%nom_bias_acc (Nsat X 1) nominal bias of the pseudorange error for accuracy
%p_sat (Nsat X 1) a priori probability of satellite fault
%p_const (Nsat X 1) a priori probability of constellation fault
%if p_sat or p_const are one scalar, it is assumed that it applies to all satellites or constellations
%rho_j is the fraction of the integrity budget given to exclusion mode j.
%For FD, rho_j =1


if nargin<11
    rho_j = 1;
end

%%%%%%%%%%%% Determine subsets and associated probabilities %%%%%%%%%%%%%%%

if nargin<13
[subsets, pfault_inst, pfault_exp, p_not_monitored]= determine_subsets_v5(G, p_sat, p_const,tau_sat, tau_const, prm.p_thres, prm.fc_thres, prm.p_exc_thres);
alloc = pfault_exp*0;
end

%Find unique subsets and consolidate probabilities

[subsets_unique, pfault_inst_unique, pfault_exp_unique] = find_unique_subsets(subsets, pfault_inst, pfault_exp);

%vector indicating how many subsets are represented by a row
nmodes = ones(length(pfault_inst_unique),1);



%%%%%%%%%%%%% Compute subset position solutions, sigmas, and biases %%%%%%%
if~isempty(subsets)
       two_out = (sum(1-subsets_unique,2)==2);
    if (prm.fast&&~opt_flag&&sum(two_out)>0) 
       %bound all two-out subsets instead of computing them explicitly
        two_out = (sum(1-subsets_unique,2)==2);
        all_other =~two_out;  
        pfault_inst_o = pfault_inst_unique(all_other);
        pfault_exp_o = pfault_exp_unique(all_other);
        nmodes_o = nmodes(all_other);

        [sigma_o, sigma_ss_o, bias_o, bias_ss_o, s1vec, s2vec, s3vec, ~, chi2_o] = compute_subset_solutions(...
             G, sigpr2_int, sigpr2_acc, nom_bias_int, nom_bias_acc, subsets_unique(all_other,:), zeros(size(G,1),1));
         
        pfault_inst_two_out =sum(pfault_inst_unique(two_out));       
        pfault_exp_two_out =sum(pfault_exp_unique(two_out));
        nmodes_two_out = sum(nmodes(two_out));

        [sigma_two_out, sigma_ss_two_out, bias_two_out, bias_ss_two_out, ~, chi2_two_out] = max_dual_subset_cov(G,sigpr2_int,nom_bias_int, nom_bias_acc);

        
        
        sigma = [sigma_o;sigma_two_out];
        sigma_ss = [sigma_ss_o;sigma_ss_two_out];
        bias = [bias_o;bias_two_out];
        bias_ss = [bias_ss_o;bias_ss_two_out];
        pfault_inst_for_pl = [ pfault_inst_o;pfault_inst_two_out];
        pfault_exp_for_pl = [ pfault_exp_o;pfault_exp_two_out];
        nmodes = [nmodes_o;nmodes_two_out];
        chi2 = [chi2_o;chi2_two_out];
        x = 0*chi2;


    else
     [sigma, sigma_ss, bias, bias_ss, s1vec, s2vec, s3vec, x, chi2] = compute_subset_solutions(...
             G, sigpr2_int, sigpr2_acc, nom_bias_int, nom_bias_acc, subsets_unique, zeros(size(G,1),1));
         
    pfault_inst_for_pl = pfault_inst_unique;
    pfault_exp_for_pl = pfault_exp_unique;
    end
else
    vpl = Inf;
    hpl = Inf;
    sig_acc = Inf;
    emt = Inf;
    return
end


%%%%%%%%%%%%% Adjust all-in-view position coefficients %%%%%%%%%%%%%%%%%%%%

if opt_flag
    nsets=size(subsets_unique,1);

    [ s3opt ] = compute_adjusted_position_1D(pfault_inst_unique,sigma(:,3),sigma_ss(:,3), s3vec, sigpr2_acc,prm.sig_acc_max_vert);
    [ s2opt ] = compute_adjusted_position_1D(pfault_inst_unique,sigma(:,2),sigma_ss(:,2), s2vec, sigpr2_acc,prm.sig_acc_max_hor2);
    [ s1opt ] = compute_adjusted_position_1D(pfault_inst_unique,sigma(:,1),sigma_ss(:,1), s1vec, sigpr2_acc,prm.sig_acc_max_hor1);

    s1vec(1,:) = s1opt;  
    s2vec(1,:) = s2opt;  
    s3vec(1,:) = s3opt;

    sigma(1,1) = sqrt((s1vec(1,:).^2)* sigpr2_int);
    sigma(1,2) = sqrt((s2vec(1,:).^2)* sigpr2_int);
    sigma(1,3) = sqrt((s3vec(1,:).^2)* sigpr2_int);

    bias(1,1)   = abs(s1vec(1,:))*nom_bias_int;
    bias(1,2)   = abs(s2vec(1,:))*nom_bias_int;
    bias(1,3)   = abs(s3vec(1,:))*nom_bias_int;

    delta_s1vec = s1vec - ones(nsets,1)*s1vec(1,:);
    delta_s2vec = s2vec - ones(nsets,1)*s2vec(1,:);
    delta_s3vec = s3vec - ones(nsets,1)*s3vec(1,:);

    sigma_ss(:,1) = sqrt((delta_s1vec.^2)* sigpr2_acc);
    sigma_ss(:,2) = sqrt((delta_s2vec.^2)* sigpr2_acc);
    sigma_ss(:,3) = sqrt((delta_s3vec.^2)* sigpr2_acc);

    bias_ss(:,1)   = abs(delta_s1vec)*nom_bias_acc;
    bias_ss(:,2)   = abs(delta_s2vec)*nom_bias_acc;
    bias_ss(:,3)   = abs(delta_s3vec)*nom_bias_acc;
end

%%%%%%%%%%%%%%%%%%%% Compute sigma accuracy %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sigma3_acc = sqrt((s3vec(1,:).^2)*sigpr2_acc);


%%%%%%%%%%%%% Filter out modes that cannot be monitored and adjust%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%% integrity budget %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin<13
     [sigma, sigma_ss, bias, bias_ss,...
     pfault_inst_for_pl, pfault_exp_for_pl, p_not_monitored, ~] = filter_out_subsets(sigma, sigma_ss, bias,...
     bias_ss, pfault_inst_for_pl, pfault_exp_for_pl, p_not_monitored, x, chi2,nmodes);
end

%%%%%%%%%%%%%%%%%%%%%%% Compute test thresholds %%%%%%%%%%%%%%%%%%%%%%%%%%%
[T1, T2, T3] = compute_test_thresholds(sigma_ss, bias_ss, prm.pfa_vert/prm.N_es_cont, prm.pfa_hor/prm.N_es_cont,nmodes);


%%%%%%%%%%%%%%%%%%%%%% Allocate integrity budget across vertical and horizontal %%%%%%%%%%%%%%%%%%%%%% %%%
pfault_inst_for_pl(1) = 2;
if prm.max_pmd_flag == 1
     allocmax =min((pfault_exp_for_pl./pfault_inst_for_pl)./prm.N_es_int,1);
else
     allocmax =min((pfault_exp_for_pl./pfault_inst_for_pl),1);
end
allocmax(1)=1;
phmi_vert = rho_j*(prm.phmi_vert - prm.phmi_vert/(prm.phmi_vert+prm.phmi_hor)*p_not_monitored);
phmi_hor  = rho_j*(.5*prm.phmi_hor -.5*prm.phmi_hor/(prm.phmi_vert+prm.phmi_hor)*p_not_monitored);


%%%%%%%%%%%%%%%%%%%%% Compute Protection Levels %%%%%%%%%%%%%%%%%%%

if prm.vpl_target~=Inf
   if prm.fast
       [vpl]    = compute_protection_level_v5(sigma(:,3), bias(:,3) + T3, pfault_inst_for_pl,  phmi_vert/prm.N_es_int, prm.pl_tol,allocmax,prm.vpl_target);
   end
   if ~prm.fast||vpl>prm.vpl_target
       [vpl,  ~]    = compute_protection_level(sigma(:,3), bias(:,3) + T3, pfault_inst_for_pl,  phmi_vert/prm.N_es_int, prm.pl_tol,allocmax);
   end
else
   vpl = 1e6;
end

if prm.hpl_variant==0
    
    if prm.fast
    [hpl1]    = compute_protection_level_v5(sigma(:,1), bias(:,1) + T1, pfault_inst_for_pl,  phmi_hor/prm.N_es_int, prm.pl_tol, allocmax, prm.hpl_target/sqrt(2));
    [hpl2]    = compute_protection_level_v5(sigma(:,2), bias(:,2) + T2, pfault_inst_for_pl,  phmi_hor/prm.N_es_int, prm.pl_tol, allocmax, prm.hpl_target/sqrt(2));
     hpl = sqrt(hpl1^2 + hpl2^2);
    end
    if (~prm.fast)||hpl>prm.hpl_target
        [hpl1, ~]    = compute_protection_level(sigma(:,1), bias(:,1) + T1, pfault_inst_for_pl,  phmi_hor/prm.N_es_int, prm.pl_tol, allocmax);
        [hpl2, ~]    = compute_protection_level(sigma(:,2), bias(:,2) + T2, pfault_inst_for_pl,  phmi_hor/prm.N_es_int, prm.pl_tol, allocmax);   
        hpl = sqrt(hpl1^2 + hpl2^2);
    end
end

%%%%% Alternate HPL computations  %%%%
if prm.hpl_variant==1
%%%% Blanch ION GNSS 2023 %%%%
pfault_inst_for_pl(1)=1;
      if prm.fast       
        pfault_inst_for_pl(1)=1;
        [hpl]    = compute_hor_protection_level_v5(sigma(:,1),sigma(:,2), bias(:,1) + T1,...
                       bias(:,2) + T2, pfault_inst_for_pl,  phmi_hor/prm.N_es_int, prm.pl_tol, allocmax,prm.hpl_target);
      end
      if (~prm.fast)||hpl>prm.hpl_target
       [hpl,~]    = compute_hor_protection_level(sigma(:,1),sigma(:,2), bias(:,1) + T1,...
                       bias(:,2) + T2, pfault_inst_for_pl,  phmi_hor/prm.N_es_int, prm.pl_tol, allocmax);
      end
end


if prm.hpl_variant==2

%%%% Racelis et al ION GNSS 2022 with factor of two correction %%%
    T_2D = sqrt((bias(:,1) + T1).^2 + (bias(:,2) + T2).^2);
    sigma_2D = sqrt(sigma(:,1).^2 + sigma(:,2).^2);
    pfault_inst_for_pl(1)=1;
    if prm.fast
    [hpl]    = compute_protection_level_v5(sigma_2D, T_2D, 2*pfault_inst_for_pl,  phmi_hor/prm.N_es_int, prm.pl_tol, allocmax,prm.hpl_target);
    end
    if (~prm.fast)||hpl>prm.hpl_target
         [hpl, ~]    = compute_protection_level(sigma_2D, T_2D, 2*pfault_inst_for_pl,  phmi_hor/prm.N_es_int, prm.pl_tol, allocmax);
    end

end


%%%%%%%%%%%%%%%%%%% Compute Effective Monitor threshold %%%%%%%%%%%%%%%%%%%

emt = compute_emt(T3, pfault_inst_unique, prm.p_emt);

if ~isempty(sigma3_acc)
    sig_acc = sigma3_acc(1);
else
    sig_acc = Inf;
    emt = Inf;
end
%End

