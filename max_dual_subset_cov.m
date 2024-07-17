function [sigma, sigma_ss, bias, bias_ss, x, chi2] = max_dual_subset_cov(G,sigpr2_int,nom_bias_int, nom_bias_acc, subset)

%*************************************************************************
%*     Copyright c 2024 The board of trustees of the Leland Stanford     *
%*                      Junior University. All rights reserved.          *
%*     This script file may be distributed and used freely, provided     *
%*     this copyright notice is always kept with it.                     *
%*                                                                       *
%*     Questions and comments should be directed to Juan Blanch at:      *
%*     blanch@stanford.edu                                               *
%*************************************************************************
%max_dual_subset_cov returns a bound on the 2-out subset sigmas
%It is based on the general formula described in 
%Blanch, J., and Walter, T., "Fast Protection Levels for Fault Detection With an Application to Advanced RAIM," 
%in IEEE Transactions on Aerospace and Electronic Systems, vol. 57, no. 1, pp. 55-65, Feb. 2021, doi: 10.1109/TAES.2020.3011997.
%Created Juan Blanch June 6, 2024

if nargin<5
    subset =1:size(G,1);
end
%sigma = determine_nminusk_sigmas(G,sigpr2_int,1,subset);

%remove columns from unobservable clocks
idob = ~(sum(abs(G),1)==0);
idob(1:3)=1;
G = G(:,idob);


N = length(subset);

%close all
W = diag(1./sigpr2_int);
cov0 = inv(G'*W*G);

S = cov0*G'*W;

P = W-W*G*S;
S = S(:,subset);
P =P(subset,subset);
P_ii= diag(P);
Pnorm = P./(sqrt(P_ii)*sqrt(P_ii'));

s3 = S(3,:);
snorm3 = s3'./sqrt(P_ii);
sigma2_3_0 = cov0(3,3);

s2 = S(2,:);
snorm2 = s2'./sqrt(P_ii);
sigma2_2_0 = cov0(2,2);

s1 = S(1,:);
snorm1 = s1'./sqrt(P_ii);
sigma2_1_0 = cov0(1,1);

[snorm3_max2,~] = maxk(snorm3.^2,2);
[snorm2_max2,~] = maxk(snorm2.^2,2);
[snorm1_max2,~] = maxk(snorm1.^2,2);

snorm3_sum = sum(snorm3_max2);
snorm2_sum = sum(snorm2_max2);
snorm1_sum = sum(snorm1_max2);

%Pnorm_max = max(max(abs(Pnorm-eye(N))));

pnormmax = max(max(abs(Pnorm-eye(N))));

sigma_ss2 = [snorm1_sum snorm2_sum snorm3_sum]/(1-pnormmax) ;

sigma2 = [sigma2_1_0 sigma2_2_0 sigma2_3_0] + sigma_ss2;

sigma = sqrt(sigma2);

sigma_ss = sqrt(sigma_ss2);

b = sqrt((nom_bias_int.^2)'*(1./sigpr2_int));

b_acc = sqrt((nom_bias_acc.^2)'*(1./sigpr2_int));

bias = sigma*b;

bias_ss = sigma_ss*b_acc;

x = zeros(1,size(G,2));

chi2 = 0;


