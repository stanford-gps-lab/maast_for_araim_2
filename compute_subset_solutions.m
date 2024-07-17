function [sigma, sigma_ss, bias, bias_ss, s1vec, s2vec, s3vec, x, chi2, So] = compute_subset_solutions(...
             G, sigpr2_int, sigpr2_acc, nom_bias_int, nom_bias_acc, subsets, pr_meas)
%*************************************************************************
%*     Copyright c 2012 The board of trustees of the Leland Stanford     *
%*                      Junior University. All rights reserved.          *
%*     This script file may be distributed and used freely, provided     *
%*     this copyright notice is always kept with it.                     *
%*                                                                       *
%*     Questions and comments should be directed to Juan Blanch at:      *
%*     blanch@stanford.edu                                               *
%*************************************************************************        
         
%COMPUTE_SUBSET_SOLUTIONS returns the standard deviation of the position error for each of the subsets
%G is the geometry matrix in ENU coordinates
%sig2pr_int (Nsat X 1) nominal variance of the pseudorange error for integrity
%sig2pr_acc (Nsat X 1) nominal variance of the pseudorange error for accuracy
%nom_bias_int (Nsat X 1) nominal bias of the pseudorange error for integrity
%nom_bias_acc (Nsat X 1) nominal bias of the pseudorange error for accurac
%subsets is an S*N matrix corresponding to all failure modes.  S is
%  arbitrary and depends on the failure assumptions.  The first line of
%  subsets is a row of ones corresponding to the fault-free mode

%Created by Juan Blanch 28/05/2012

N=size(G,1);
nsets=size(subsets,1);
W = diag(1./sigpr2_int);
GtW =G'*W;
invcov0 = GtW*G;

not_zero = find(sum(abs(invcov0))>0);   

invcov0 = invcov0(not_zero,not_zero);
G =  G(:,not_zero);
GtW = GtW(not_zero,:);
s1vec = ones(nsets,N)*Inf;
s2vec = ones(nsets,N)*Inf;
s3vec = ones(nsets,N)*Inf;
x     = zeros(nsets,size(G,2));
chi2 = ones(nsets,1)*Inf;
So = ones(length(not_zero),N)*Inf;
if cond(invcov0)<1e7
cov0    = inv(invcov0);
So = cov0*GtW;
P = W-GtW'*So;
%chi2  = ones(nsets,1)*Inf;
for i=1:nsets      
    %indb = find(subsets(i,:)==0);
    indb = (subsets(i,:)==0);
    indb_c = (subsets(i,:)==1);
    if sum(indb)<3
    Sj = So(:,indb);    
    %Pjj = P(indb,indb);
    %covindb = inv(P(indb,indb));
    %deltaCov = Sj*(covindb*(Sj')); 
    deltaCov = Sj*(P(indb,indb)\(Sj'));
    % = Sj*(P(indb,indb)*(Sj'));    
    S = 0*So;
    S(:,indb_c) = So(:,indb_c) + deltaCov*GtW(:,indb_c);

    else
    S = compute_s_coefficients(G, W, GtW, invcov0, indb, cov0);
    end
    s1vec(i,:)= S(1,:);
    s2vec(i,:)= S(2,:);
    s3vec(i,:)= S(3,:);
    if nargin==7
    x(i,:)    = (S*pr_meas)';
%     if(i==1)
%         So = S;
%     end
    chi2(i,:) = sum(((pr_meas-G*x(i,:)').^2).*subsets(i,:)'./sigpr2_int);
    end
end

%chi2 = sum(((pr_meas*ones(1,nsets)-G*x').^2).*subsets'./(sigpr2_int*ones(1,nsets)))';
end
sigma(:,1) = sqrt((s1vec.^2)* sigpr2_int);
sigma(:,2) = sqrt((s2vec.^2)* sigpr2_int);
sigma(:,3) = sqrt((s3vec.^2)* sigpr2_int);

bias(:,1)   = abs(s1vec)*nom_bias_int;
bias(:,2)   = abs(s2vec)*nom_bias_int;
bias(:,3)   = abs(s3vec)*nom_bias_int;

delta_s1vec = s1vec - ones(nsets,1)*s1vec(1,:);
delta_s2vec = s2vec - ones(nsets,1)*s2vec(1,:);
delta_s3vec = s3vec - ones(nsets,1)*s3vec(1,:);

sigma_ss(:,1) = sqrt((delta_s1vec.^2)* sigpr2_acc);
sigma_ss(:,2) = sqrt((delta_s2vec.^2)* sigpr2_acc);
sigma_ss(:,3) = sqrt((delta_s3vec.^2)* sigpr2_acc);

bias_ss(:,1)   = abs(delta_s1vec)*nom_bias_acc;
bias_ss(:,2)   = abs(delta_s2vec)*nom_bias_acc;
bias_ss(:,3)   = abs(delta_s3vec)*nom_bias_acc;




%End 