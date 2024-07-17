function subsets_k = determine_k_subsets(n,k)
%*************************************************************************
%*     Copyright c 2013 The board of trustees of the Leland Stanford     *
%*                      Junior University. All rights reserved.          *
%*     This script file may be distributed and used freely, provided     *
%*     this copyright notice is always kept with it.                     *
%*                                                                       *
%*     Questions and comments should be directed to Juan Blanch at:      *
%*     blanch@stanford.edu                                               *
%*************************************************************************

%Created 14 August 2013 by Juan Blanch

%DETERMINE_K_SUBSETS determines all the subsets of size k out of n
%The output is a matrix where each line corresponds to a subset.
%If subsets_k(i,j) = 1, satellite j is in subset i, otherwise
%subsets_k(i,j)=0;

if k==0
    subsets_k = zeros(1,n);
end

if k==1
    subsets_k = eye(n);
end

if k>=2
    subsets_k = zeros(nchoosek(n,k),n);
    index =1;
    for i=1:n-k+1
        subsets_n_i_k_1 = determine_k_subsets(n-i, k-1);
        subsets_i =[ zeros(size(subsets_n_i_k_1,1),i-1) ones(size(subsets_n_i_k_1,1),1) subsets_n_i_k_1];
        subsets_k(index:(index+size(subsets_i,1)-1),:) = subsets_i;
        index =index+size(subsets_i,1); 
    end
end


%End




