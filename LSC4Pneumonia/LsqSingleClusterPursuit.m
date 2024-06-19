function Cluster = LsqSingleClusterPursuit(A,Gamma,n0,epsilon,t,reject)
% function [C,v] = LsqSingleClusterPursuit(A,Gamma,n0,epsilon,t,reject)

% ========================= Acknowledgement ============================
% This code is based on the code of SingleClusterPursuit algorithm by
% Dr. Daniel Mckenzie, with the parameter 'reject' in Mckenzie's 
% SingleClusterPursuit is removed. It was written by Zhaiming Shen under 
% Dr. Ming-Jun Lai's supervision in 2022. 
% ======================================================================

% INPUT
% =================================
% A .................... Adjacency matrix of data converted to graph form. 
% Gamma ................ VECTOR. Labelled data within cluster of interest 
% n0 ................... (estimated) size of C_a
% epsilon .............. Omega_a will be of size (1+epsilon)n0(a)
% t .................... Depth of random walk 
%
% OUTPUT
% ================================
% Cluster...... VECTOR. The elements in the cluster of interest.
%

% ========= Initialization ================= %
n = size(A,1); % number of vertices
degvec = sum(A,2);
Dinv = spdiags(1./degvec,0,n,n);
DinvA = Dinv*A;
L = speye(n,n) - DinvA;

% ============ Call the subroutines ============= %
Omega = RandomWalkThresh(A,Gamma,n0,epsilon,t);
[Cluster,~] = LsqClusterPursuit(L,Gamma,Omega,n0,reject);
end

