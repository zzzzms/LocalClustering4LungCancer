function Cluster = main_CS_LCE(A,Gamma,n0,epsilon,t,reject)

% ================================================================= %
% This program implements the LCE Algorithm.

% ========================= Acknowledgement =============================
% It is modified based on Dr. Daniel Mckenzie's original code. 
% Zhaiming Shen. April 2023
% =======================================================================


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
[Cluster,~] = CS_LCE(L,Gamma,Omega,n0,reject);
end