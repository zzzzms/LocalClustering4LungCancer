function Cluster = LsqSingleClusterPursuit(A,Gamma,n0,epsilon,t,reject)
% function [C,v] = LsqSingleClusterPursuit(A,Gamma,n0,epsilon,t,reject)


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

