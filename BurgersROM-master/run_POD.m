function [POD] = run_POD(snapshots,r,M)
% function [Phi,omega,lambda,b,Xdmd] = DMD(snapshots,r,dt)
% Computes the Dynamic Mode Decomposition of X1, X2
%
% INPUTS: 
% snapshots = data matrix
% r = target rank of SVD
% dt = time step advancing X1 to X2 (X to X')
% M = FEM mass matrix
%
% OUTPUTS:
% POD, the POD modes


%% POD
% [U,S,V] = svd(snapshots,'econ');
% POD = U(:,1:r); % POD modes

%-------------------------------------------------------------------------
%  Weight the snapshot data
%-------------------------------------------------------------------------

[Q,R] = mgs_weighted(snapshots,M);

[U,Sigma,V] = svd(R,0);
% U(:,1:r) is V_r = W_r
POD = Q*U(:,1:r);

% fprintf('Generated POD basis functions\n')

end
