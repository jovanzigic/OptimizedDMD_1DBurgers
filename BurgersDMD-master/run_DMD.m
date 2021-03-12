function [relerrDMD_r,relerr1_r,ROM_error_DMD,ROM_error_optDMD,tOD,tDMD] = run_DMD(X,r,dt,M,t)
%%
% The DMD method provides spatiotemporal decomposition of data into a set
% of dynamic modes that are derived from snapshots of a system in time.
%
% INPUTS: 
% snapshots = data matrix
% r = target rank of SVD
% dt = time step advancing X1 to X2 (X to X')
% M = FEM mass matrix
%
% OUTPUTS:
% Phi, the DMD modes
% omega, the continuous-time DMD eigenvalues
% lambda, the discrete-time DMD eigenvalues
% b, a vector of magnitudes of modes Phi
% Xdmd, the data matrix reconstructed by Phi, omega, b

%% Step 1: 
% 'm' state snapshots from k = 1:m are arranged into two large matrices:
% X1 = snapshots from k = 1:m-1
% X2 = snapshots from k = 2:m

tic
% Columns of X1 and X2 are state snapshots 
X1 = X(:, 1:end-1);
X2 = X(:, 2:end);

% %% Step 2: Weight snapshot data by mass matrix
% 
% % R1 and R2 are the weighted snapshots
% [Q1,R1] = mgs_weighted(X1,M);
% [Q2,R2] = mgs_weighted(X2,M);

%% Step 3a: Compute the DMD modes
% the DMD modes are the eigenvectors of A

[U, S, V] = svd(X1, 'econ'); % SVD of X1
r = min(r, size(U,2)); % choose rank of reduced model

% truncate to rank-r
U_r = U(:, 1:r);
S_r = S(1:r, 1:r);
V_r = V(:, 1:r);

% A_r is r x r projection of A onto POD modes
A_r = U_r' * X2 * V_r / S_r; % low-rank dynamics
[W_r, D] = eig(A_r); % eigendecomposition of A_r
Phi = X2 * V_r / S_r * W_r; % DMD modes

% DMD spectra
lambda = diag(D); % discrete-time eigenvalues
omega = log(lambda)/dt; % continuous-time eigenvalues

%% Step 4: Compute DMD mode amplitudes b

x1 = X1(:, 1); % initial state snapshot
b = Phi\x1; % x1 = Phi*b, b = the initial amplitudes of each mode

%% Step 5: DMD reconstruction
% Xdmd = x(t), the approximate solution of future states

time_dynamics = zeros(r, length(X1));
t_vec = (0:length(X1)-1)*dt; % time vector
for i = 1:length(X1)
    time_dynamics(:,i) = (b.*exp(omega*t_vec(i)));
end
Xdmd = Phi * time_dynamics;
tDMD = toc;

relerrDMD_r = norm(Xdmd-X1,'fro')/norm(X1,'fro');
% % fprintf('Relative Error of Standard DMD NLLS Problem is %g\n',relerrDMD_r);


%% Step 3b, 4b, 5b: Compute Optimized DMD modes, amplitudes and reconstruct X

% fit to unprojected data
tic
imode = 1;
[w1,e1,b1] = optdmd(X,t,r,imode);

% reconstructed values
Xd1 = w1*diag(b1)*exp(e1*t);
relerr1_r = norm(Xd1-X,'fro')/norm(X,'fro');
tOD = toc;

% fprintf('Relative Error of Optimized DMD NLLS Problem is %g\n',relerr1_r);

%% Step 6: Plot DMD eigenvalues
% figure
% hold on
% scatter(real(omega),imag(omega),'bd','LineWidth',2)
% scatter(real(e1),imag(e1),'rx','LineWidth',2)
% box on
% grid on
% legend('Standard DMD eigenvalues','Optimized DMD eigenvalues')
% title(['Modes for Rank ' num2str(r), ' Approximation'])

%% Step 7: Plot Modes and Error
Nx = 161;
Nt = 200;

x = linspace(0,1,Nx);
t = linspace(0,10,Nt);

% Standard DMD

ROM_error_DMD = simDMD...
    (Xdmd,X1,x,t,Nt,M,'Standard DMD Model',r,'z_{DMD}', 'z-z_{DMD}');
% fprintf('At P_1, p_1 = %g, p_2 = %g, and p_3 = %g,\n',epsilon,q1,q2);
% fprintf('ROM error of DMD is %g\n',ROM_error_DMD);

% Optimized DMD

ROM_error_optDMD = simDMD...
    (Xd1(:, 1:end-1),X1,x,t,Nt,M,'Optimized DMD Model',r,'z_{OD}', 'z-z_{OD}');
% fprintf('At P_1, p_1 = %g, p_2 = %g, and p_3 = %g,\n',epsilon,q1,q2);
% fprintf('ROM error of optDMD is %g\n',ROM_error_optDMD);

end

function [ROM_error_DMD] = simDMD(Xdmd,X1,x,t,Nt,M,name,r,model1,model2)

realXdmd = real(Xdmd);
% figure
% mesh(x,t,realXdmd'); 
% view([1.0 .1 .1])
% title([name, ', order = ' , num2str(r) ] )
% xlabel('x'); ylabel('t'); zlabel(model1)
% xlim([0 1]) 
% ylim([0 10])
% zlim([-0.1 0.5])

error = X1 - realXdmd;
% figure
% mesh(x,t,error')
% view([1.0 .1 .1])
% xlabel('x'); ylabel('t'); zlabel(model2)
% title(['Error in ', name, ', order = ' , num2str(r)] )
% xlim([0 1]) 
% ylim([0 10])
% zlim([-0.25 0.25])

delta_t = (t(2)-t(1))/Nt;

ROM_error_DMD = 0.5*error(:,1)'*M*error(:,1);
  for i=2:Nt-1
    ROM_error_DMD = ROM_error_DMD + error(:,i)'*M*error(:,i);
  end
  ROM_error_DMD = ROM_error_DMD + 0.5*error(:,Nt)'*M*error(:,Nt);
  
  ROM_error_DMD = sqrt(delta_t*ROM_error_DMD); % (\int_0^T \| error(t,.) \|^2 dt)^{1/2}

end
