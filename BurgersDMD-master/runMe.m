%  A script to compute a POD model for Burgers equation.  This script drives
%  the following three steps:
%
%   (1)  compute the simulation data:  using burgers_1d_periodic.
%
%   (2)  use the simulation to compute a POD basis: using generate_pod.m
%
%   (3)  run test_model.m that builds the low-dimensional dynamical
%        system, simulates it, reconstructs an approximation to the solution
%        to Burgers equation, then computes the error in the approximation
%        (difference between FEM solution and POD solution).
%
%
%  Copyright 2019, Jeff Borggaard, Virginia Tech
%
%  Permission is hereby granted, free of charge, to any person obtaining 
%  a copy of this software and associated documentation files (the "Software"),
%  to deal in the Software without restriction, including without limitation 
%  the rights to use, copy, modify, merge, publish, distribute, sublicense, 
%  and/or sell copies of the Software, and to permit persons to whom the 
%  Software is furnished to do so, subject to the following conditions:
%  
%  The above copyright notice and this permission notice shall be included 
%  in all copies or substantial portions of the Software.
%
%  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS 
%  OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, 
%  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL 
%  THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER 
%  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING 
%  FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER 
%  DEALINGS IN THE SOFTWARE.
%
%%

clear, close all

% choose rank loop
N = 25;
space = 5;
n = N/space;

% choose reynolds loop
ct = 1;

% array of rank approximations for a reynolds number
dim = zeros(1,n);
for i = 1:n
    dim(1,i) = space*i ;
end

% collection of errors for plotting
PODerrF_full = zeros(ct,n);
DMDerrF_full = zeros(ct,n);
ODerrF_full = zeros(ct,n);
PODerr2_full = zeros(ct,n);
DMDerr2_full = zeros(ct,n);
ODerr2_full = zeros(ct,n);

iter1 = 0;
timeOD = zeros(1,ct*n);
timeDMD = zeros(1,ct*n);
timePOD = zeros(1,ct*n);
timespace = linspace(1,ct*n,ct*n);

% Reynolds numbers re
for re = 1000:1000:(1000*ct)
    
    %count iteration
    iter1 = iter1 + 1;
    
    % Create the simulation data at parameter "p1"
    R  = re;
    q1 = 0.5;
    q2 = 0.0;
    Nx = 80; 
    dx = 1/Nx;
    Nt = 201;
    dt = 1/Nt;
    epsilon = 1/R;
    p1 = [ epsilon; q1; q2 ]; 
    
    % initialize array of errors
    PODerrF = zeros(1,n);
    DMDerrF = zeros(1,n);
    ODerrF = zeros(1,n);
    PODerr2 = zeros(1,n);
    DMDerr2 = zeros(1,n);
    ODerr2 = zeros(1,n);
    
    iter2 = 0;
    % Rank approximation rk
    for rk = 1:n
        % stopwatch
%         tic
        
        % count iteration
        iter2 = iter2 + 1;

        % the number of POD basis functions to use
        r_dim  = dim(1,rk);   

        % the 80 and 201 below are spatial and temporal FEM discretization parameters
        [z,x,t,e_conn] = burgers_1d_periodic(epsilon,q1,q2,Nx,Nt);

          %  Comment out the eye candy below for parametric experiments...
          figure
          mesh(x,t,z')
        %   view([.7 -.8 .6])
          view([1.0 .1 .1])
          pause(0.001)
          xlabel('x'); ylabel('t'); zlabel('z')
          title('Finite Element Solution')
          xlim([0 1]) 
          ylim([0 10])
          zlim([-0.1 0.5])

        % Use the simulation data to compute a POD basis of dimension r_dim
        [M] = compute_mass_matrix(x,e_conn);

        % [POD1,lam1] = generate_pod(z,r_dim,M);
        
        tic
        [POD] = run_POD(z,r_dim,M);

        % Assemble the nonlinear ROM, perform simulation, and compute the error

        [ROM_error_POD,relerrPOD_r] = test_model(POD,r_dim,epsilon,q1,q2,M,z);
        tPOD = toc;
        
        %   fprintf('Inputs: epsilon = %g, q1 = %g, q2 = %g, r_dim = %g, Nx = %g, Nt = %g \n',...
        %       epsilon,q1,q2,r_dim,Nx,Nt);

        % Use the POD modes to compute a DMD bases of dimension r_dim

        [relerrDMD_r,relerr1_r,ROM_error_DMD,ROM_error_optDMD,tOD,tDMD] = run_DMD(z,r_dim,dt,M,t);

        %   Collect error norms
        PODerrF(1,rk) = relerrPOD_r;
        DMDerrF(1,rk) = relerrDMD_r;
        ODerrF(1,rk) = relerr1_r;
        PODerr2(1,rk) = ROM_error_POD;
        DMDerr2(1,rk) = ROM_error_DMD;
        ODerr2(1,rk) = ROM_error_optDMD;
        
        % For 5 Re, 10 rk
        timeOD(1, (iter2 - 1)*5 + iter1) = tOD;
        timeDMD(1, (iter2 - 1)*5 + iter1) = tDMD;
        timePOD(1, (iter2 - 1)*5 + iter1) = tPOD;
        
%         % For 10 Re, 5 rk
%         timeOD(1, (iter2 - 1)*10 + iter1) = tOD;
%         timeDMD(1, (iter2 - 1)*10 + iter1) = tDMD;
%         timePOD(1, (iter2 - 1)*10 + iter1) = tPOD;
    end

    PODerrF_full(iter1,:) = PODerrF;
    DMDerrF_full(iter1,:) = DMDerrF;
    ODerrF_full(iter1,:) = ODerrF;
    PODerr2_full(iter1,:) = PODerr2;
    DMDerr2_full(iter1,:) = DMDerr2;
    ODerr2_full(iter1,:) = ODerr2;
    
    % Error Plots for single Re
    figure('Name','Relative Error of ROM')
    semilogy(dim,PODerrF,'r-o',dim,DMDerrF,'b-+',dim,ODerrF,'m-x')
    title(['Relative ROM Error $ = \frac{\|Z-Z_r \|_F}{\|Z \|_F} $, Re = ', num2str(re)], 'Interpreter','latex')
    legend('POD', 'DMD', 'OD')
    xlim([dim(1,1) dim(1,n)])
    xlabel('Dimension (Order of Approximation)')
    ylabel('Size of Error')
    
    figure('Name','L2 Error of ROM')
    semilogy(dim,PODerr2,'r-o',dim,DMDerr2,'b-+',dim,ODerr2,'m-x')
    title(['$L_2$ ROM Error $\approx (\int_0^T \| r(t) \|^2 dt)^{\frac{1}{2}} $, Re = ', num2str(re)],'Interpreter','latex')
    legend('POD', 'DMD', 'OD')
    xlim([dim(1,1) dim(1,n)])
    xlabel('Dimension (Order of Approximation)')
    ylabel('Size of Error')
    
end

% % Error Plots for 5 Re

% figure('Name','Relative Error of ROM')
% semilogy(dim,PODerrF_full(1,:),'r-o',dim,DMDerrF_full(1,:),'b-+',dim,ODerrF_full(1,:),'m-x',...
%     dim,PODerrF_full(2,:),'r-o',dim,PODerrF_full(3,:),'r-o',dim,PODerrF_full(4,:),'r-o',dim,PODerrF_full(5,:),'r-o',...
%     dim,DMDerrF_full(2,:),'b-+',dim,DMDerrF_full(3,:),'b-+',dim,DMDerrF_full(4,:),'b-+',dim,DMDerrF_full(5,:),'b-+',...
%     dim,ODerrF_full(2,:),'m-x',dim,ODerrF_full(3,:),'m-x',dim,ODerrF_full(4,:),'m-x',dim,ODerrF_full(5,:),'m-x')
% title('Relative ROM Error $ = \frac{\|Z-Z_r \|_F}{\|Z \|_F} $', 'Interpreter','latex')
% legend('POD','DMD','OD')
% xlim([dim(1,1) dim(1,end)])
% xlabel('Dimension (Order of Approximation)')
% ylabel('Size of Error')
% 
% figure('Name','L2 Error of ROM')
% semilogy(dim,PODerr2_full(1,:),'r-o',dim,DMDerr2_full(1,:),'b-+',dim,ODerr2_full(1,:),'m-x',...
%     dim,PODerr2_full(2,:),'r-o',dim,PODerr2_full(3,:),'r-o',dim,PODerr2_full(4,:),'r-o',dim,PODerr2_full(5,:),'r-o',...
%     dim,DMDerr2_full(2,:),'b-+',dim,DMDerr2_full(3,:),'b-+',dim,DMDerr2_full(4,:),'b-+',dim,DMDerr2_full(5,:),'b-+',...
%     dim,ODerr2_full(2,:),'m-x',dim,ODerr2_full(3,:),'m-x',dim,ODerr2_full(4,:),'m-x',dim,ODerr2_full(5,:),'m-x')
% title('$L_2$ ROM Error $\approx (\int_0^T \| r(t) \|^2 dt)^{\frac{1}{2}} $','Interpreter','latex')
% legend('POD','DMD','OD')
% xlim([dim(1,1) dim(1,end)])
% xlabel('Dimension (Order of Approximation)')
% ylabel('Size of Error')

% % Time plot avec DMD

% figure('Name','Computation Time')
% plot(timespace,timePOD(1,:),'r-o',timespace,timeDMD(1,:),'b-+',timespace,timeOD(1,:),'m-x')
% title('Computation Time')
% legend('POD', 'DMD', 'OD', 'Location','northwest')
% xlim([timespace(1,1) timespace(1,end)])
% xlabel('Loop, sorted by (i) Rank of Approximation (ii) Re')
% ylabel('Seconds')

% % Error Plots sans DMD for 5 Re

% figure('Name','Relative Error of ROM')
% semilogy(dim,PODerrF_full(1,:),'r-o',dim,ODerrF_full(1,:),'b-x',...
%     dim,PODerrF_full(2,:),'r-o',dim,PODerrF_full(3,:),'r-o',dim,PODerrF_full(4,:),'r-o',dim,PODerrF_full(5,:),'r-o',...
%     dim,ODerrF_full(2,:),'b-x',dim,ODerrF_full(3,:),'b-x',dim,ODerrF_full(4,:),'b-x',dim,ODerrF_full(5,:),'b-x')
% title('Relative ROM Error $ = \frac{\|Z-Z_r \|_F}{\|Z \|_F} $', 'Interpreter','latex')
% legend('POD','OD')
% xlim([dim(1,1) dim(1,end)])
% xlabel('Dimension (Order of Approximation)')
% ylabel('Size of Error')

% figure('Name','L2 Error of ROM')
% semilogy(dim,PODerr2_full(1,:),'r-o',dim,ODerr2_full(1,:),'b-x',...
%     dim,PODerr2_full(2,:),'r-o',dim,PODerr2_full(3,:),'r-o',dim,PODerr2_full(4,:),'r-o',dim,PODerr2_full(5,:),'r-o',...
%     dim,ODerr2_full(2,:),'b-x',dim,ODerr2_full(3,:),'b-x',dim,ODerr2_full(4,:),'b-x',dim,ODerr2_full(5,:),'b-x')
% title('$L_2$ ROM Error $\approx (\int_0^T \| r(t) \|^2 dt)^{\frac{1}{2}} $','Interpreter','latex')
% legend('POD','OD')
% xlim([dim(1,1) dim(1,end)])
% xlabel('Dimension (Order of Approximation)')
% ylabel('Size of Error')

% % Time plot sans DMD
% figure('Name','Computation Time')
% plot(timespace,timePOD(1,:),'r-o',timespace,timeOD(1,:),'b-x')
% title('Computation Time')
% legend('POD', 'OD', 'Location','northwest')
% xlim([timespace(1,1) timespace(1,end)])
% xlabel('Loop, sorted by (i) Rank of Approximation (ii) Re')
% ylabel('Seconds')

% % Error Plots for 10 Re
% figure('Name','L2 Error of ROM')
% semilogy(dim,PODerr2_full(1,:),'r-o',dim,DMDerr2_full(1,:),'b-+',dim,ODerr2_full(1,:),'m-x',...
%     dim,PODerr2_full(2,:),'r-o',dim,PODerr2_full(3,:),'r-o',dim,PODerr2_full(4,:),'r-o',dim,PODerr2_full(5,:),'r-o',...
%     dim,DMDerr2_full(2,:),'b-+',dim,DMDerr2_full(3,:),'b-+',dim,DMDerr2_full(4,:),'b-+',dim,DMDerr2_full(5,:),'b-+',...
%     dim,ODerr2_full(2,:),'m-x',dim,ODerr2_full(3,:),'m-x',dim,ODerr2_full(4,:),'m-x',dim,ODerr2_full(5,:),'m-x',...
%     dim,PODerr2_full(6,:),'r-o',dim,DMDerr2_full(6,:),'b-+',dim,ODerr2_full(6,:),'m-x',...
%     dim,PODerr2_full(7,:),'r-o',dim,PODerr2_full(8,:),'r-o',dim,PODerr2_full(9,:),'r-o',dim,PODerr2_full(10,:),'r-o',...
%     dim,DMDerr2_full(7,:),'b-+',dim,DMDerr2_full(8,:),'b-+',dim,DMDerr2_full(9,:),'b-+',dim,DMDerr2_full(10,:),'b-+',...
%     dim,ODerr2_full(7,:),'m-x',dim,ODerr2_full(8,:),'m-x',dim,ODerr2_full(9,:),'m-x',dim,ODerr2_full(10,:),'m-x')
% title('$L_2$ ROM Error $\approx (\int_0^T \| r(t) \|^2 dt)^{\frac{1}{2}} $','Interpreter','latex')
% legend('POD','DMD','OD')
% xlim([dim(1,1) dim(1,end)])
% xlabel('Dimension (Order of Approximation)')
% ylabel('Size of Error')
% 
% figure('Name','Relative Error of ROM')
% semilogy(dim,PODerrF_full(1,:),'r-o',dim,DMDerrF_full(1,:),'b-+',dim,ODerrF_full(1,:),'m-x',...
%     dim,PODerrF_full(2,:),'r-o',dim,PODerrF_full(3,:),'r-o',dim,PODerrF_full(4,:),'r-o',dim,PODerrF_full(5,:),'r-o',...
%     dim,DMDerrF_full(2,:),'b-+',dim,DMDerrF_full(3,:),'b-+',dim,DMDerrF_full(4,:),'b-+',dim,DMDerrF_full(5,:),'b-+',...
%     dim,ODerrF_full(2,:),'m-x',dim,ODerrF_full(3,:),'m-x',dim,ODerrF_full(4,:),'m-x',dim,ODerrF_full(5,:),'m-x',...
%     dim,PODerrF_full(6,:),'r-o',dim,DMDerrF_full(6,:),'b-+',dim,ODerrF_full(6,:),'m-x',...
%     dim,PODerrF_full(7,:),'r-o',dim,PODerrF_full(8,:),'r-o',dim,PODerrF_full(9,:),'r-o',dim,PODerrF_full(10,:),'r-o',...
%     dim,DMDerrF_full(7,:),'b-+',dim,DMDerrF_full(8,:),'b-+',dim,DMDerrF_full(9,:),'b-+',dim,DMDerrF_full(10,:),'b-+',...
%     dim,ODerrF_full(7,:),'m-x',dim,ODerrF_full(8,:),'m-x',dim,ODerrF_full(9,:),'m-x',dim,ODerrF_full(10,:),'m-x')
% title('Relative ROM Error $ = \frac{\|Z-Z_r \|_F}{\|Z \|_F} $', 'Interpreter','latex')
% legend('POD','DMD','OD')
% xlim([dim(1,1) dim(1,end)])
% xlabel('Dimension (Order of Approximation)')
% ylabel('Size of Error')

% save('Rom_POD_DMD_OD_5re_10n.mat')