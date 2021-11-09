% Math Camp HW
% Author: Paul Weifeng Dai

clear;
clc;
close;
%% Parameters

global beta nkgrid iter_max err_tol kgrid

% exo parameter
beta = .98;
alpha = 1/3;
delta = 1;
z = 1;

% steady state
nkgrid = 200;
k_ss = (alpha*z/(1/beta-(1-delta)))^(1/(1-alpha));
kmin = 0.05*k_ss;
kmax = 1.2*k_ss;

kgrid = linspace(kmin,kmax,nkgrid);

% tolerance and iteration setting
iter_max = 10^3;
err_tol = 10^-5;
infty = 10^10;
%% VFI
% initial guess
Vg1 = zeros(nkgrid,1);
Vg2 = log(z*kgrid.^alpha-delta*kgrid)/(1-beta);
Vg2 = Vg2';


[kk,kklead] = ndgrid(kgrid,kgrid);
c = z*kk.^alpha+(1-delta)*kk - kklead;
u = (c<=0).*(-infty) + (c>0).*log(c);

disp('VFI with Guess 1');
[opt_k1] = VFI(u,Vg1);
disp('VFI with Guess 2');
[opt_k2] = VFI(u,Vg2);


%% Plot PC
figure(1);
subplot(1,2,1);
plot(kgrid,kgrid,'b','LineWidth',3.0);hold on;
plot(kgrid,opt_k1,'r','LineWidth',3.0);hold on;
legend('45 degree line','Policy function (Using Guess 1)');

subplot(1,2,2);
plot(kgrid,kgrid,'b','LineWidth',3.0);hold on;
plot(kgrid,opt_k2,'r','LineWidth',3.0);hold on;
legend('45 degree line','Policy function (Using Guess 2)');

%% Function: VFI
function [opt_k] = VFI(u,V)

global beta nkgrid iter_max err_tol kgrid

iter = 0;
err = 10^9;

tic
while iter<iter_max && err>err_tol
Vtemp = u + beta*permute(repmat(V,[1,nkgrid]),[2,1]);
[Vnew,opt_k_idx] = max(Vtemp,[],2);
err = max(abs(Vnew-V),[],'all');
V = Vnew;
end

if iter<iter_max && err<err_tol
    disp('VFI is done')
else
    disp('Not converged')
end

opt_k = kgrid(opt_k_idx);
toc
end
