% Math Camp HW
% Author: Paul Weifeng Dai

clear;
clc;
close;
%% Parameters

N = 5000;
T = 1000;
y = 0;
p = .3;
q = .2;

rng('default');

%% Q5.2(2)
emat0 = rand([N,T]); 
% +1 with prob p, 0 with prob 1-p-q, -1 with prob q
emat = (emat0<p)*1 + (emat0>=p).*(emat0<1-p-q+p)*0 + (emat0>=1-p-q+p)*(-1);

ymat = zeros([N,T]);
for col = 2:T
    ymat(:,col) = ymat(:,col-1)+emat(:,col);
end


mu = T*(p-q);
sigma2 = T*((p+q)-(p-q)^2);
ymin = min(ymat(:,T));
ymax = max(ymat(:,T));
ygrid = linspace(ymin,ymax,200);
ypdf = normpdf(ygrid,mu,sigma2^.5);

%% Markov Chain

%% Q5.2(1)
state = 200;
trans = zeros(state,state);
for i = 1:state
    for j = 1:state
        trans(i,j) = (j==i+1)* p + (j==i)*(1-p-q) +(j==i-1) *q;
    end
end

trans(1,1) = 1-p;
trans(state,state) = 1-q;

%% Q5.2(2)
pi0 = zeros(state,1);
pi0(1,1)= 1;

piT = trans'^T*pi0;

%% Q5.2(3)


figure(1);
subplot(2,2,1);
for row = 1:10
plot(1:T,ymat(row,:),'LineWidth',2.0);hold on;
end

subplot(2,2,2);
[f,xi] =  ksdensity(ymat(:,T));
plot(xi,f,'r','Linewidth',3.0);hold on;
plot(ygrid,ypdf,'b-','Linewidth',3.0);
legend('kernel density', 'theoretical');

subplot(2,2,3);
plot(1:length(piT),piT,'g-','Linewidth',3.0);xlabel('state')



