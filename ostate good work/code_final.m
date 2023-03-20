%% One dimensional problem
% March 15, 2023
%--- CKF, CKF-RD, ECKF-RD

clc, clear
close all
%% Given parameters
Q = 10;
R = 1;
p = [0.2,0.2 0.3 0.5]; %0.2 delayed probability
del=[3  , 6  ,2 ,    2]; %2..6
M=800;
n = 1;
ny = 1;

%% CQKF point and weight generation
n1 = 2;
[cp, W] = cqkf_p(n, n1);
[~, col] = size(W);
w1 = sqrt(diag(W));
%% GHF points and weights
% alpha = 3;
% [cp, W] = gh_pt(alpha);
% [row, col] = size(cp);
% w1 = sqrt(diag(W));
%% UKF points and weights
% kappa = 3-n;
% cp = sqrt(n+kappa)*[zeros(n,1) eye(n) -eye(n)]; %sqrt(n+kappa)*[0 1 -1];
% W = [kappa/(n+kappa) 1/(2*(n+kappa))*ones(1, 2*n)];
% [row, col] = size(W);
% w1 = sqrt(diag(W));
%% User defined function
phi = @(x1, x2) 0.5*x1 + 25*x1./(1+x1.^2) +8*cos(1.2*x2);
gam = @(x3) (x3.^2)./20;


mmax = 200;
kmax = 100;
ep = zeros(mmax, kmax); % error

%% error
for i = 1:length(del)
g{i}.err = zeros(4, kmax); %CKF
end

for m = 1:mmax
    fprintf('MC run:%d\n',m)
    %% Initial estimated value
    xp = zeros(n, kmax);
    pp = 1;
    xp(1) = 0;
   for i=1:length(del)
    %% Truth generation
    x = zeros(n, kmax); z = zeros(ny, kmax);
    y = zeros(ny, kmax);
    x(1) = normrnd(xp(1), sqrt(pp));
    z(1) = gam(x(1)) + normrnd(0, sqrt(R));
    y(1) = z(1);
    xi1 = rand(1, kmax) <= p(i);
    for k = 1:kmax-1
        x(k+1) = phi(x(k), k) + normrnd(0, sqrt(Q));
        z(k+1) = gam(x(k+1)) + normrnd(0, sqrt(R));
        y(k+1) = (1-xi1(k+1))*z(k+1) + xi1(k+1)*z(k);
    end
     [er gamt] = Filter(n, ny, Q, R, cp, W, col,M, phi, gam, y, xp, pp, kmax,xi1,z, x,del(i),p(i));
     g{i}.err= g{i}.err+er;

   end


end
%%%======== RMSE Calculation ==========
for i=1:length(del)
    f{i}.rmse = sqrt(g{i}.err/mmax); 
end

f{1}.plotStyle1 = {'r-.','g-.','b-.','m-.'};
f{1}.plotStyle2 = {'k:','r-.','b-.','m-.'};
%  g1=f2;
 save(strcat('./rmse.mat'),'f','del','p')
figure(1)
box on 
hold on
for i=1:4
    plot(f{1}.rmse(i,:), f{1}.plotStyle2{i}, 'Linewidth', 1.5),
end
legend('CKF-OD','ECKF-OD','PF-OD','EPF-OD ');

%% --- CKF is performing surprisingly better than GHF and CQKF.
count = nnz(gamt)/kmax