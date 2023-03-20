function drawRMSEdifferentthreshold
close all;
clear all;
clc; 
set(0,'defaulttextinterpreter','latex');
set(0, 'defaultAxesTickLabelInterpreter','latex');
set(0, 'defaultLegendInterpreter','latex');
load('rmse.mat') %threshold = 2  82%
% load('differentETlatency.mat')
f{1}.plotStyle1 = {'r-.','g-.','b-.','m-.'};
f{1}.plotStyle2 = {'k-','bx-','k-.','b-.'};
% load('differentETlatency.mat')
 plotStyle = {'b-','r-','k-'};
 lin =[0.5 0.5 1 1]
 figure(1) 
subplot(2,2,1)
box on 
for i=1:4
    plot(f{1}.rmse(i,:), f{1}.plotStyle2{i}, 'Linewidth', lin(i)),hold on
end
 xlabel('Time(s)');
 ylabel('RMSE');
 legend('CKF-OD','ECKF-OD','PF-OD ','EPF-OD ');
% title('with fixed latency $\alpha=0.2$ and $\xi=3$')
ylim([0 22])

subplot(2,2,2)
box on 
for i=1:4
    plot(f{2}.rmse(i,:), f{1}.plotStyle2{i}, 'Linewidth', 0.8),hold on
end
 grid, xlabel('Time(s)');
 ylabel('RMSE');
legend('CKF-OD','ECKF-OD','PF-OD ','EPF-OD ');
title('with fixed latency $\alpha=0.2$ and $\xi=6$')
ylim([0 21])

subplot(2,2,3)
box on 
for i=1:4
    plot(f{3}.rmse(i,:), f{1}.plotStyle2{i}, 'Linewidth', 0.6),hold on
end
 grid, xlabel('Time(s)');
 ylabel('RMSE');ylim([0 15])

legend('CKF-OD','ECKF-OD','PF-OD ','EPF-OD ');
title('with fixed threshold  $\xi=2$ and $\alpha=0.3$ ' )
subplot(2,2,4)
box on 
for i=1:4
    plot(f{4}.rmse(i,:), f{1}.plotStyle2{i}, 'Linewidth', 0.6),hold on
end
 grid, xlabel('Time(s)');
 ylabel('RMSE');
legend('CKF-OD','ECKF-OD','PF-OD ','EPF-OD ');
title('with fixed threshold  $\xi=2$ and $\alpha=0.5$')
 ylim([5 15]);
% set(gca,'FontSize',10);
%  set(gca,'FontName','Times New Roman','FontWeight','normal','FontSize',10);
grid on;

f = figure(1);
f.Position = [500 100 560 400];
% 
filename = 'RMSESCP1.eps';
print('-depsc2', filename, '-r100');
ans = filename % return the filename to org-mode