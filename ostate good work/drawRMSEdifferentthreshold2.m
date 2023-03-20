function drawRMSEdifferentthreshold2
close all;
clear all;
clc; 
set(0,'defaulttextinterpreter','latex');
set(0, 'defaultAxesTickLabelInterpreter','latex');
set(0, 'defaultLegendInterpreter','latex');
load('rmse.mat') %threshold = 2  82%

% load('differentETlatency.mat')
f{1}.plotStyle1 = {'r-.','g-.','b-.','m-.'};
f{1}.plotStyle2 = {'k-','r-.','b--.','m-.'}; figure(1)
subplot(2,1,1)
box on 
 plot(f{1}.rmse(3,:), f{1}.plotStyle2{1}, 'Linewidth', 0.6),hold on
for i=1:2
    plot(f{i}.rmse(4,:), f{1}.plotStyle2{i+1}, 'Linewidth', 0.6),hold on
end
 grid, xlabel('Time(s)');
 ylabel('RMSE');
legend('PF-OD','EPF-OD ($\xi=3$) ','EPF-OD ($\xi=6$) ');
% title('EPF-OD  wiyh deffirent threshold')
ylim([0 15])
subplot(2,1,2)
box on 
 plot(f{1}.rmse(3,:), f{1}.plotStyle2{1}, 'Linewidth', 0.6),hold on
for i=3:4
    plot(f{i}.rmse(4,:), f{1}.plotStyle2{i-1}, 'Linewidth', 0.6),hold on
end
 grid, xlabel('Time(s)');
 ylabel('RMSE');
legend('PF-OD','EPF-OD ($\alpha=0.3$) ','EPF-OD ($\alpha=0.5$) ');
% title('EPF-OD with deffernt latency ')
ylim([0 15])
% ylim([5 20]);
% set(gca,'FontSize',10);
%  set(gca,'FontName','Times New Roman','FontWeight','normal','FontSize',10);
grid on;

f = figure(1);
f.Position = [500 100 560 400];
% 
filename = 'RMSE1.eps';
print('-depsc2', filename, '-r100');
ans = filename % return the filename to org-mode