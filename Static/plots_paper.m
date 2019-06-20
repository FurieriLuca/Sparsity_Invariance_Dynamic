
%% The figure of performance comparsion in our ECC paper
%clc;clear;close all
% First plot all
%load PaperEx1.mat
lead=0:1:16;
figure;
h1 = plot(lead(1,:),J(:,1),'-rs','LineWidth',1.2,'MarkerSize',8); hold on;
h2 = plot(lead(1,:),J(:,3),'-g+','LineWidth',1.2,'MarkerSize',8); hold on;
h3 = plot(lead(1,:),J(:,5),'-bo','LineWidth',1.2,'MarkerSize',8); hold on;
h4 = plot(lead(1,:),J(:,7),'-kv','LineWidth',1.2,'MarkerSize',8);

% Then format
xlim([0 16]); ylim([16 23])
xlabel('Number of nodes that have full information: $L$','Interpreter','latex');
ylabel('$\mathcal{H}_2$ norm','Interpreter','latex');
set(gca,'TickLabelInterpreter','latex','fontsize',13);

% Set legend
h = legend([h1,h4],'block-diagonal','sparsity invirance with $T = S$', 'sparsity invirance with $T = T_{new}$','centralized');
set(h,'box','off','orientation','vertical','FontSize',11,...
    'Position',[0.1 0.4 0.7 0.03],'Interpreter','latex')

set(gcf,'Position',[250 150 450 350]);
print(gcf,['figs/Performance'],'-painters','-depsc2','-r600')



