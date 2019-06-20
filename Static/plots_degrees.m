
%% The figure of performance comparsion in our ECC paper
%clc;clear;close all
% First plot all
%load PaperEx1.mat
figure(2);
h1 = plot(lead(1,:),degrees(:,3),'-rs','LineWidth',1.2,'MarkerSize',10); hold on;
h2 = plot(lead(1,:),degrees(:,1),'-gp','LineWidth',1.2,'MarkerSize',10); hold on;
h3 = plot(lead(1,:),degrees(:,4),'-kv','LineWidth',1.2,'MarkerSize',10); hold on;
%h3 = plot(lead(1,:),J(:,5),'-bo','LineWidth',1.2,'MarkerSize',8); hold on;
%h4 = plot(lead(1,:),J(:,7),'-kv','LineWidth',1.2,'MarkerSize',8);

% Then format
xlim([0 16]); ylim([1 16])
yticks([1:1:16])
xlabel('\#nodes with full information: $L$','Interpreter','latex');
ylabel('separation degree $r$','Interpreter','latex');
set(gca,'TickLabelInterpreter','latex','fontsize',13);

% Set legend
h = legend([h1,h2,h3],'block-diagonal','sparsity invariance with $T = S$', 'centralized');
set(h,'box','off','orientation','vertical','FontSize',13,...
    'Position',[0.1 0.4 0.7 0.03],'Interpreter','latex')

set(gcf,'Position',[250 150 450 350]);
print(gcf,['figs/Performance'],'-painters','-depsc2','-r600')



