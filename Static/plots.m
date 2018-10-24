lead=0:1:L;

figure(1)
hold on;
grid on;
xticks(0:1:L)

%{
for(i=1:size(x_columns,1)/2)
    plot(time_vec_1(1,:),x_columns(i,:)-20,'--r','LineWidth',2);
    plot(time_vec_1(1,:),x_columns(i,:)-20,'ro','LineWidth',2,'MarkerSize', 20);
end
%}

plot(lead(1,:),J(:,1),'-kv','LineWidth',3,'MarkerSize',15);
%plot(leaders_vec(1,:),J(:,2),'--ks','LineWidth',3,'MarkerSize',15);
plot(lead(1,:),J(:,3),'-bo','LineWidth',4,'MarkerSize',15);
%plot(leaders_vec(1,:),J(:,4),'--bd','LineWidth',3,'MarkerSize',15);

plot(lead(1,:),J(:,5),'-gx','LineWidth',3,'MarkerSize',15);
%plot(leaders_vec(1,:),J(:,6),'--g*','LineWidth',3,'MarkerSize',15);

plot(lead(1,:),J(:,7),'--r','LineWidth',3,'MarkerSize',15);
%plot(leaders_vec(1,:),J(:,6),'--g','LineWidth',2);

set(gca,'fontsize',30,'fontweight','bold')
xlabel('L','fontweight','bold','fontsize',30);
ylabel('H2 norm','fontweight','bold','fontsize',30);


%legend('block-diagonal','GD from block-diagonal','spars. inv. with T = S', 'GD from spars. inv. with T=S', 'spars. inv. with T = max. cliques','GD from spars. inv. with T = max. cliques','centralized')
legend('block-diagonal','spars. inv. with T = S', 'spars. inv. with T = T_{new}','centralized')

%plot(leaders_vec(1,:),J(:,1),'bx','LineWidth',2);

figure(2)
hold on;
grid on;
xticks(0:1:L)
yticks(0:2:2*L)



plot(lead(1,:),degrees(:,1),'-bo','LineWidth',3,'MarkerSize',15);
plot(lead(1,:),degrees(:,2),'-gx','LineWidth',3,'MarkerSize',15);
set(gca,'fontsize',30,'fontweight','bold')
xlabel('L','fontweight','bold','fontsize',30);
ylabel('separation','fontweight','bold','fontsize',30)
legend('spars. inv. with T = S','spars. inv. with T = max. cliques');


