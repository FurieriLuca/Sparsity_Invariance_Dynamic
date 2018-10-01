leaders_vec=0:1:L;


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

plot(leaders_vec(1,:),J(:,1),'-k','LineWidth',2);
plot(leaders_vec(1,:),J(:,2),'--k','LineWidth',2);
plot(leaders_vec(1,:),J(:,3),'-b','LineWidth',2);
plot(leaders_vec(1,:),J(:,4),'--b','LineWidth',2);
plot(leaders_vec(1,:),J(:,5),'--r','LineWidth',2);



legend('block-diagonal','grad. desc. from block-diagonal','sparsity invariance', 'grad. desc. from sparsity invariance', 'centralized')

%plot(leaders_vec(1,:),J(:,1),'bx','LineWidth',2);


