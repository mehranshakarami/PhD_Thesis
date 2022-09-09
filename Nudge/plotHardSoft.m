%%%% hard and soft
clear all %#ok<*CLALL>
close all
global  x_store h del_p p0 x_star1 x_star2 delta_bar a b N n delta eta lam_hat dmax d maxRate price price_time;

load 'Data\PEVsData.mat';
%% price
p0=0.3*ones(n,1);
del_p=0.1;
delta_bar=0.15; %<min{delta}-del_p
if delta_bar>=min(delta)-del_p
    error('Delta bar is too big!');
end
rho=0.3-del_p;
%% finding desired behaviors
options = optimoptions('quadprog','Display','off');
ps0=0.1*[3.23    3.23    3.23    3.23    3.23    3.23  ...
    3.23    3.23    3.23    3.23    3.23 3.06 ...
    2.98    2.93    2.91    2.90    2.898    2.899  ...
    2.903    2.92    2.94    2.97    3.03    3.23]'; % the price for x_star1

if norm(ps0-p0)>delta_bar
    error('P_star out of the ball!');
end
x_star1=zeros(n,1);
for i=1:N
    x0=[];
    x_star1 =x_star1+quadprog(2*a(i)*eye(n),b(i)*ones(n,1)+ps0,[],[],ones(n,1)',d(i),zeros(n,1),maxRate(i)*ones(n,1),x0,options);
end
%% dynamics
%%%generation of price
price_time=0:0.001:1;
price=zeros(numel(price_time),n);
for i=1:numel(price_time)
delp_temp=1-2*rand(n,1);
delp_normlized=delp_temp/norm(delp_temp);
price(i,:)=transpose(p0+del_p*rand(1)*delp_normlized);
end
gamma0=0.7*rand(N,1);

x_store=zeros(n,N);
[t_hard1,state_hard1]=ode45(@hard,[0 0.003],[p0;gamma0],odeset('Maxstep',1e-5));%,odeset('Events',@eventsFcn));%,odeset('Maxstep',1e-2));
[t_hard2,state_hard2]=ode45(@hard,[0.003 1],state_hard1(end,:)',odeset('Maxstep',1e-3));
t_hard=[t_hard1;t_hard2];
state_hard=[state_hard1;state_hard2];
p_hat_hard=state_hard(:,1:n);
gamma_hard=state_hard(:,n+1:n+N);

x_store=zeros(n,N);
[t_soft1,state_soft1]=ode45(@soft,[0 0.003],[p0+0.06*ones(n,1);gamma0],odeset('Maxstep',1e-5));%,odeset('Events',@eventsFcn));%,odeset('Maxstep',1e-2));
[t_soft2,state_soft2]=ode45(@soft,[0.003 1],state_soft1(end,:)',odeset('Maxstep',1e-3));
t_soft=[t_soft1;t_soft2];
state_soft=[state_soft1;state_soft2];
p_hat_soft=state_soft(:,1:n);
gamma_soft=state_soft(:,n+1:n+N);
%% plots
%%% calculations
x_store=zeros(n,N); % for increasing the speed to use as intial condition
x_er_norm_hard=zeros(1,numel(t_hard)); % norm of agg. error
x_er_norm_soft=zeros(1,numel(t_soft));
x_sum_hard=zeros(n,numel(t_hard));
x_sum_soft=zeros(n,numel(t_soft));
p_hat_dist_p0_hard=zeros(numel(t_hard),1); % distance of p_hat to p0
p_hat_dist_p0_soft=zeros(numel(t_soft),1);
for k=1:numel(t_hard)
    for i=1:N
        lam= gamma_hard(k,i)*p_hat_hard(k,:)'+(1-gamma_hard(k,i))*lam_hat(:,i);
        x0=x_store(:,i);
        x_sum_hard(:,k) =x_sum_hard(:,k)+quadprog(2*a(i)*eye(n),b(i)*ones(n,1)+lam,[],[],ones(n,1)',d(i),zeros(n,1),maxRate(i)*ones(n,1),x0,options);
    end
    x_er_norm_hard(k)=norm(x_sum_hard(:,k)-x_star1);
    p_hat_dist_p0_hard(k)=sqrt((p_hat_hard(k,:)-0.3)*(p_hat_hard(k,:)-0.3)');
end
for k=1:numel(t_soft)
    for i=1:N
        lam= gamma_soft(k,i)*p_hat_soft(k,:)'+(1-gamma_soft(k,i))*lam_hat(:,i);
        x0=x_store(:,i);
        x_sum_soft(:,k) =x_sum_soft(:,k)+quadprog(2*a(i)*eye(n),b(i)*ones(n,1)+lam,[],[],ones(n,1)',d(i),zeros(n,1),maxRate(i)*ones(n,1),x0,options);
    end
    x_er_norm_soft(k)=norm(x_sum_soft(:,k)-x_star1);
    p_hat_dist_p0_soft(k)=sqrt((p_hat_soft(k,:)-0.3)*(p_hat_soft(k,:)-0.3)');
end
%%
%%% x_star
figure
box on
hold on
grid on
plot(x_star1,'linewidth',2);

xlabel('Charging horizon','Interpreter','latex','fontsize',12)
ylabel('$x^*$(kW)','Interpreter','latex','fontsize',12)


%%% distance of phat to p0 and gamma
%set(groot,'defaultAxesTickLabelInterpreter','latex');  
figure
subplot 211
box on
hold on
grid on
plot(t_hard,p_hat_dist_p0_hard,'linewidth',2)
plot(t_soft,p_hat_dist_p0_soft,'r','linewidth',2)
plot(t_hard,ones(size(t_hard)).*delta_bar,'k--','linewidth',1.5)
plot(t_hard,ones(size(t_hard)).*rho,'--','color',[0, 0.5, 0],'linewidth',1.5)
xlabel('Time','Interpreter','latex','fontsize',12)
ylabel('$\|\hat{p}-p_0\|$(\$/kWh)','Interpreter','latex','fontsize',12)
% yticks([0.15 0.2])
% yticklabels({'$\bar{\delta}=0.15$','$\rho=0.2$'});
h=legend('Hard nudge','Soft nudge','$\bar{\delta}$','$\rho$');
set(h,'Interpreter','latex','fontsize',12)
h.NumColumns = 4;
% ylim([0 0.28])

subplot 212
box on
hold on
grid on
plot(t_hard,gamma_hard*ones(N,1)/N,'linewidth',2)
plot(t_soft,gamma_soft*ones(N,1)/N,'r','linewidth',2)
xlabel('Time','Interpreter','latex','fontsize',12)
ylabel('$\frac{1}{N}\sum\limits_{i\in\mathcal{I}}\gamma_i$','Interpreter','latex','fontsize',12)
h=legend('Hard nudge','Soft nudge');
set(h,'Interpreter','latex','fontsize',12)
ylim([0 1.03])

%%% agg behvaior and tracking error
figure
subplot(2,2,1);
box on
hold on
grid on
plot(t_hard,x_sum_hard,'linewidth',2);
title('Hard nudge','Interpreter','latex','fontsize',12)
ylabel('$\sum\limits_{i\in\mathcal{I}}x_i$(kW)','Interpreter','latex','fontsize',12)
xlabel('Time','Interpreter','latex','fontsize',12)

subplot(2,2,2);
box on
hold on
grid on
plot(t_soft,x_sum_soft,'linewidth',2);
title('Soft nudge','Interpreter','latex','fontsize',12)
ylabel('$\sum\limits_{i\in\mathcal{I}}x_i$(kW)','Interpreter','latex','fontsize',12)
xlabel('Time','Interpreter','latex','fontsize',12)

subplot(2,2,[3,4]);
box on
hold on
grid on
plot(t_hard,x_er_norm_hard,'linewidth',2);
plot(t_soft,x_er_norm_soft,'r','linewidth',2);
ylabel('$\|\sum\limits_{i\in\mathcal{I}}x_i-x^*\|$(kW)','Interpreter','latex','fontsize',12)
xlabel('Time','Interpreter','latex','fontsize',12)
h=legend('Hard nudge','Soft nudge');
set(h,'Interpreter','latex','fontsize',12)




