clear all %#ok<*CLALL>
close all
global  x_store h del_p p0 x_star1 x_star2 delta_bar a b N n delta eta lam_hat dmax d maxRate;

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
    3.23    3.23    3.23    3.23    3.23 3.08 ...
    3.0    2.95    2.93    2.92    2.918    2.919  ...
    2.923    2.94    2.96    2.99    3.05    3.23]'; % the price for x_star1
ps1=0.1*[3.24    3.24    3.24    3.24    3.24    3.24  ...
    3.24    3.24    3.24    3.1    3.06 3.03 ...
    3    2.99    2.98    2.975    1.02*2.918    1.02*2.919  ...
    1.02*2.923    1.02*2.94    1.02*2.96    1.02*2.99    1.02*3.05    3.24]'; % the price for x_star2
if norm(ps1-p0)>delta_bar || norm(ps0-p0)>delta_bar
    error('P_star out of the ball!');
end
x_star1=zeros(n,1);
x_star2=zeros(n,1);
for i=1:N
    x0=[];
    x_star1 =x_star1+quadprog(2*a(i)*eye(n),b(i)*ones(n,1)+ps0,[],[],ones(n,1)',d(i),zeros(n,1),maxRate(i)*ones(n,1),x0,options);
    x_star2 =x_star2+quadprog(2*a(i)*eye(n),b(i)*ones(n,1)+ps1,[],[],ones(n,1)',d(i),zeros(n,1),maxRate(i)*ones(n,1),x0,options);
end
%% dynamics
gamma0=0.7*rand(N,1);
% K_initial=zeros(n);
K_initial=0;
x_store=zeros(n,N);
% p0+0.06*ones(n,1)
[t,state]=ode45(@adap,[0 2],[p0+0.06*ones(n,1);gamma0;K_initial(:)],odeset('Maxstep',1e-3));%,odeset('Events',@eventsFcn));%,odeset('Maxstep',1e-2));
p_hat=state(:,1:n);
gamma=state(:,n+1:n+N);
K=state(:,n+N+1:end);

figure
plot(t,K)
%% plots
%%% calculations
x_store=zeros(n,N); % for increasing the speed to use as intial condition
x_er_norm=zeros(1,numel(t)); % norm of agg. error
x_sum=zeros(n,numel(t));
xs=x_sum; % x_star at each time
p_hat_dist_p0=zeros(numel(t),1); % distance of p_hat to p0
for k=1:numel(t)
    for i=1:N
        lam= gamma(k,i)*p_hat(k,:)'+(1-gamma(k,i))*lam_hat(:,i);
        x0=x_store(:,i);
        x_sum(:,k) =x_sum(:,k)+quadprog(2*a(i)*eye(n),b(i)*ones(n,1)+lam,[],[],ones(n,1)',d(i),zeros(n,1),maxRate(i)*ones(n,1),x0,options);
    end
    xs(:,k)=x_star(t(k));
    x_er_norm(k)=norm(x_sum(:,k)-xs(:,k));
    p_hat_dist_p0(k)=sqrt((p_hat(k,:)-0.3)*(p_hat(k,:)-0.3)');
end
%%
%%% x_star1 and x_star2, x_star
figure
subplot 211
box on
hold on
grid on
plot(x_star1,'linewidth',2);
plot(x_star2,'r','linewidth',2);
h=legend('$m$','$s$');
set(h,'Interpreter','latex','fontsize',12)
xlabel('Charging horizon','Interpreter','latex','fontsize',12)
ylabel('Power demand(kW)','Interpreter','latex','fontsize',12)

subplot 212
box on
hold on
grid on
plot(t,xs,'linewidth',2);
ylabel('$x^*$(kW)','Interpreter','latex','fontsize',12)
xlabel('Time','Interpreter','latex','fontsize',12)
% ylabel('$x^*$','Interpreter','latex','fontsize',12)

%%% distance of phat to p0

figure
subplot(2,2,1:2)
box on
hold on
grid on
plot(t,p_hat_dist_p0,'linewidth',2)
plot(t,ones(size(t)).*delta_bar,'k--','linewidth',1.5)
plot(t,ones(size(t)).*rho,'--','color',[0, 0.5, 0],'linewidth',1.5)
xlabel('Time','Interpreter','latex','fontsize',12)
ylabel('$\|\hat{p}-p_0\|$(\$/kWh)','Interpreter','latex','fontsize',12)
h=legend('Adaptive nudge','$\bar{\delta}$','$\rho$');
set(h,'Interpreter','latex','fontsize',12)
% ylim([0 0.28])
h.NumColumns = 3;

subplot(2,2,3)
box on
hold on
grid on
plot(t,gamma*ones(N,1)/N,'linewidth',2)
xlabel('Time','Interpreter','latex','fontsize',12)
ylabel('$\frac{1}{N}\sum\limits_{i\in\mathcal{I}}\gamma_i$','Interpreter','latex','fontsize',12)
set(h,'Interpreter','latex','fontsize',12)
ylim([0 1.03])

subplot(2,2,4)
box on
hold on
grid on
plot(t,K,'linewidth',2)
xlabel('Time','Interpreter','latex','fontsize',12)
ylabel('$k$','Interpreter','latex','fontsize',12)

%%% agg behvaior and tracking error
figure
subplot 211
box on
hold on
grid on
plot(t,x_sum,'linewidth',2);
ylabel('$\sum\limits_{i\in\mathcal{I}}x_i$(kW)','Interpreter','latex','fontsize',12)
xlabel('Time','Interpreter','latex','fontsize',12)

subplot 212
box on
hold on
grid on
plot(t,x_er_norm,'linewidth',2);
ylabel('$\|\sum\limits_{i\in\mathcal{I}}x_i-x^*\|$(kW)','Interpreter','latex','fontsize',12)
xlabel('Time','Interpreter','latex','fontsize',12)





%% adaptive nudge
function ds=adap(t,s)
global  x_store x_star1 h x_star2 del_p p0 delta_bar a b N n delta eta lam_hat d maxRate;
ep=2e-5;
K0=10;
sigma=1e5;
tau=1;
% if ep>1/(1.5*norm(x_star1-x_star2)*(1+(0.5*N/0.004)))
%     error('ep is out of interval!')
% elseif sigma<2*(1.5*norm(x_star1-x_star2))*(1+(0.5*N/0.004))
%     error('sigma is out of interval!')
% elseif K0<sqrt(n)/(0.5*N/0.006)
%     error('K0 is out of interval!')
% end
p=s(1:n);
gamma=s(n+1:n+N);
% K=reshape(s(n+N+1:end),[n,n]);
K=s(n+N+1);
options = optimoptions('quadprog','Display','off');

delp_temp=1-2*rand(n,1);
delp_normlized=delp_temp/norm(delp_temp);
price=p0+del_p*rand(1)*delp_normlized;

x_sum=zeros(n,1);
dgamma=zeros(N,1);
for i=1:N
    %     psi=2*exp(-log(2)*norm(price-p)/delta(i))-1;
    %     psi=2*sech(log(2+sqrt(3))*norm(price-p)/delta(i))-1;
    psi=-tanh(h(i)*(norm(price-p)-delta(i)));
    if gamma(i)<1-1e-5 && gamma(i)>1e-5
        dgamma(i)=eta(i)*psi;
    else
        dgamma(i)=eta(i)*psi+eta(i)*max(0,psi*(2*gamma(i)-1))*(1-2*gamma(i));
    end
    x0=x_store(:,i);
    lam= gamma(i)*p+(1-gamma(i))*lam_hat(:,i);
    x_store(:,i)=quadprog(2*a(i)*eye(n),b(i)*ones(n,1)+lam,[],[],ones(n,1)',d(i),zeros(n,1),maxRate(i)*ones(n,1),x0,options);
    x_sum =x_sum+x_store(:,i);
end
t
[xstar,dxs]=x_star(t);
if norm(p0-p)<delta_bar-1e-5
    dp=x_sum-xstar+K*dxs;
else
    dp=x_sum-xstar+K*dxs+((delta_bar/norm(p0-p)-1)/ep)*(p-p0);
end


% if norm(p0-p)<delta_bar-1e-5      %%for hard-adaptive nudge instead of soft
%     dp=x_sum-xstar+K*dxs;
%  else
%     dp=x_sum-xstar+K*dxs-max([0,(x_sum-xstar+K*dxs)'*(p-p0)/(norm(p-p0))])*(p-p0)/(norm(p-p0));
% end


if norm(K,'fro')<K0
    sigma_s=0;
elseif norm(K,'fro')>=K0 && norm(K,'fro')<=2*K0
    sigma_s=sigma*(norm(K,'fro')/K0-1);
else
    sigma_s=sigma;
end


% dK=tau*((x_sum-xstar)*dxs'-sigma_s*K);
dK=tau*((x_sum-xstar)'*dxs-sigma_s*K);
ds=[dp;dgamma;dK(:)];
end

%% x_star
function [x,dx]=x_star(t)
global x_star1 x_star2
g=0.5*(1+cos(3*t));
dg=-0.5*3*sin(3*t);
x=g*x_star1+(1-g)*x_star2;
if nargout>1
    dx=dg*(x_star1-x_star2);
end
end

%% events function to stop simulation
function [position,isterminal,direction] = eventsFcn(t,state)
global n N
position = norm(state(n+1:n+N)-ones(N))<1e-3; % The value that we want to be zero
isterminal = 1;  % Halt integration 
direction = 0;   % The zero can be approached from either direction
end