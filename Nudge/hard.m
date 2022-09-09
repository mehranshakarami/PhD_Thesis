function ds=hard(t,s)
global  x_store x_star1 h x_star2 del_p p0 delta_bar a b N n delta eta lam_hat d maxRate price price_time;

p=s(1:n);
gamma=s(n+1:n+N);

options = optimoptions('quadprog','Display','off');

price_local=transpose(interp1(price_time,price,t));

x_sum=zeros(n,1);
dgamma=zeros(N,1);
for i=1:N
    %     psi=2*exp(-log(2)*norm(price-p)/delta(i))-1;
    %     psi=2*sech(log(2+sqrt(3))*norm(price-p)/delta(i))-1;
    psi=-tanh(h(i)*(norm(price_local-p)-delta(i)));
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

%%%%%%%%%%%%%%%%%%%%%%%%%%% soft %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% if norm(p0-p)<delta_bar-1e-5
%     dp=x_sum-x_star1;
% else
%     dp=x_sum-x_star1+((delta_bar/norm(p0-p)-1)/ep)*(p-p0);
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%% hard %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if norm(p0-p)<delta_bar-1e-5      
    dp=x_sum-x_star1;
 else
    dp=x_sum-x_star1-max([0,(x_sum-x_star1)'*(p-p0)/(norm(p-p0))])*(p-p0)/(norm(p-p0));
end


ds=[dp;dgamma];
end
