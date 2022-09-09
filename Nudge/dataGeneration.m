function dataGeneration
N=10; %number of vehicels;
n=24; %size of action
%% EVs' cost 
a_min=0.004;
a_max=0.006;
a = a_min+(a_max-a_min).*rand(N,1); % a_min <= a <= a_max
b_min=0.065;
b_max=0.085;
b = b_min+(b_max-b_min).*rand(N,1); % b_min <= b <= b_max
maxRate_min=8;
maxRate_max=10; 
maxRate=maxRate_min+(maxRate_max-maxRate_min).*rand(N,1); % maxRate_min <= maxRate <= maxRate_max
d_min=25;
d_max=35;
d=d_min+(d_max-d_min).*rand(N,1); % d_min <= d <= d_max
%% EVs' model
delta_min=0.3;
delta_max=0.5;
delta=delta_min+(delta_max-delta_min).*rand(N,1);
eta_min=3;
eta_max=5; 
eta=eta_min+(eta_max-eta_min).*rand(N,1);
h_min=2;
h_max=5;
h=h_min+(h_max-h_min).*rand(N,1);
lam_hat_min=0.1;
lam_hat_max=0.5;
lam_hat=lam_hat_min+(lam_hat_max-lam_hat_min).*rand(n,N); %columns for agents
%% export data
save ('Data/PEVsData.mat','N','n','a','b','maxRate','d','delta','eta','lam_hat','h');
end