% Dual pricing is modeled by considering two separate queues for the two
% different types of users.
% c1 spots for Type1 users (those with parking time < 1hr)
% c2 spots for Type2 users (those with parking time > 1hr)
% All prices are in $/hr
%Defing all the parameters required to find social utility
global user0; %variable defining the type of user - used when calling functions
user0 = 0;
global user1;
user1 = 1;
global user2;
user2 = 2;
global c;
c=30; %total no. of parking spots
global c_1;
c_1 = 8;
global c_2;
c_2 = 22;
global mu;
mu = 1/2; % 1/2 used in paper %service time parameter for single queue scenario
global T; %mean service time of users
T = 1/mu;
global mu_1;
mu_1 = mu*(1-exp(-mu*T))/(1-exp(-mu*T)-mu*T*exp(-mu*T));
global mu_2;
mu_2 = mu/2;
global lambda;
lambda = 60/5; %rate of poisson process arrivals in single queue (per hr)
global lambda_1;
lambda_1 = lambda*(1-exp(-mu*T)); %rate of poisson process for Type1
global lambda_2;
lambda_2 = lambda*exp(-mu*T); %rate of poisson process for Type2
global P;
P = 5; %initial parking price
global P_1;
P_1 = 2; %parking price for Type1 users 
global P_2;
P_2 = 5; %parking price for Type2 users
global P_w;
P_w = 48; %price for waiting in the queue
global R;
R = 75; %reward for parking
global rho;
rho = lambda/(c*mu); %traffic intensity of single queue
global rho_1;
rho_1 = lambda_1/(c_1*mu_1); %traffic intensity of Type1 queue
global rho_2;
rho_2 = lambda_2/(c_2*mu_2); %traffic intensity of Type1 queue
global delta;
delta=0.0553; %0.0032; %0.1 %changing it to p_nb1(nb1) 
global gamma;
gamma = 0; %proportion of Type1 users in Queue2
global mu_tilda;
mu_tilda = gamma*mu_1 + (1-gamma)*mu_2;
global lambda_tilda;
lambda_tilda = delta*lambda_1 + lambda_2;
global rho_tilda;
rho_tilda = lambda_tilda/(c_2*mu_tilda);
global rho_hat;
rho_hat = delta*lambda_1/(c_2*mu_tilda);
n_b1 = floor((R*mu_1*c_1 + P_w*c_1 - P_1*c_1)/P_w); %balking level of type1 users
n_b2 = floor((R*mu_2*mu_tilda*c_2 + P_w*mu_2*c_2 - P_2*mu_tilda*c_2)/(P_w*mu_2)); %balking level of type2 users
n_b = floor((R*mu*c + P_w*c - P*c)/P_w); %balking level in single pricing model assuming the price to be P_2

global n; %max number of customers in the system %%for simulation purposes
n = 120;

%Defining the acceptance rates
global kappa; %acceptance rate for type 1 cars arriving at queue 2
global zeta; %acceptance rate for type 2 cars arriving at queue 2
global psi; %overall acceptance rate for cars arriving at queue 2
%Assigning some initial values to start with
kappa =0;
zeta=0.8;
psi=0;
utility = zeros(1,n); % array to store single price utility for different values of n
utility_1 = zeros(1,n); %array to store queue-1 utility for different values of n
%utility_2 = zeros(1,n);

for i=1:n
    utility(i) = findutility_n(user0,i);
    utility_1(i) = findutility_n(user1,i);
    %utility_2(i) = findutility_n_2(i);
end
[U_sw,f] = max(utility);
[U_sw_1,g] = max(utility_1);
%[U_sw_2,h] = max(utility_2);
n_so=f-1;
n_so_1=g-1;
%n_so_2=h-1;
%U_sw_new = U_sw_1+U_sw_2;

%Change P,P_1 to get n_b1=n_so_1 and n_b=n_so
P = floor((R*mu*c +P_w*c-n_so*P_w)/c);
P_1 = floor((R*mu_1*c_1 +P_w*c_1-n_so_1*P_w)/c_1);
n_b = floor((R*mu*c + P_w*c - P*c)/P_w);
n_b1 = floor((R*mu_1*c_1 + P_w*c_1 - P_1*c_1)/P_w); %new balking level with adjusted price

global nhat;
nhat=1;
global n2;
n2= 21; %desired n2
global n1;
n1=1;
%Implementation of "Iterative price selection"
for z=1:5

    mu_tilda = gamma*mu_1+(1-gamma)*mu_2;
    P_2= (R*mu_2*mu_tilda*c_2 + P_w*mu_2*c_2 - n2*P_w*mu_2 )/(mu_tilda*c_2);
    n1 = floor((R*mu_1*mu_tilda*c_2 + P_w*c_2*mu_1 - P_2*mu_tilda*c_2)/(P_w*mu_1));
    gamma = findgamma(n1,n2,mu_tilda);
    disp(mu_tilda);
    disp(P_2);
    disp(n1);
    disp(gamma);
end

figure(1);
bar(utility);
title('Total expected utility vs n');
ylabel('utility');

figure(2);
bar(utility_1);
title('Total expected utility of Type1 users vs n');
ylabel('utility1');
%{
figure(3);
bar(utility_2);
title('Total expected utility of Type2 users vs n');
ylabel('utility2');

function utility_n_2 = findutility_n_2(n) %Implement this in findutility_n.m

%}