% Dual pricing is modeled by considering two separate queues for the two
% different types of users.
% c1 spots for Type1 users (those with parking time < 1hr)
% c2 spots for Type2 users (those with parking time > 1hr)
% All prices are in $/hr
%Defing all the parameters required to find social utility
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

global n; %max number of customers in the system for simulation purposes
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
    utility(i) = findutility_n(i);
    utility_1(i) = findutility_n_1(i);
    %utility_2(i) = findutility_n_2(i);
end

global n2;
n2= 21; %desired n2
global n1;
n1=1;

%Implementation of "Iterative price selection"
for z=1:5

    mu_tilda = gamma*mu_1+(1-gamma)*mu_2;
    P_2= floor((R*mu_2*mu_tilda*c_2 + P_w*mu_2*c_2 - n2*P_w*mu_2 )/(mu_tilda*c_2));
    n1 = floor((R*mu_1*mu_tilda*c_2 + P_w*c_2*mu_1 - P_2*mu_tilda*c_2)/(P_w*mu_1));
    gamma = findgamma(n1,n2,mu_tilda);
    disp(mu_tilda);
    disp(P_2);
    disp(n1);
    disp(gamma);
end

[U_sw,f] = max(utility);
[U_sw_1,g] = max(utility_1);
%[U_sw_2,h] = max(utility_2);
n_so=f-1;
n_so_1=g-1;
%n_so_2=h-1;
%U_sw_new = U_sw_1+U_sw_2;
%U_sw_user = U_sw/n_so;
%U_sw_new_user = U_sw_new/(n_so_1+n_so_2);


%Change P_1 to get n_b1=n_so_1
P_1 = floor((R*mu_1*c_1 +P_w*c_1-n_so_1*P_w)/c_1);
n_b1 = floor((R*mu_1*c_1 + P_w*c_1 - P_1*c_1)/P_w); %new balking level with adjusted price
%fiding delta
p_k_n_1 = zeros(1,n_b1+1);
d_k_n_1 = zeros(1,n_b1+1);
d_k_total_1 = 0;
for i=1:n_b1+1
    d_k_n_1(i) = findd_k_1(i-1);
    d_k_total_1 = d_k_total_1 + d_k_n_1(i);
end
for i=1:n_b1+1
    p_k_n_1(i) = d_k_n_1(i)/d_k_total_1;
end
delta = p_k_n_1(n_b1+1);

%{
figure(1);
bar(utility);
title('Total expected utility vs n');
ylabel('utility');

figure(2);
bar(utility_1);
title('Total expected utility of Type1 users vs n');
ylabel('utility1');

figure(3);
bar(utility_2);
title('Total expected utility of Type2 users vs n');
ylabel('utility2');


%}

function utility_n = findutility_n(n) %function to find total expected utility per unit time obtained by the customers in the system
global lambda;
p_k_n = zeros(1,n+1);
d_k_n = zeros(1,n+1);
beta_k_n = zeros(1,n+1);
d_k_total = 0;
sigma_total=0;
for i=1:n+1
    d_k_n(i) = findd_k(i-1);
    d_k_total = d_k_total+d_k_n(i);
end
for i=1:n+1
    p_k_n(i) = d_k_n(i)/d_k_total;
end
for i=1:n+1
    beta_k_n(i) = findbeta(i-1);
end
for k=0:n-1
    sigma_total = sigma_total + (p_k_n(k+1)*beta_k_n(k+1));
end
utility_n = lambda*sigma_total;
end

function utility_n_1 = findutility_n_1(n) %function to find total expected utility per unit time obtained by the customers in the system
global lambda_1;
p_k_n_1 = zeros(1,n+1);
d_k_n_1 = zeros(1,n+1);
beta_k_n_1 = zeros(1,n+1);
d_k_total_1 = 0;
sigma_total_1 = 0;
for i=1:n+1
    d_k_n_1(i) = findd_k_1(i-1);
    d_k_total_1 = d_k_total_1 + d_k_n_1(i);
end
for i=1:n+1
    p_k_n_1(i) = d_k_n_1(i)/d_k_total_1;
end
for i=1:n+1
    beta_k_n_1(i) = findbeta_1(i-1);
end
for k=0:n-1
    sigma_total_1 = sigma_total_1 + (p_k_n_1(k+1)*beta_k_n_1(k+1));
end
utility_n_1 = lambda_1*sigma_total_1;
end
%{
function utility_n_2 = findutility_n_2(n) %function to find total expected utility per unit time obtained by the customers in the system
global lambda_2;
p_k_n_2 = zeros(1,n+1);
d_k_n_2 = zeros(1,n+1);
beta_k_n_2 = zeros(1,n+1);
d_k_total_2 = 0;
sigma_total_2 = 0;
for i=1:n+1
    d_k_n_2(i) = findd_k_2(i-1);
    d_k_total_2 = d_k_total_2 + d_k_n_2(i);
end
for i=1:n+1
    p_k_n_2(i) = d_k_n_2(i)/d_k_total_2;
end
for i=1:n+1
    beta_k_n_2(i) = findbeta_2(i-1);
end
for k=0:n-1
    sigma_total_2 = sigma_total_2 + (p_k_n_2(k+1)*beta_k_n_2(k+1));
end
utility_n_2 = lambda_2*sigma_total_2;
end
%}
function d_k = findd_k(k) 
global rho;
global c;
if k<c
    d_k = ((rho*c)^k)/factorial(k);
else
    d_k = (((rho*c)^c)*(rho^(k-c)))/factorial(c);
end
end

function d_k_1 = findd_k_1(k) 
global rho_1;
global c_1;
if k<c_1
    d_k_1 = ((rho_1*c_1)^k)/factorial(k);
else
    d_k_1 = ((rho_1*c_1)^c_1)*(rho_1^(k-c_1))/factorial(c_1);
end
end

function beta = findbeta(k) %expected utility of a user who enters the system in state k
global R;
global P_w;
%global P;
global mu;
global c;
if k==0
    beta = 0;
elseif k<c
    beta = R ;%- P/mu;
else
    beta = R - (P_w*(k-c+1))/(c*mu);% - P/mu ;%beta = R - (P_w*(k+1))/(c*mu) - P/mu ;
end
end

function beta_1 = findbeta_1(k) %expected utility of a Type1 user who enters the system in state k
global R;
global P_w;
%global P_1;
global mu_1;
global c_1;
if k==0
    beta_1 = 0;
elseif k<c_1
    beta_1 = R ;%- P_1/mu_1;
else
    beta_1 = R - (P_w*(k-c_1+1))/(c_1*mu_1);% - P_1/mu_1 ;
end
end

function d_k_2 = findd_k_2(k,n2,mu_tilda)
global delta;
global lambda_tilda;
global lambda_2;
global c_2;
rho_tilda=lambda_tilda/(c_2*mu_tilda);
rho_hat = (delta*lambda_2)/(c_2*mu_tilda);
if k<c_2
    d_k_2 = ((rho_tilda*c_2)^k)/factorial(k);
elseif k<n2
    d_k_2 = ((rho_tilda*c_2)^c_2)*((rho_tilda)^(k-c_2))/factorial(c_2);
else
    d_k_2 = ((rho_tilda*c_2)^n2)*((rho_hat)^(k-n2))/(factorial(c_2)*(c_2^(n2-c_2)));
end
end

function alpha_k_2 = findalpha_k_2(k)
global R;
global P_w;
global mu_tilda;
global c_2;
if k==0
    alpha_k_2 = 0;
elseif k<c_2
    alpha_k_2 = R ;
else
    alpha_k_2 = R - (P_w*(k-c_2+1))/(c_2*mu_tilda);
end 
end

function gamma = findgamma(n1,n2,mu_tilda)
global kappa;
p_k_n_2=findp_k_n_2(n1,n2,mu_tilda);
kappa = 1-p_k_n_2(n1+1);
global delta;
global lambda_1;
global lambda_2;
global zeta;
zeta=0;
for i=1:n2
    zeta = zeta+p_k_n_2(i);
end
gamma = (kappa*delta*lambda_1)/(kappa*delta*lambda_1+zeta*lambda_2);
end

function p_k_n_2 = findp_k_n_2(n1,n2,mu_tilda)
p_k_n_2=zeros(1,n1+1);
d_k_2=zeros(1,n1+1);
D_n_2=0;
for i=1:n1+1
    d_k_2(i) = findd_k_2(i-1,n2,mu_tilda);
    D_n_2 = D_n_2+d_k_2(i);
end
for i=1:n1+1
    p_k_n_2(i) = d_k_2(i)/D_n_2;
end
end