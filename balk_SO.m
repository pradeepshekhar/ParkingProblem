% Dual pricing is modeled by considering two separate queues for the two
% different types of users.
% c1 spots for Type1 users (those with parking time < 1hr)
% c2 spots for Type2 users (those with parking time > 1hr)
% All prices are in $/hr
%Defing all the parameters required to find social utility
global c;
c=30; %total no. of parking spots
global c_1;
c_1 = 15; 
global c_2;
c_2 = 15;
global mu;
mu = 1; % 1/2 used in paper %service time parameter for single queue scenario
global mu_1;
mu_1 = mu*(1-exp(-mu))/(1-exp(-mu)-mu*exp(-mu));
global mu_2;
mu_2 = mu/(1+mu);
global lambda;
lambda = 140/5; %rate of poisson process arrivals in single queue (per hr)
global lambda_1;
lambda_1 = lambda*(1-exp(-mu)); %rate of poisson process for Type1
global lambda_2;
lambda_2 = lambda*exp(-mu); %rate of poisson process for Type2
global P_1;
P_1 = 2; %parking price for Type1 users 
global P_2;
P_2 = 5; %parking price for Type2 users
global P_w;
P_w = 30; %price for waiting in the queue
global R;
R = 75; %reward for parking
global rho;
rho = lambda/(c*mu); %traffic intensity of single queue
global rho_1;
rho_1 = lambda_1/(c_1*mu_1); %traffic intensity of Type1 queue
global rho_2;
rho_2 = lambda_2/(c_2*mu_2); %traffic intensity of Type1 queue
n_b1 = floor((R*mu_1*c_1 + P_w*c_1 - P_1*c_1)/P_w); %balking level of type1 users
n_b2 = floor((R*mu_2*c_2 + P_w*c_2 - P_2*c_2)/P_w); %balking level of type2 users
n_b = floor((R*mu*c + P_w*c - P_2*c)/P_w); %balking level in single pricing model assuming the price to be P_2

global n; %max number of customers in the system
n = 60;

utility = zeros(1,n);
utility_1 = zeros(1,n);
utility_2 = zeros(1,n);

for i=1:n
    utility(i) = findutility_n(i);
    utility_1(i) = findutility_n_1(i);
    utility_2(i) = findutility_n_2(i);
end

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

%{


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
    %{
    sigma_beta = 0;
    for l=0:k
        sigma_beta = sigma_beta + beta_k_n(l+1);
    end
    sigma_total = sigma_total + (p_k_n(k+1)*sigma_beta);
    %}
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


function d_k = findd_k(k) 
global rho;
global c;
if k<c
    d_k = ((rho*c)^k)/factorial(k);
else
    d_k = ((rho*c)^c)*(rho^(k-c))/factorial(c);
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

function d_k_2 = findd_k_2(k) 
global rho_2;
global c_2;
if k<c_2
    d_k_2 = ((rho_2*c_2)^k)/factorial(k);
else
    d_k_2 = ((rho_2*c_2)^c_2)*(rho_2^(k-c_2))/factorial(c_2);
end
end

function beta = findbeta(k) %expected utility of a user who enters the system in state k
global R;
global P_w;
global P_2;
global mu;
global c;
if k==0
    beta = 0;
elseif k<c
    beta = R - P_2/mu;
else
    beta = R - (P_w*(k-c+1))/(c*mu) - P_2/mu ;%beta = R - (P_w*(k+1))/(c*mu) - P_2/mu ;
end
end

function beta_1 = findbeta_1(k) %expected utility of a Type1 user who enters the system in state k
global R;
global P_w;
global P_1;
global mu_1;
global c_1;
if k==0
    beta_1 = 0;
elseif k<c_1
    beta_1 = R - P_1/mu_1;
else
    beta_1 = R - (P_w*(k-c_1+1))/(c_1*mu_1) - P_1/mu_1 ;
end
end

function beta_2 = findbeta_2(k) %expected utility of a Type2 user who enters the system in state k
global R;
global P_w;
global P_2;
global mu_2;
global c_2;
if k==0
    beta_2 = 0;
elseif k<c_2
    beta_2 = R - P_2/mu_2;
else
    beta_2 = R - (P_w*(k-c_2+1))/(c_2*mu_2) - P_2/mu_2 ;
end    
end