% Dual pricing is modeled by considering two separate queues for the two
% different types of users.
% c1 spots for Type1 users (those with parking time < 1hr)
% c2 spots for Type2 users (those with parking time > 1hr)
% All prices are in $/hr

%Defing all the parameters required to find social utility
c=30; %total no. of parking spots

% parking spots in lot 1 and 2
c_1 = 10;
c_2 = 20;
assert(c == c_1+c_2);

% service rate in single-queue case. Value in Ratcliff paper: 1/2
mu = 1/2;

% mean service time of users in single-queue case
T = 1/mu;

% service rate for type 1 and 2 vehicles
mu_1 = mu*(1-exp(-mu*T))/(1-exp(-mu*T)-mu*T*exp(-mu*T));
mu_2 = mu/2;

% arrivals in single queue (per hr)
lambda = 60/5; 

% arrival rate of type 1 and 2vehicles
lambda_1 = lambda*(1-exp(-mu*T)); 
lambda_2 = lambda*exp(-mu*T); %rate of poisson process for Type2

%single queue initial parking price (note: why do we need this?)
P = 5; 

% parking price for lots 1 and 2 (note: per hour?)
P_1 = 2; 
P_2 = 5; 

% Waiting cost
P_w = 48; %price for waiting in the queue

% Reward for parking, single queue
R = 75; %reward for parking

% Create for indexing structures
type_1_idx = 1;
lot_1_idx = 1;
type_2_idx = 2;
lot_2_idx = 2;
original_idx = 3;


% Create an array of structures describing type 1 and 2 users
type_1 = struct('reward', R, 'waiting_cost', P_w, 'arrival_rate', ...
    lambda_1, 'departure_rate', mu_1);
type_2 = struct('reward', R, 'waiting_cost', P_w, 'arrival_rate', ...
    lambda_2, 'departure_rate', mu_2);
original = struct('reward', R, 'waiting_cost', P_w, 'arrival_rate',...
    lambda, 'departure_rate', mu);
users = [type_1, type_2, original];

% Create an array of structures describing lots 1 and 2
lot_1 = struct('num_spots', c_1, 'parking_cost', P_1);
lot_2 = struct('num_spots', c_2, 'parking_cost', P_2);
original = struct('num_spots', c, 'parking_cost', P);
lots = [lot_1, lot_2, original];

% Erlang loss for lot 1 (note: shouldn't this be determined on the fly?)
delta=0.0553; %0.0032; %0.1 %changing it to p_nb1(nb1) 

%initial proportion of type 1 users in lot 2
gamma = 0; 

% average service time from lot 2(note: fix this)
mu_tilda = gamma*mu_1 + (1-gamma)*mu_2;

% average arrical time to lot 2
lambda_tilda = delta*lambda_1 + lambda_2;

% average traffic intensity at lot 2
rho_tilda = lambda_tilda/(c_2*mu_tilda);

% Todo: describe rho hat
rho_hat = delta*lambda_1/(c_2*mu_tilda);


%balking level of type1 users
n_b1 = floor((R*mu_1*c_1 + P_w*c_1 - P_1*c_1)/P_w); 
%balking level of type2 users
n_b2 = floor((R*mu_2*mu_tilda*c_2 + P_w*mu_2*c_2 - P_2*mu_tilda*c_2)/(P_w*mu_2));
%balking level in single pricing model assuming the price to be P_2
n_b = floor((R*mu*c + P_w*c - P*c)/P_w); 

%max number of customers in the system (for simulation purposes)
n_max = 120; 

% array to store single price utilityand queue 1 utility for different values of n
utility = zeros(1,n_max); 
utility_1 = zeros(1,n_max);
%utility_2 = zeros(1,n);

for i=1:n_max
    utility(i) = findutility_n(i, users(original_idx),...
        lots(original_idx));
    utility_1(i) = findutility_n(i, users(type_1_idx),...
        lots(lot_1_idx));
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

% probability of acceptance of cars of types 1 and 2 arriving at queue 2 
% assigning some initial guess of the correct values
kappa = 0; 
zeta = 0.8; 


% balking level at lot 1
n_hat=16;

% desired balking level of type 2 users at lot 2
n2= 19;

% resulting balking level of type 1 users at lot 2 (initial guess)
n1=1;

delta = finddelta(n_hat, users(type_1_idx), lots(lot_1_idx));

% Update our overall arrival rate for lot 2 based on Erlang loss from lot 1
lambda_tilda = delta*lambda_1 + lambda_2;

%Implementation of "Iterative price selection"
for z=1:10

    mu_tilda = gamma*mu_1+(1-gamma)*mu_2;
    P_2= (R*mu_2*mu_tilda*c_2 + P_w*mu_2*c_2 - n2*P_w*mu_2 )/(mu_tilda*c_2);
    n1 = floor((R*mu_1*mu_tilda*c_2 + P_w*c_2*mu_1 - P_2*mu_tilda*c_2)/(P_w*mu_1));
    [gamma, kappa, zeta] = findgamma(n1, n2, mu_tilda, delta, lambda_1, lambda_2, c_2);
    disp(mu_tilda);
    disp(P_2);
    disp(n1);
    disp(gamma);
end

%total utility; first include total utility from lot 1
utility_total = findutility_n(n_hat, users(type_1_idx), lots(lot_1_idx));

%disp(utility_total);

% Note: the utility below is incorrect
p_k_n_2 = findp_k_n_2(n1, n2, mu_tilda, delta, lambda_1, lambda_2, c_2);
for i=1:n2+1
    utility_total=utility_total+(lambda_tilda)*(p_k_n_2(i))*...
        (findalpha(i-1, users(type_2_idx), lots(lot_2_idx))); 
end
for i=n2+2:n1+1
    utility_total=utility_total+(delta*lambda_1)*(p_k_n_2(i))*...
        (findalpha(i-1, users(type_2_idx), lots(lot_2_idx)));
end

%{
figure(1);
bar(utility);
title('Total expected utility vs n');
ylabel('utility');

figure(2);
bar(utility_1);
title('Total expected utility of Type1 users vs n');
ylabel('utility1');

function utility_n_2 = findutility_n_2(n) %Implement this in findutility_n.m

%}