% Dual pricing is modeled by considering two separate queues for the two
% different types of users.
% c1 spots for Type1 users (those with parking time < 1hr)
% c2 spots for Type2 users (those with parking time > 1hr)
% All prices are in $/hr

%Defining all the parameters required to find social utility
c=30; %total no. of parking spots

% parking spots in lot 1 and 2
c_1 = 10;
c_2 = 20;
assert(c == c_1+c_2);

% service rate in single-queue case. Value in Ratliff paper: 1/2
mu = 1/2;

% mean service time of users in single-queue case
T = 1;%/mu;

% service rate for type 1 and 2 vehicles
mu_1 = mu*(1-exp(-mu*T))/(1-exp(-mu*T)-mu*T*exp(-mu*T));
mu_2 = mu/(1+mu*T);

% arrivals in single queue (per hr)
lambda = 60/5;

% arrival rate of type 1 and 2 vehicles
lambda_1 = lambda*(1-exp(-mu*T)); 
lambda_2 = lambda*exp(-mu*T); %rate of poisson process for Type2

%Traffic intensity of single queue
rho = lambda/(c*mu); 
%Traffic intensity of Type1 queue
rho_1 = lambda_1/(c_1*mu_1); 

%single queue initial parking price 
%not required as it is determined based on optimality later
P = 0;

% parking price for lots 1 and 2 (note: per hour?)
P_1 = 0; 
P_2 = 0; 

% Waiting cost
P_w = 48; 

% Reward for parking: single queue, Type1 and Type2
R = 75;
%R_1 = 50; R_2 = 100;
R_1 = R*lambda*(mu_2)/(lambda_1*mu_2+lambda_2*mu_1); R_2 = R*lambda*(mu_1)/(lambda_1*mu_2+lambda_2*mu_1); 

% Create for indexing structures
type_1_idx = 1;
lot_1_idx = 1;
type_2_idx = 2;
lot_2_idx = 2;
original_idx = 3;
mixed_idx = 4;

% Create an array of structures describing type 1 and 2 users
type_1 = struct('reward', R_1, 'waiting_cost', P_w, 'arrival_rate', ...
    lambda_1, 'departure_rate', mu_1);
type_2 = struct('reward', R_2, 'waiting_cost', P_w, 'arrival_rate', ...
    lambda_2, 'departure_rate', mu_2);
original = struct('reward', R, 'waiting_cost', P_w, 'arrival_rate',...
    lambda, 'departure_rate', mu);
mixed = struct('reward_1',R_1,'reward_2',R_2,'waiting_cost', P_w); %note: implement later
users = [type_1, type_2, original];

% Create an array of structures describing lots 1 and 2
% Note - P, P_1 and P_2 are changed dynamically to get the required balking
% level. So, why create a structure with them? 
lot_1 = struct('num_spots', c_1, 'parking_cost', P_1);
lot_2 = struct('num_spots', c_2, 'parking_cost', P_2);
original = struct('num_spots', c, 'parking_cost', P);
lots = [lot_1, lot_2, original];

% Erlang loss for lot 1 (note: shouldn't this be determined on the fly?)
delta=0; %initialized to 0. Will be changed to p_nb1(nb1) later on the fly

%initial proportion of type 1 users in lot 2
gamma = 0;

% average service time from lot 2
mu_tilda = (mu_1*mu_2)/(gamma*mu_2 + (1-gamma)*mu_1);

% average arrival time to lot 2
lambda_tilda = delta*lambda_1 + lambda_2;

% average traffic intensity at lot 2
rho_tilda = lambda_tilda/(c_2*mu_tilda);

% is this ok?- effective load of type1 users at lot2
rho_hat = delta*lambda_1/(c_2*mu_tilda);


%balking level of type1 users
%n_b1 = floor((R*mu_1*c_1 + P_w*c_1 - P_1*c_1)/P_w); 
%balking level of type2 users
n_b2 = floor((R*mu_2*mu_tilda*c_2 + P_w*mu_2*c_2 - P_2*mu_tilda*c_2)/(P_w*mu_2));
%balking level in single pricing model
%n_b = floor((R*mu*c + P_w*c - P*c)/P_w); 

%max number of customers in the system (for simulation purposes)
n_max = 60; 

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
[U_sw,n_so] = max(utility);
[U_sw_1,n_so_1] = max(utility_1);
%[U_sw_2,n_so_2] = max(utility_2);
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
n_hat=0;

% desired balking level of type 2 users at lot 2
n2=0;

n_hat_start = 1;
n2_start = 1;
n_hat_range = 40;
n2_range = 40;

%Matrix to store uitilities for all options of (n_hat,n2)
utility_total = zeros(n_hat_range+1,n2_range+1);
utility_11 = zeros(n_hat_range+1,n2_range+1);
utility_21 = zeros(n_hat_range+1,n2_range+1);
utility_22 = zeros(n_hat_range+1,n2_range+1);
P_2 = zeros(n_hat_range+1,n2_range+1);
n1 = zeros(n_hat_range+1,n2_range+1);

%Cycling through different pairs of (n_hat,n2) and finding utility
for n_hat=n_hat_start:n_hat_start+n_hat_range
    for n2=n2_start:n2_start+n2_range
        x = n_hat-n_hat_start+1;
        y = n2-n2_start+1;
        
        % resulting balking level of type 1 users at lot 2 (initial guess)
        n1(x,y)=5;

        %Delta found using n_hat
        delta = finddelta(n_hat, users(type_1_idx), lots(lot_1_idx));

        % Update our overall arrival rate for lot 2 based on Erlang loss from lot 1
        lambda_tilda = delta*lambda_1 + lambda_2;
        rho_hat = delta*lambda_1/(c_2*mu_tilda);

        %Implementation of "Iterative price selection"
        for z=1:10
            
            mu_tilda = (mu_1*mu_2)/(gamma*mu_2 + (1-gamma)*mu_1);
            P_2(x,y) = (R_2*mu_2*mu_tilda*c_2 + P_w*mu_2*c_2 - n2*P_w*mu_2 )/(mu_tilda*c_2);
            %note: is n1 formulation correct? Observed n1<n2.. so no gain in utility!
            n1(x,y) = floor((R_1*mu_tilda*c_2 + P_w*c_2 - (P_2(x,y))*c_2*mu_tilda)/(P_w));
            %n1(x,y) = floor((R_1*mu_1*c_2 + P_w*c_2 - (P_2(x,y))*c_2)/(P_w));
            %n1(x,y) = floor((R_1*mu_1*mu_tilda*c_2 + P_w*c_2*mu_1 - (P_2(x,y))*mu_tilda*c_2)/(P_w*mu_1));
            [gamma, kappa, zeta] = findgamma(n1(x,y), n2, mu_tilda, delta, lambda_1, lambda_2, c_2);
            %{ 
            disp(mu_tilda);
            disp(P_2);
            disp(n1);
            disp(gamma);
            %}
        end

        %total utility; first include total utility from lot 1
        utility_total(x,y) = findutility_n(n_hat, users(type_1_idx), lots(lot_1_idx));
        utility_11(x,y) = utility_total(x,y);
        % Note: the utility below is incorrect?
        p_k_n_2 = findp_k_n_2(n1(x,y), n2, mu_tilda, delta, lambda_1, lambda_2, c_2);
        for i=1:n2
            utility_total(x,y)=utility_total(x,y)+(lambda_tilda)*(p_k_n_2(i))*...
                (findalpha_tilda(i-1, mixed, lots(lot_2_idx), mu_tilda, gamma)); 
        end
        utility_22(x,y) = utility_total(x,y)-utility_11(x,y);
        % Note: alpha used below is not correct?
        for i=n2+1:n1(x,y)
            utility_total(x,y)=utility_total(x,y)+(delta*lambda_1)*(p_k_n_2(i))*...
                (findalpha_21(i-1, users(type_1_idx), lots(lot_2_idx), mu_tilda));
        end
        utility_21(x,y) = utility_total(x,y)-utility_11(x,y)-utility_22(x,y);
        gamma=0;
    end
end

[U_sw_total_max] = max(utility_total(:));
[row, col] = find(utility_total == U_sw_total_max);
U_11_max = utility_11(row,col);
U_21_max = utility_21(row,col);
U_22_max = utility_22(row,col);
P_2_max = P_2(row,col);
n_hat_max = n_hat_start+row-1;
n1_max = n1(row,col);
n2_max = n2_start+col-1;
%change P_1 to achieve n_b1 = n_hat*
P_1 = floor((R_1*mu_1*c_1 +P_w*c_1-n_hat_max*P_w)/c_1);
figure(1);
mesh(utility_total);

%{
figure(2);
bar(utility);
title('Total expected utility vs n');
ylabel('utility');
%}