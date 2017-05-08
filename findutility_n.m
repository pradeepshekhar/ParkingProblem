%Function to find Utilities of Single priced queue and Queue 1
function utility_n = findutility_n(n, user_params, lot_params) 
% function to find total expected utility per unit time obtained by the
% customers in the system

lambda = user_params.arrival_rate;


% Initialize the variables we'll need
p_k_n = zeros(1,n+1);
d_k_n = zeros(1,n+1);
alpha_k_n = zeros(1,n+1);
d_k_total = 0;
sigma_total=0;

for i=1:n+1
    d_k_n(i) = findd_k(i-1, user_params, lot_params);
    d_k_total = d_k_total+d_k_n(i);
end
for i=1:n+1
    p_k_n(i) = d_k_n(i)/d_k_total;
end
for i=1:n+1
    alpha_k_n(i) = findalpha(i-1, user_params, lot_params);
end
for k=0:n-1
    sigma_total = sigma_total + (p_k_n(k+1)*alpha_k_n(k+1));
end
utility_n = lambda*sigma_total;

end
