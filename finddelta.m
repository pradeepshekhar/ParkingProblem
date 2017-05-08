%Fuction to calculate the rate at which type 1 cars will defect from lot 1
function delta = finddelta(n_hat, user_params, lot_params)

p_k_n = zeros(1,n_hat+1);
d_k_n = zeros(1,n_hat+1);
d_k_total_1 = 0;

for i=1:n_hat+1
    d_k_n(i) = findd_k(i-1, user_params, lot_params);
    d_k_total_1 = d_k_total_1 + d_k_n(i);
end
for i=1:n_hat+1
    p_k_n(i) = d_k_n(i)/d_k_total_1;
end
delta = p_k_n(n_hat+1);
end
