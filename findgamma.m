function [gamma, kappa, zeta] = findgamma(n1, n2, mu_tilda, delta, lambda_1, lambda_2, c_2)
%Finding stationary distribution of Queue2
p_k_n_2 = findp_k_n_2(n1, n2, mu_tilda, delta, lambda_1, lambda_2, c_2);
kappa = 1-p_k_n_2(n1+1);

numer = 0;
for i = 1:(n2)
    numer = numer+p_k_n_2(i);
end

denom = numer;
for i = (n2+1):(n1)
    denom = denom+p_k_n_2(i);
end

zeta = numer;

gamma = (kappa*delta*lambda_1)/(kappa*delta*lambda_1+zeta*lambda_2);
end