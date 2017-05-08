%Fuction for finding the stationary distribution of Queue 2
function p_k_n_2 = findp_k_n_2(n1, n2, mu_tilda, delta, lambda_1, lambda_2, c_2)
p_k_n_2=zeros(1,n1+1);
d_k_2=zeros(1,n1+1);
D_n_2=0;
for i=1:n1+1
    d_k_2(i) = findd_k_2(i-1,n2,mu_tilda, delta, lambda_1, lambda_2, c_2);
    D_n_2 = D_n_2+d_k_2(i);
end
for i=1:n1+1
    p_k_n_2(i) = d_k_2(i)/D_n_2;
end
end

function d_k_2 = findd_k_2(k, n2, mu_tilda, delta, lambda_1, lambda_2, c_2)

lambda_tilda = delta*lambda_1 + lambda_2;
rho_tilda=lambda_tilda/(c_2*mu_tilda);
rho_hat = (delta*lambda_1)/(c_2*mu_tilda);
if k<c_2
    d_k_2 = ((rho_tilda*c_2)^k)/factorial(k);
elseif k<n2
    d_k_2 = ((rho_tilda*c_2)^c_2)*((rho_tilda)^(k-c_2))/factorial(c_2);
else
    d_k_2 = ((rho_tilda*c_2)^n2)*((rho_hat)^(k-n2))/(factorial(c_2)*(c_2^(n2-c_2)));
end
end