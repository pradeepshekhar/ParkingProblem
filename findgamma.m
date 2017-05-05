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