function utility_n = findutility_n(type,n) %function to find total expected utility per unit time obtained by the customers in the system
global lambda;
global lambda_1;
global user0;
global user1;
p_k_n = zeros(1,n+1);
d_k_n = zeros(1,n+1);
alpha_k_n = zeros(1,n+1);
d_k_total = 0;
sigma_total=0;
if type==0
    for i=1:n+1
        d_k_n(i) = findd_k(user0,i-1);
        d_k_total = d_k_total+d_k_n(i);
    end
    for i=1:n+1
        p_k_n(i) = d_k_n(i)/d_k_total;
    end
    for i=1:n+1
        alpha_k_n(i) = findalpha(user0,i-1);
    end
    for k=0:n-1
        sigma_total = sigma_total + (p_k_n(k+1)*alpha_k_n(k+1));
    end
    utility_n = lambda*sigma_total;
elseif type==1
    for i=1:n+1
        d_k_n(i) = findd_k(user1,i-1);
        d_k_total = d_k_total + d_k_n(i);
    end
    for i=1:n+1
        p_k_n(i) = d_k_n(i)/d_k_total;
    end
    for i=1:n+1
        alpha_k_n(i) = findalpha(user1,i-1);
    end
    for k=0:n-1
        sigma_total = sigma_total + (p_k_n(k+1)*alpha_k_n(k+1));
    end
    utility_n = lambda_1*sigma_total;    
end
end
