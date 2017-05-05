function d_k = findd_k(type,k) 
global rho;
global rho_1;
global c;
global c_1;
if type==0
    if k<c
        d_k = ((rho*c)^k)/factorial(k);
    else
        d_k = (((rho*c)^c)*(rho^(k-c)))/factorial(c);
    end
elseif type==1
    if k<c_1
        d_k = ((rho_1*c_1)^k)/factorial(k);
    else
        d_k = ((rho_1*c_1)^c_1)*(rho_1^(k-c_1))/factorial(c_1);
    end    
else %type2-hard to implement here because it requires extra input parameters compared to others
end
end
