function d_k = findd_k(k, user_params, lot_params) 

c = lot_params.num_spots;
rho = user_params.arrival_rate/(c*user_params.departure_rate);


if k<c
    d_k = ((rho*c)^k)/factorial(k);
else
    d_k = (((rho*c)^c)*(rho^(k-c)))/factorial(c);
end
 
%type2-hard to implement here because it requires extra input parameters compared to others
