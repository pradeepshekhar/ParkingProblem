%Expected utility of a user who enters the system in state k, excluding 
%the price of parking (for Single priced queue and Queue 1)
function alpha = findalpha(k, user_params, lot_params) 

c = lot_params.num_spots;
R = user_params.reward;
P_w = user_params.waiting_cost;
mu = user_params.departure_rate;

if k<c
    alpha = R ;%- P/mu;
else
    alpha = R - (P_w*(k-c+1))/(c*mu);    % - P/mu ;%alpha = R - (P_w*(k+1))/(c*mu) - P/mu ;
end