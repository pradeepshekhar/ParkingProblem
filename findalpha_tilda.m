function alpha = findalpha_tilda(k, user_params, lot_params, mu_tilda, gamma) 

c = lot_params.num_spots;
R_1 = user_params.reward_1;
R_2 = user_params.reward_2;
P_w = user_params.waiting_cost;

if k<c
    alpha = gamma*R_1+(1-gamma)*R_2;
else
    alpha = gamma*R_1+(1-gamma)*R_2 - (P_w*(k-c+1))/(c*mu_tilda);
end