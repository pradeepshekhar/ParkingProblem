function alpha = findalpha_21(k, user_params, lot_params, mu_tilda) 

c = lot_params.num_spots;
R_1 = user_params.reward;
P_w = user_params.waiting_cost;

if k<c
    alpha = R_1;
else
    alpha = R_1 - (P_w*(k-c+1))/(c*mu_tilda);
end