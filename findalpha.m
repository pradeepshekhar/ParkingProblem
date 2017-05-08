function alpha = findalpha(k, user_params, lot_params) 
% expected utility of a user who enters the system in state k, not counting
% the price of parking

c = lot_params.num_spots;
R = user_params.reward;
P_w = user_params.waiting_cost;
mu = user_params.departure_rate;

%if k==0 %  Note(DL): Why??
%    alpha = 0;
if k<c
    alpha = R ;%- P/mu;
else
    alpha = R - (P_w*(k-c+1))/(c*mu);    % - P/mu ;%alpha = R - (P_w*(k+1))/(c*mu) - P/mu ;
end