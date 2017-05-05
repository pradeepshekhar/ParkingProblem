function alpha = findalpha(type,k) %expected utility of a user who enters the system in state k
global R;
global P_w;
global mu;
global mu_1;
global mu_tilda;
global c;
global c_1;
global c_2;
if type==0
    if k==0
        alpha = 0;
    elseif k<c
        alpha = R ;%- P/mu;
    else
        alpha = R - (P_w*(k-c+1))/(c*mu);    % - P/mu ;%alpha = R - (P_w*(k+1))/(c*mu) - P/mu ;
    end
elseif type==1
    if k==0
        alpha = 0;
    elseif k<c_1
        alpha = R ;%- P_1/mu_1;
    else
        alpha = R - (P_w*(k-c_1+1))/(c_1*mu_1);   % - P_1/mu_1 ;
    end
else %type2 user
    if k==0
        alpha = 0;
    elseif k<c_2
        alpha = R ;
    else
        alpha = R - (P_w*(k-c_2+1))/(c_2*mu_tilda);
    end 
end
end