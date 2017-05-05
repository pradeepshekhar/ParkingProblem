function delta = finddelta(nhat)
global user1;
p_k_n = zeros(1,nhat+1);
d_k_n = zeros(1,nhat+1);
d_k_total_1 = 0;
for i=1:nhat+1
    d_k_n(i) = findd_k(user1,i-1);
    d_k_total_1 = d_k_total_1 + d_k_n(i);
end
for i=1:nhat+1
    p_k_n(i) = d_k_n(i)/d_k_total_1;
end
delta = p_k_n(nhat+1)
end
