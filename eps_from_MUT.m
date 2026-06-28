function epsMUT = eps_from_MUT(A_cal,B_cal,D_cal,S11_MUT,N)

epsMUT = zeros(N,1);
for k = 1:N
    L = log(S11_MUT(k));
    epsMUT(k) = (A_cal(k)*L + B_cal(k)) / (L + D_cal(k));
end


end