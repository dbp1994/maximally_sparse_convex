function [L2Err,L1Err,SE,bpdn_x,bpdn_supp,avg_FP,avg_FN] = L1_opt(x,A,b_noisy,eps)

x0 = A'*b_noisy;
% x_idx = zeros(size(x,1),1);
bpdn_x = l1qc_logbarrier(x0, A,[], b_noisy,eps);
bpdn_x(find(abs(bpdn_x)<.05)) = 0;
bpdn_supp = find(bpdn_x~=0);

[L2Err,L1Err,SE,avg_FP,avg_FN] = matric_cal(bpdn_x,x);



end

