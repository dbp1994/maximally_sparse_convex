clc;
clear
close all

a = [1 -1.047 .81];
b = [1 .8];
N = 1000;
K = 50;

addpath('L1opt');
addpath(genpath('sparsify_0_5'))
addpath('ISD_v1.1')
addpath('SeDuMi_1_3')
options_AIHT = {'stopTol',1e-9,'step_size',0,'maxIter',100,'Acc',0};

intrvl_t = [5,35];
intrvl_x = [-1,1];
% sigma_w = [.1,.12,.15,.18,.2,.22,.25,.27,.3];
sigma_w = [.2];

trials = 2;

%L1 Optimization
L2Err_BPDN = zeros(trials,1); L1Err_BPDN = zeros(trials,1);SE_BPDN = zeros(trials,1);
bpdn_FN = zeros(trials,1);bpdn_FP = zeros(trials,1);
avg_L2Err_BPDN = zeros(length(sigma_w),1);
avg_L1Err_BPDN = zeros(length(sigma_w),1);
avg_SE_BPDN = zeros(length(sigma_w),1);
avg_bpdn_FP = zeros(length(sigma_w),1);
avg_bpdn_FN = zeros(length(sigma_w),1);

% L0 Optimization
L2Err_OMP = zeros(trials,1); L1Err_OMP = zeros(trials,1);SE_OMP = zeros(trials,1);
omp_FN = zeros(trials,1);omp_FP = zeros(trials,1);
avg_L2Err_OMP = zeros(length(sigma_w),1);
avg_L1Err_OMP = zeros(length(sigma_w),1);
avg_SE_OMP = zeros(length(sigma_w),1);
avg_omp_FN = zeros(length(sigma_w),1);
avg_omp_FP = zeros(length(sigma_w),1);

%IRLS1
L2Err_IRLSl1 = zeros(trials,1); L1Err_IRLSl1 = zeros(trials,1);SE_IRLSl1 = zeros(trials,1);
irlsl1_FP = zeros(trials,1);irlsl1_FN = zeros(trials,1);
avg_L2Err_IRLSl1 = zeros(length(sigma_w),1);
avg_L1Err_IRLSl1 = zeros(length(sigma_w),1);
avg_SE_IRLSl1 = zeros(length(sigma_w),1);
avg_irlsl1_FP = zeros(length(sigma_w),1);
avg_irlsl1_FN = zeros(length(sigma_w),1);

% AIHT 
L2Err_AIHT = zeros(trials,1); L1Err_AIHT = zeros(trials,1);SE_AIHT = zeros(trials,1);
aiht_FP = zeros(trials,1);aiht_FN = zeros(trials,1);
avg_L2Err_AIHT = zeros(length(sigma_w),1);
avg_L1Err_AIHT = zeros(length(sigma_w),1);
avg_SE_AIHT = zeros(length(sigma_w),1);
avg_aiht_FP = zeros(length(sigma_w),1);
avg_aiht_FN = zeros(length(sigma_w),1);

% ISD
L2Err_ISD = zeros(trials,1); L1Err_ISD = zeros(trials,1);SE_ISD = zeros(trials,1);
isd_FP = zeros(trials,1);isd_FN = zeros(trials,1);
avg_L2Err_ISD = zeros(length(sigma_w),1);
avg_L1Err_ISD = zeros(length(sigma_w),1);
avg_SE_ISD = zeros(length(sigma_w),1);
avg_SE_ISD = zeros(length(sigma_w),1);

% Gradient Descent
L2Err_GD = zeros(trials,1); L1Err_GD = zeros(trials,1);SE_GD = zeros(trials,1);
gd_FP = zeros(trials,1);gd_FN = zeros(trials,1);
avg_L2Err_GD = zeros(length(sigma_w),1);
avg_L1Err_GD = zeros(length(sigma_w),1);
avg_SE_GD = zeros(length(sigma_w),1);
avg_gd_FP = zeros(length(sigma_w),1);
avg_gd_FN = zeros(length(sigma_w),1);



% Lp
p = .7;
lambda_lp = .5;
L2Err_IRLSlp = zeros(trials,1); L1Err_IRLSlp = zeros(trials,1);
SE_IRLSlp = zeros(trials,1); irlslp_FP = zeros(trials,1);irlslp_FN = zeros(trials,1);
avg_L2Err_IRLSlp = zeros(length(sigma_w),1);
avg_L1Err_IRLSlp = zeros(length(sigma_w),1);
avg_SE_IRLSlp = zeros(length(sigma_w),1);
avg_irlslp_FP = zeros(length(sigma_w),1);
avg_irlslp_FN = zeros(length(sigma_w),1);


% IMSC_tan
L2Err_IMSC_tan = zeros(trials,1); L1Err_IMSC_tan = zeros(trials,1);
SE_IMSC_tan = zeros(trials,1); imsc_tan_FP = zeros(trials,1);imsc_tan_FN = zeros(trials,1);
avg_L2Err_IMSC_tan = zeros(length(sigma_w),1);
avg_L1Err_IMSC_tan = zeros(length(sigma_w),1);
avg_SE_IMSC_tan = zeros(length(sigma_w),1);
avg_imsc_tan_FP = zeros(length(sigma_w),1);
avg_imsc_tan_FN = zeros(length(sigma_w),1);

%IMSC_log
L2Err_IMSC_log = zeros(trials,1); L1Err_IMSC_log = zeros(trials,1);
SE_IMSC_log = zeros(trials,1); imsc_log_FP = zeros(trials,1);imsc_log_FN = zeros(trials,1);
avg_L2Err_IMSC_log = zeros(length(sigma_w),1);
avg_L1Err_IMSC_log = zeros(length(sigma_w),1);
avg_SE_IMSC_log = zeros(length(sigma_w),1);
avg_imsc_log_FP = zeros(length(sigma_w),1);
avg_imsc_log_FN = zeros(length(sigma_w),1);

lambda_val = zeros(length(sigma_w),1);
tic
for l = 1:length(sigma_w)
    for t = 1:trials
        % Generate a sparse signal x and time interval t
        [x,K,supp] = random_sig(intrvl_t,intrvl_x,K,N);
        % Generate banded Toeplitz Matrix A and B
        [A,B] = gen_toeplitz(a,b,N);
        % Number of non-zeros for AIHT
        M = K;
        % Generate Random Noise
        w = normrnd(0,sigma_w(l),[N,1]);
        % Option of ISD
        opts.xs = x;
        opts.sigma = sigma_w(l);
        opts.tol=norm(w);

        H = A\B;
%         H = normc(H);
        HtH = H'*H;
        [~,eign_val] = eig(HtH);
        eigen_min = min(diag(eign_val));
        R = eigen_min*eye(N);
        lambda_n = 3*sigma_w(l)*norm(H(:,1));
        a_n =  eigen_min/lambda_n;
        y = H*x+w;
        
%         %% L1 Optimization L1 - Magic Toolbox
%         [L2Err_BPDN(t),L1Err_BPDN(t),SE_BPDN(t),bpdn_x,bpdn_supp,bpdn_FP(t),bpdn_FN(t)] = L1_opt(x,H,y,norm(w));
%         %% L0 Optimization - OMP
%         [L2Err_OMP(t),L1Err_OMP(t),SE_OMP(t),omp_x,omp_supp,omp_FP(t),omp_FN(t)] = OMP_algo(H,x,y,norm(w));
%         %% Accelerated Iterative Hard Threshold
%         [L2Err_AIHT(t),L1Err_AIHT(t),SE_AIHT(t),aiht_x,aiht_supp,aiht_FP(t),aiht_FN(t)] = AIHT(y,x,H,N,M,norm(w),options_AIHT);
%         %% Iterative Support Detection
% %       [L2Err_ISD(t),L1Err_ISD(t),SE_ISD(t),isd_x,isd_supp,isd_FP(t),isd_FN(t)] = Threshold_ISD_1D(A,y,opts);
%         %% Lp Optimization
%          [L2Err_IRLSlp(t),L1Err_IRLSlp(t),SE_IRLSlp(t),irlslp_x,irlslp_supp,irlslp_FP(t),irlslp_FN(t)] = IRLS_lp(H,y,x,p);
%         
%         %% Gradient Descent
%          [L2Err_GD(t),L1Err_GD(t),SE_GD(t),gd_x,gd_supp,gd_FP(t),gd_FN(t)] = grad_descent(H,y,x,lambda_n,a_n);
        %% IMSC Log- IRLS1 initialization
        %     % Step 1 : Initialization Find L1 the norm solution
%         [L2Err_IRLSl1(t),L1Err_IRLSl1(t),SE_IRLSl1(t),irlsl1_x,irlsl1_supp,irlsl1_FP(t),irlsl1_FN(t)] = IRLS1(H,y,x,lambda_n,sigma_w(l));
%             Ks = zeros(100,length(x));
%             Ks(1,:) = 1;
%             x_i = irlsl1_x;
%             for i  = 2:100
%                 Ks(i,find(x_i~=0)) = 1;
%                 supp_Ks = find(Ks(i,:));
%                 if length(find(Ks(i-1,:))) <= length(find(Ks(i,:)))
%                     break;
%                 end
%                     H_i = H(:,supp_Ks);
%                     HitHi = H_i'*H_i;
%                     [~,eign_val] = eig(HitHi);
%                     eigen_min = min(diag(eign_val));
%                     R_i = SDP_rn(HitHi,eigen_min);
% %                     R_i = eigen_min*eye(length(supp_Ks));
%                     lambda_n = 3*sigma_w(l)*norm(H_i(:,1));
%                     a_n =  eigen_min/lambda_n; 
%                     if i == 2
%                         lambda = lambda_n;
%                     else
%                         lambda = lambda_n/eigen_min;
%                     end
%                   [L2Err_IMSC_log(t),L1Err_IMSC_log(t),SE_IMSC_log(t),imsc_log_x,imsc_log_supp,imsc_log_FP(t),imsc_log_FN(t)] = IMSC_log(H_i,y,x,a_n,lambda,norm(w),supp_Ks);
%                 x_i = imsc_log_x;     
%               
%             end
        % IMSC atan- IRLS1 initialization
            % Step 1 : Initialization Find L1 the norm solution
              
        [L2Err_IRLSl1(t),L1Err_IRLSl1(t),SE_IRLSl1(t),irlsl1_x,irlsl1_supp,irlsl1_FP(t),irlsl1_FN(t)] = IRLS1(H,y,x,lambda_n,sigma_w(l));
%         [L2Err_IRLSl1(t),L1Err_IRLSl1(t),SE_IRLSl1(t),irlsl1_x,irlsl1_supp,irlsl1_FP(t),irlsl1_FN(t)] = L1_opt(x,H,y,norm(w));
            Ks = zeros(100,length(x));
            Ks(1,:) = 1;
            x_i = irlsl1_x;
            for i  = 2:100
                Ks(i,find(x_i~=0)) = 1;
                supp_Ks = find(Ks(i,:));
                if length(find(Ks(i-1,:))) <= length(find(Ks(i,:)))
                    break;
                end
                    H_i = H(:,supp_Ks);
                    HitHi = H_i'*H_i;
                    [~,eign_val] = eig(HitHi);
                    eigen_min = min(diag(eign_val));
%                     R_i = SDP_rn(HitHi,eigen_min);
%                     r_n = diag(R_i);
                    R_i = eigen_min*eye(length(supp_Ks));
                    lambda_n = 3*sigma_w(l)*norm(H_i(:,1));
%                     a_n =  r_n./lambda_n; 
%                     lambda = lambda_n./r_n;
%                     a_n =  eigen_min/lambda_n; 
                    a_n = 1; 
                    if i == 2
                        lambda = lambda_n;
                    else
                        lambda = lambda_n/eigen_min;
%                         lambda = lambda_n./r_n;
                    end
   
                    [L2Err_IMSC_tan(t),L1Err_IMSC_tan(t),SE_IMSC_tan(t),imsc_tan_x,imsc_tan_supp,imsc_tan_FP(t),imsc_tan_FN(t)] = IMSC_tan(H_i,y,x,a_n,lambda,norm(w),supp_Ks);
                x_i = imsc_tan_x;
                
            end
    end
    
    avg_L2Err_IRLSl1(l) = mean(L2Err_IRLSl1);
    avg_L1Err_IRLSl1(l) = mean(L1Err_IRLSl1);
    avg_SE_IRLSl1(l) = mean(SE_IRLSl1);
     avg_irlsl1_FP(l) = mean(irlsl1_FP);
      avg_irlsl1_FN(l) = mean(irlsl1_FN);
    
    avg_L2Err_IMSC_tan(l) = mean(L2Err_IMSC_tan);
    avg_L1Err_IMSC_tan(l) = mean(L1Err_IMSC_tan);
    avg_SE_IMSC_tan(l) = mean(SE_IMSC_tan);
     avg_imsc_tan_FP(l) = mean(imsc_tan_FP);
      avg_imsc_tan_FN(l) = mean(imsc_tan_FN);
    
    avg_L2Err_IMSC_log(l) = mean(L2Err_IMSC_log);
    avg_L1Err_IMSC_log(l) = mean(L1Err_IMSC_log);
    avg_SE_IMSC_log(l) = mean(SE_IMSC_log);
    avg_imsc_log_FP(l) = mean(imsc_log_FP);
    avg_imsc_log_FN(l) = mean(imsc_log_FN);

    avg_L2Err_IRLSlp(l) = mean(L2Err_IRLSlp);
    avg_L1Err_IRLSlp(l) = mean(L1Err_IRLSlp);
    avg_SE_IRLSlp(l) = mean(SE_IRLSlp);
    avg_irlslp_FP(l) = mean(irlslp_FP);
    avg_irlslp_FN(l) = mean(irlslp_FN);
    
    avg_L2Err_BPDN(l) = mean(L2Err_BPDN);
    avg_L1Err_BPDN(l) = mean(L1Err_BPDN);
    avg_SE_BPDN(l) = mean(SE_BPDN);
     avg_bpdn_FP(l) = mean(bpdn_FP);
      avg_bpdn_FN(l) = mean(bpdn_FN);
    
    avg_L2Err_OMP(l) = mean(L2Err_OMP);
    avg_L1Err_OMP(l) = mean(L1Err_OMP);
    avg_SE_OMP(l) = mean(SE_OMP);
    avg_omp_FP(l) = mean(omp_FP);
    avg_omp_FN(l) = mean(omp_FN);
    
    avg_L2Err_AIHT(l) = mean(L2Err_AIHT);
    avg_L1Err_AIHT(l) = mean(L1Err_AIHT);
    avg_SE_AIHT(l) = mean(SE_AIHT);
    avg_aiht_FP(l) = mean(aiht_FP);
    avg_aiht_FN(l) = mean(aiht_FN);
    
%     avg_L2Err_ISD(l) = mean(L2Err_ISD);
%     avg_L1Err_ISD(l) = mean(L1Err_ISD);
%     avg_SE_ISD(l) = mean(SE_ISD);
    
    avg_L2Err_GD(l) = mean(L2Err_GD);
    avg_L1Err_GD(l) = mean(L1Err_GD);
    avg_SE_GD(l) = mean(SE_GD);
    avg_gd_FP(l) = mean(gd_FP);
    avg_gd_FN(l) = mean(gd_FN);
     
    lambda_val(l) = lambda_n;
end

L2Err_All = [avg_L2Err_AIHT;avg_L2Err_BPDN;avg_L2Err_OMP;avg_L2Err_GD;avg_L2Err_IMSC_log;avg_L2Err_IMSC_tan;avg_L2Err_IRLSl1;avg_L2Err_IRLSlp];
L1Err_All = [avg_L1Err_AIHT;avg_L1Err_BPDN;avg_L1Err_OMP;avg_L1Err_GD;avg_L1Err_IMSC_log;avg_L1Err_IMSC_tan;avg_L1Err_IRLSl1;avg_L1Err_IRLSlp];
SE_All = [avg_SE_AIHT;avg_SE_BPDN;avg_SE_GD;avg_SE_OMP;avg_SE_IMSC_log;avg_SE_IMSC_tan;avg_SE_IRLSl1;avg_SE_IRLSlp];
FP_All = [avg_aiht_FP;avg_bpdn_FP;avg_gd_FP;avg_omp_FP;avg_imsc_log_FP;avg_imsc_tan_FP;avg_irlsl1_FP;avg_irlslp_FP];
FN_All = [avg_aiht_FN;avg_bpdn_FN;avg_gd_FN;avg_omp_FN;avg_imsc_log_FN;avg_imsc_tan_FN;avg_irlsl1_FN;avg_irlslp_FN];


toc
aa = 1;
% plot(lambda_val,avg_L2Err_IRLSlp)
% hold on
% plot(lambda_val,avg_L2Err_BPDN)
% hold on
% plot(lambda_val,avg_L2Err_OMP)
% hold on
% plot(lambda_val,avg_L2Err_AIHT)
% hold on
% % plot(lambda_val,avg_L2Err_ISD)
% % hold on
% plot(lambda_val,avg_L2Err_GD)
% legend('MSC','L1','OMP','AIHT','ISD','GD');
% grid on
% title('L2-Error Comparison');
% xlabel('lambda');
% ylabel('Normalized L2- Error');

