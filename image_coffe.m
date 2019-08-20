clc;
clear;
close all;

var = .02;
lambda_n = 3*sqrt(var);

I = imread('cameraman.tif');
[IA1,IH1,IV1,ID1]= dwt2(I,'haar');
[IA2,IH2,IV2,ID2]= dwt2(IA1,'haar');
[IA3,IH3,IV3,ID3]= dwt2(IA2,'haar');

e = sqrt(length(I(:)))*sqrt(var);

% Vectorizatin of matrices.
x_A3 = IA3(:);
x_H3 = IH3(:);
x_V3 = IV3(:);
x_D3 = ID3(:);

x_H2 = IH2(:);
x_V2 = IV2(:);
x_D2 = ID2(:);
imshow(I);
I_n = imnoise(I,'gaussian',0,var);
imshow(I_n);

[A1,H1,V1,D1]= dwt2(I_n,'haar');
[A2,H2,V2,D2]= dwt2(A1,'haar');
[A3,H3,V3,D3]= dwt2(A2,'haar');
col3 = size(H3,2);
col2 = size(H2,2);
col1 = size(H1,2);


% Vectorizatin of matrices.
y_H3 = H3(:);
y_V3 = V3(:);
y_D3 = D3(:);

y_H2 = H2(:);
y_V2 = V2(:);
y_D2 = D2(:);

y = {y_H3,y_V3,y_D3,y_H2,y_V2,y_D2};
x = {x_H3,x_V3,x_D3,x_H2,x_V2,x_D2};


for j = 1:length(y)
    
    H = eye(length(y{j}));
    % IMSC - IRLS1 initialization
    % Step 1 : Initialization Find L1 the norm solution
    [L2Err_IRLSl1,L1Err_IRLSl1,SE_IRLSlp,irlsl1_x{j},irlsl1_supp{j},irlsl1_FP,irlsl1_FN] = IRLS1(H,y{1,j},x{1,j},lambda_n,e);
    Ks = zeros(100,length(y{j}));
    Ks(1,:) = 1;
    x_i = irlsl1_x{1,j};
    irlslp_x{j} = x_i;
    for i  = 2:100
        Ks(i,find(x_i~=0)) = 1;
        supp_Ks = find(Ks(i,:));
        if length(find(Ks(i-1,:))) <= length(find(Ks(i,:)))
            break;
        end
        H_i = H(:,supp_Ks);
        [~,eign_val] = eig(H_i'*H_i);
        eigen_min = min(diag(eign_val));
        
        R_i = eigen_min*eye(length(supp_Ks));
%         lambda_n = 3*sigma_w(l)*norm(H_i(:,1));
        a_n =  eigen_min/lambda_n;
        
        [L2Err_IRLSlp,L1Err_IRLSlp,SE_IRLSlp,irlslp_x{j},irlslp_supp{j},irlslp_FP,irlslp_FN] = IMSC_tan(H_i,y{1,j},x{1,j},a_n,lambda_n,e/100,supp_Ks);
        
        x_i = irlslp_x{1,j};
    end
end

aa = 1;

[y_A2]= idwt2(A3,vec2mat(irlslp_x{1,1},col3),vec2mat(irlslp_x{1,2},col3),vec2mat(irlslp_x{1,3},col3),'haar');
[y_A1]= idwt2(y_A2,vec2mat(irlslp_x{1,4},col2),vec2mat(irlslp_x{1,5},col2),vec2mat(irlslp_x{1,6},col2),'haar');
[Y_J]= idwt2(y_A1,H1,V1,D1,'haar');