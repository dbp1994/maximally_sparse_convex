clc;
clear;
close all

load('Result_data_100_trials_plot.mat')

%% L2 - Error
figure(1)
% plot(lambda_val,avg_L2Err_IRLSlp,'-*')

plot(lambda_val,avg_L2Err_IMSC_tan,'-*')
hold on
plot(lambda_val,avg_L2Err_IMSC_log,'-*')


title('L2 Error comparison');
xlabel('\lambda');
ylabel('L2 Error');
grid on;
legend('IMSC-atan','IMSC-log')
axis tight

%% L1 - Error
figure(2)
plot(lambda_val,avg_L1Err_IMSC_tan,'-*')
hold on

plot(lambda_val,avg_L1Err_IMSC_log,'-*')


title('L1 Error comparison');
xlabel('\lambda');
ylabel('L1 Error');
grid on;
legend('IMSC-atan','IMSC-log')
axis tight

%% SE - Error
figure(3)
plot(lambda_val,avg_SE_IMSC_tan,'-*')
hold on
plot(lambda_val,avg_SE_IMSC_log,'-*')


title('Support Error comparison');
xlabel('\lambda');
ylabel('Support Error');
grid on;
legend('IMSC-atan','IMSC-log')
axis tight