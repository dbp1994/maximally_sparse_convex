clc;
clear;
close all

load('Result_data_50_trials.mat')

%% Original Signal and Corrupted Signal
figure(1)
subplot(2,1,1)
plot(x)
title('Original Signal');
xlabel('sample time');
ylabel('Amplitude');
grid on;
axis tight
hold on
subplot(2,1,2)
plot(y)
title('Corrupted Signal \sigma = .2');
xlabel('sample time');
ylabel('Amplitude');
grid on;
axis tight

%% Original Signal and IMSC-atan recovered signal
figure(2)
plot(x)
hold on
plot(imsc_tan_x)
title('Original signal and IMSC-atan recovered signal');
xlabel('sample time');
ylabel('Amplitude');
legend('Actual signal','IMSC-atan')
grid on;
axis tight


%% Original Signal and IMSC-log recovered signal
figure(3)
plot(x)
hold on
plot(imsc_log_x)
title('Original signal and IMSC-log recovered signal');
xlabel('sample time');
ylabel('Amplitude');
legend('Actual signal','IMSC-log')
grid on;
axis tight

%% Original Signal and IRL1 recovered signal
figure(4)
plot(x)
hold on
plot(irlsl1_x)
title('Original signal and IRL1 recovered signal');
xlabel('sample time');
ylabel('Amplitude');
legend('Actual signal','IRL1')
grid on;
axis tight

%% Original Signal and IRLp recovered signal
figure(5)
plot(x)
hold on
plot(irlslp_x)
title('Original signal and IRLp recovered signal');
xlabel('sample time');
ylabel('Amplitude');
legend('Actual signal','IRLp')
grid on;
axis tight

%% Original Signal and IRLp recovered signal
figure(6)
plot(x)
hold on
plot(aiht_x)
title('Original signal and AIHT recovered signal');
xlabel('sample time');
ylabel('Amplitude');
legend('Actual signal','AIHT')
grid on;
axis tight