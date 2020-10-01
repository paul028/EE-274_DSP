%% EE 274 Digital Signal Processing 1 Lab Activity 1
% Name: Paul Vincent S. Nonat
%
%% F. Audio File Formats
%The following exercise will demonstrate the effects of using quantization 
%and sampling on audio signals. 
% 
% # Load music1.flac provided in UVLe folder. Can also be downloaded here
% # Using the MATLAB functions you have created in parts A-E, quantize, up/
%downsample using the following configurations:

%%
[y,fs] = audioread ('Sample_BeeMoved_96kHz24bit.flac')
info =audioinfo('Sample_BeeMoved_96kHz24bit.flac')
t = 0:seconds(1/fs):seconds(info.Duration);
t = t(1:end-1);

plot(t,y)
title('Original Audio')
xlabel('Time')
ylabel('Audio Signal')


%%
Q=16;
q=2/(Q+1);
target_sampling =48000;
x1=downsample(y,fs/target_sampling);
xq= q*(floor(x1/q));
plot(xq)
soundsc(xq,target_sampling)

%% X1 R=10, B=16, downsampled by 2
R=10
B=16
target_sampling =48000
x1=downsample(y,fs/target_sampling)
%x1_1d = x1(:,[1])
y1= adc_NU(x1,R,B)

%% X2 R=10, B=8, downsampled by 2
R=10
B=8
target_sampling =48000
x2=downsample(y,fs/target_sampling)


%% X3 R=10, B=4, downsampled by 2
R=10
B=4
target_sampling =48000
x3=downsample(y,fs/target_sampling)

%% X4 R=10, B=16, downsampled by 6
R=10
B=16
target_sampling =16000
x4=downsample(y,fs/target_sampling)

%% X5 R=10, B=8, downsampled by 6
R=10
B=8
target_sampling =16000
x5=downsample(y,fs/target_sampling)

%% X6 R=10, B=4, downsampled by 6
R=10
B=4
target_sampling =16000
x6=downsample(y,fs/target_sampling)

%% X7 R=10, B=16, downsampled by 12
R=10
B=16
target_sampling =8000
x7=downsample(y,fs/target_sampling)

%% X8 R=10, B=8, downsampled by 12
R=10
B=8
target_sampling =8000
x8=downsample(y,fs/target_sampling)

%% X9 R=10, B=4, downsampled by 12
R=10
B=4
target_sampling =8000
x9=downsample(y,fs/target_sampling)


%%
function y = adc_NU(x, R, B)
level = [0:R/(2^B):R-R/(2^B)];
temp = [-Inf,(level(2:end)-R/(2^(B+1))),Inf];
y = zeros(1,length(x));
i=1
y=(x >= temp(i)).*(x < temp(i+1)).*level(i)
for i = 2:length(level)
    y = y + (x >= temp(i)).*(x < temp(i+1)).*level(i);
end
end