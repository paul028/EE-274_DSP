%% EE 274 Digital Signal Processing 1 Lab Activity 1
% Name: Paul Vincent S. Nonat
%
%% E. Quantization
%Sampling is done by periodically obtaining samples from a continuous time
%signal. The period also known as the sampling period is the reciprocal of
%sampling frequency $F_s$. Using up-sampling and down-sampling, information
%can be added or removed from a discrete time signal.
%
% # Load *signal1.wave* file in your workspace
% # Using *[y,fs] = audioread()*, import the audio and sampling rate
% information in your workspace.
%
[y,fs] = audioread('signal1.wav');
%info=audioinfo('signal1.wav');
t = 0:seconds(1/fs):seconds(info.Duration);
t = t(1:end-1);

plot(t,y)
title('Original Audio')
xlabel('Time')
ylabel('Audio Signal')
sound(y,fs)
%%
%y up-sampled by 2
M=2
y1= upsample(y,M)
t1 = 0:seconds(1/(M*fs)):seconds(info.Duration);
t1 = t1(1:end-1);
plot(t1,y1)
title('y Up-sampled by 2')
xlabel('Time')
ylabel('y1=y(n/2)')
sound(y1,fs)

%%
%y down-sampled by 2
M=2
y2= downsample(y,M)
t2 = 0:seconds(1/(fs/M)):seconds(info.Duration);
t2 = t2(1:end-1);
plot(t2,y2)
title('y down-sampled by 2')
xlabel('Time')
ylabel('y2=y(2n)')
sound(y2,fs)

%%
%y2 up-sampled by 2
M=2
y3 = upsample(y2,M)
t3 = 0:seconds(1/fs):seconds(info.Duration);
t3 = t3(1:end-1);
plot(t3,y3)
title('y2 up-sampled by 2')
xlabel('Time')
ylabel('y3=y2(n/2)')
sound(y3,fs)

%%
%y3, up-sampled by 2
M=2
y4 = upsample(y3,M)
t4 = 0:seconds(1/(fs*M)):seconds(info.Duration);
t4 = t4(1:end-1);
plot(t4,y4)
title('y3 up-sampled by 2')
xlabel('Time')
ylabel('y4=y3(n/2)')
sound(y4,fs)

%yes however some information will be loss as the upsampler cannot predict 
%the value inbetween the downsampled signales. In effect, there are loss 
%information in the upsampled signals $y_3$ and $y_4$.
