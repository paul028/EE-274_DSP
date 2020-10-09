%% Paul Vincent S. Nonat 2018-21366
% EE 274 Digital Signal Processing 1 Lab Activity 2
%

%% A. DIFFERENCE EQUATIONS IN MATLAB
%Discrete time systems can be represented using block diagrams and
%difference equation. Shown below is an example of a discrete time ARMA
%system:
%%
% 
% $$y[n] = x[n] + 5x[n-1] + 2y[n-1]$$
% 
%%
% 
% <<ARMA_system.PNG>>
% 

%% B. SYSTEM RESPONSE CALCULATION
% # Using recursion
% We can simulate the recursion using the for loop function. FOr the code
% below, the simulation runs on the entire domain of n:
% $y[n] = x[n] + 5x[n-1] +2y[n-1]$
x = randn(1,5); %random input signal x
y = zeros(1,length(x)); %initialize output signal y
for n = 1:length(x) %5 time indices
if n<2
y(n) = x(n);
else
y(n) = x(n) + 5*x(n-1) + 2*y(n-1);
end
end
figure;
subplot 211
stem(1:5,x); title('input signal');
subplot 212
stem(1:5,y); title('output signal');

%%
% # Using the impulse response and DT convolution
% We can also simulate response using the *impz ()* function and * conv ();*
b = [1 5]; a = [1 -2];
h = impz(b,a); %impulse response
y = conv(h,x);
figure;
subplot 211
stem(1:5,x(1:5)); title('input signal');
subplot 212
stem(1:5,y(1:5)); title('output signal');

%% C. EXERCISES
% Given the following discrete time systems:
% # System 1
% <<system1.PNG>>
% # System 2
% <<system2.PNG>>
% # System 3
% <<system3.PNG>>
% # System 4
% <<system4.PNG>>
%
% For each system, determine the following using MATLAB:

%Create a function file (M-file) that returns a vector y from a given input x. 
%(Make sure the lengths are the same and are using the same sampling period). 
%For systems 1&2, perform the recursive method. For systems 3&4, use the 
%impulse response method. System 4, should have an extra input L. Assume 
%zero initial conditions.
%
% *Format:* y = dt_1(x), y = dt_2(x), y = dt_3(x), y = dt_4(x,L)
%
% Investigate the output response from the given input signals. Use
% *figure()* and *subplot()* to compare the domain plot of x versus y. For
% system 4, try L = 100.
%
%Compare the input and output signals by listening using *soundsc(y,fs)*.
%For system 4, try the following values of L = 50, 100, 400.
%
% Answer the following questions (Use MATLAB to support your answer):
% # Is the system BIBO stable? (Hint: use impz ())
% # Is the system causal?
% # Is the system FIR or IIR?
% What does the system do? (Based from #2&#3, you may also use freqz (b,a)
%to have an insight on the filter / frequency response)

%%
% Input Signals
[x1,fs1] = audioread('inputs/x1.wav');
info=audioinfo('inputs/x1.wav');
t1 = 0:seconds(1/fs1):seconds(info.Duration);
t1 = t1(1:end-1);

[x2,fs2] = audioread('inputs/x2.wav');
info=audioinfo('inputs/x2.wav');
t2 = 0:seconds(1/fs2):seconds(info.Duration);
t2 = t2(1:end-1);

[x3,fs3] = audioread('inputs/x3.wav');
info=audioinfo('inputs/x1.wav');
t3 = 0:seconds(1/fs3):seconds(info.Duration);
t3 = t3(1:end-1);

[x4,fs4] = audioread('inputs/x4.wav');
info=audioinfo('inputs/x4.wav');
t4 = 0:seconds(1/fs4):seconds(info.Duration);
t4 = t4(1:end-1);

[x5,fs5] = audioread('inputs/x5.wav');
info=audioinfo('inputs/x5.wav');
t5 = 0:seconds(1/fs5):seconds(info.Duration);
t5 = t5(1:end-1);

%%
% System 1
% <<system1.PNG>>
y1_1=dt_1(x1)
figure;
subplot 211
stem(x1)
title('input signal (x1)')
subplot 212
stem(y1_1)
title('output signal (y1_1)')

y1_2=dt_1(x2)
figure;
subplot 211
stem(x2)
title('input signal (x2)')
subplot 212
stem(y1_2)
title('output signal (y1_2)')

y1_3=dt_1(x3)
figure;
subplot 211
stem(x3)
title('input signal (x3)')
subplot 212
stem(y1_3)
title('output signal (y1_3)')

y1_4=dt_1(x4)
figure;
subplot 211
stem(x4)
title('input signal (x4)')
subplot 212
stem(y1_4)
title('output signal (y1_4)')

y1_5=dt_1(x5)
figure;
subplot 211
stem(x5)
title('input signal (x5)')
subplot 212
stem(y1_5)
title('output signal (y1_5)')

% Answer to # 4

%%
% System 2
% <<system2.PNG>>
y2_1=dt_2(x1)
figure;
subplot 211
stem(x1)
title('input signal (x1)')
subplot 212
stem(y2_1)
title('output signal (y2_1)')

y2_2=dt_2(x2)
figure;
subplot 211
stem(x2)
title('input signal (x2)')
subplot 212
stem(y2_2)
title('output signal (y2_2)')

y2_3=dt_2(x3)
figure;
subplot 211
stem(x3)
title('input signal (x3)')
subplot 212
stem(y2_3)
title('output signal (y2_3)')

y2_4=dt_2(x4)
figure;
subplot 211
stem(x4)
title('input signal (x4)')
subplot 212
stem(y2_4)
title('output signal (y2_4)')

y2_5=dt_2(x5)
figure;
subplot 211
stem(x5)
title('input signal (x5)')
subplot 212
stem(y2_5)
title('output signal (y2_5)')

% Answer to # 4

%%
% System 3
% <<system3.PNG>>
y3_1=dt_3(x1)
figure;
subplot 211
stem(x1)
title('input signal (x1)')
subplot 212
stem(y3_1)
title('output signal (y3_1)')

y3_2=dt_3(x2)
figure;
subplot 211
stem(x2)
title('input signal (x2)')
subplot 212
stem(y3_2)
title('output signal (y3_2)')

y3_3=dt_3(x3)
figure;
subplot 211
stem(x3)
title('input signal (x3)')
subplot 212
stem(y3_3)
title('output signal (y3_3)')

y3_4=dt_3(x4)
figure;
subplot 211
stem(x4)
title('input signal (x4)')
subplot 212
stem(y3_4)
title('output signal (y3_4)')

y3_5=dt_3(x5)
figure;
subplot 211
stem(x5)
title('input signal (x5)')
subplot 212
stem(y3_5)
title('output signal (y3_5)')

% Answer to # 4

%%
% System 4
% <<system4.PNG>>
L=100
y4_1=dt_4(x1,L)
figure;
subplot 211
stem(x1)
title('input signal (x1)')
subplot 212
stem(y4_1)
title('output signal (y4_1)')

y4_2=dt_4(x2,L)
figure;
subplot 211
stem(x2)
title('input signal (x2)')
subplot 212
stem(y4_2)
title('output signal (y4_2)')

y4_3=dt_3(x3,L)
figure;
subplot 211
stem(x3)
title('input signal (x3)')
subplot 212
stem(y4_3)
title('output signal (y4_3)')

y4_4=dt_4(x4,L)
figure;
subplot 211
stem(x4)
title('input signal (x4)')
subplot 212
stem(y4_4)
title('output signal (y4_4)')

y4_5=dt_4(x5)
figure;
subplot 211
stem(x5)
title('input signal (x5)')
subplot 212
stem(y4_5)
title('output signal (y4_5)')

% Answer to # 4

