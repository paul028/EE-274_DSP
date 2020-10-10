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
% Load Input Signals Input Signals
[x1,fs1] = audioread('inputs/x1.wav');
%divide x1 signal into two channel for computational convenience
x1_c1 = x1(:,1); 
x1_c2 = x1(:,2);
[x2,fs2] = audioread('inputs/x2.wav');

[x3,fs3] = audioread('inputs/x3.wav');

[x4,fs4] = audioread('inputs/x4.wav');

[x5,fs5] = audioread('inputs/x5.wav');

%%
%
%Simulate System 1
% <<system1.PNG>>
y1_c1=dt_1(x1_c1)
y1_c2=dt_1(x1_c2)
y1_2=dt_1(x2)
y1_3=dt_1(x3)
y1_4=dt_1(x4)
y1_5=dt_1(x5)

%%
%Simulate System 2
% <<system2.PNG>>
y2_c1=dt_2(x1_c1)
y2_c2=dt_2(x1_c2)
y2_2=dt_2(x2)
y2_3=dt_2(x3)
y2_4=dt_2(x4)
y2_5=dt_2(x5)

%%
%Simulate System 3
y3_c1=dt_3(x1_c1)
y3_c2=dt_3(x1_c2)
y3_2=dt_3(x2)
y3_3=dt_3(x3)
y3_4=dt_3(x4)
y3_5=dt_3(x5)

%%
%Simulate System 4
L1=50
L2=100
L3=400

y4_c1_L1=dt_4(x1_c1,L1)
y4_c1_L2=dt_4(x1_c1,L2)
y4_c1_L3=dt_4(x1_c1,L3)

y4_c2_L1=dt_4(x1_c1,L1)
y4_c2_L2=dt_4(x1_c1,L2)
y4_c2_L3=dt_4(x1_c1,L3)

y4_2_L1=dt_4(x2,L1)
y4_2_L2=dt_4(x2,L2)
y4_2_L3=dt_4(x2,L3)

y4_3_L1=dt_4(x3,L1)
y4_3_L2=dt_4(x3,L2)
y4_3_L3=dt_4(x3,L3)

y4_4_L1=dt_4(x4,L1)
y4_4_L2=dt_4(x4,L2)
y4_4_L3=dt_4(x4,L3)

y4_5_L1=dt_4(x5,L1)
y4_5_L2=dt_4(x5,L2)
y4_5_L3=dt_4(x5,L3)