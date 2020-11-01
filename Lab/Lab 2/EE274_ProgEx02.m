%% Paul Vincent S. Nonat 2018-21366
% EE 274 Digital Signal Processing 1 Lab Activity 2
%

%% A. DIFFERENCE EQUATIONS IN MATLAB
% Discrete time systems can be represented using block diagrams and
% difference equation. Shown below is an example of a discrete time ARMA
% system:
%%
% 
% $$y[n] = x[n] + 5x[n-1] + 2y[n-1]$$
% 
%%
% 
% <<ARMA_SYSTEM.PNG>>
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

%% 
% Create a function file (M-file) that returns a vector y from a given input x. 
% (Make sure the lengths are the same and are using the same sampling period). 
% For systems 1&2, perform the recursive method. For systems 3&4, use the 
% impulse response method. System 4, should have an extra input L. Assume 
% zero initial conditions.
%
%  *Format:* y = dt_1(x), y = dt_2(x), y = dt_3(x), y = dt_4(x,L)
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
% Simulate System 1
% <<system1.png>>
y1_c1=dt_1(x1_c1);
y1_c2=dt_1(x1_c2);
%%
y1_2=dt_1(x2);
y1_3=dt_1(x3);
y1_4=dt_1(x4);
y1_5=dt_1(x5);
%%
% combine y1_1 channels
y1_1 = [y1_c1(:),y1_c2(:)];
%%
% Simulate System 2
% <<system2.png>>
y2_c1=dt_2(x1_c1);
y2_c2=dt_2(x1_c2);
%%
y2_2=dt_2(x2);
y2_3=dt_2(x3);
y2_4=dt_2(x4);
y2_5=dt_2(x5);

%%
%combine y2_1
y2_1 = [y2_c1(:),y2_c2(:)];

%%
% Simulate System 3
% <<system3.png>>
y3_c1=dt_3(x1_c1);
y3_c2=dt_3(x1_c2);
y3_2=dt_3(x2);
y3_3=dt_3(x3);
y3_4=dt_3(x4);
y3_5=dt_3(x5);

%%
% combine y3_1
y3_1 = [y3_c1(:),y3_c2(:)];

%%
%Simulate System 4
% <<system4.png>>
L1=50;
L2=100;
L3=400;

y4_c1_L1=dt_4(x1_c1,L1);
y4_c1_L2=dt_4(x1_c1,L2);
y4_c1_L3=dt_4(x1_c1,L3);

y4_c2_L1=dt_4(x1_c1,L1);
y4_c2_L2=dt_4(x1_c1,L2);
y4_c2_L3=dt_4(x1_c1,L3);

y4_2_L1=dt_4(x2,L1);
y4_2_L2=dt_4(x2,L2);
y4_2_L3=dt_4(x2,L3);

y4_3_L1=dt_4(x3,L1);
y4_3_L2=dt_4(x3,L2);
y4_3_L3=dt_4(x3,L3);

y4_4_L1=dt_4(x4,L1);
y4_4_L2=dt_4(x4,L2);
y4_4_L3=dt_4(x4,L3);

y4_5_L1=dt_4(x5,L1);
y4_5_L2=dt_4(x5,L2);
y4_5_L3=dt_4(x5,L3);

%%
%combine y4_1_L1
%combine y4_1_L2
%combine y4_1_L3
y4_1_L1 = [y4_c1_L1(:),y4_c2_L1(:)];
y4_1_L2 = [y4_c1_L2(:),y4_c2_L2(:)];
y4_1_L3 = [y4_c1_L3(:),y4_c2_L3(:)];



%% 
% 2.Investigate the output response from the given input signals. Use
% *figure()* and *subplot()* to compare the domain plot of x versus y. For
% system 4, try L = 100.
%% 
% 3. Compare the input and output signals by listening using *soundsc(y,fs)*.
% For system 4, try the following values of L = 50, 100, 400.
%%
% System 1 Input X1
figure;
subplot 211
stem(1:length(x1),x1);
title(' input signal x1');
subplot 212
stem(1:length(y1_1),y1_1);
title('Output signal y1_1');
%%
% System 1 X1
soundsc(x1,fs1)
soundsc(y1_1,fs1)
%%
% Observation:
% Base from the graph of System 1 with input X1, there is a small 
% difference between the input and output of audio x1 although the
% difference cannot be distinguish through hearing.


%%
% System 1 Input X2
figure;
subplot 211
stem(1:length(x2),x2);
title(' input signal x2')
subplot 212
stem(1:length(y1_2),y1_2);
title('Output signal y1_2')
%%
% System 1 X2
soundsc(x2,fs2)
soundsc(y1_2,fs2)
%%
% Observation:
% Base from the graph of System 1 with input X2, there is a
% small difference between the input and output of speech x2 although the
% difference cannot be distinguish through hearing.

%%
% System 1 Input X3
figure;
subplot 211
stem(1:length(x3),x3);
title(' input signal x3')
subplot 212
stem(1:length(y1_3),y1_3);
title('Output signal y1_3')
%%
% System 1 X3
soundsc(x3,fs3)
soundsc(y1_3,fs3)
%%
% Observation:
% Base from the graph of System 1 with input X3, there is a
% small difference between the input and output of audio x3 although the
% difference cannot be distinguish through hearing.

%%
% System 1 Input X4
figure;
subplot 211
stem(1:length(x4),x4);
title(' input signal x4')
subplot 212
stem(1:length(y1_4),y1_4);
title('Output signal y1_4')
%%
% System 1 X4
soundsc(x4,fs4)
soundsc(y1_4,fs4)
%%
% Observation:
% Base from the graph of System 1 with input X4, there is a
% small difference between the input and output of audio x4 although the
% difference cannot be distinguish through hearing.

%%
% System 1 Input X5
figure;
subplot 211
stem(1:length(x5),x5);
title(' input signal x5')
subplot 212
stem(1:length(y1_5),y1_5);
title('Output signal y1_5')
%%
% System 1 X5
soundsc(x5,fs5)
soundsc(y1_5,fs5)
%%
% Observation:
% Base from the graph of System 1 with input X5, there is a
% small difference between the input and output of audio x5 although the
% difference cannot be distinguish through hearing.

%%
% System 2 Input X1
figure;
subplot 211
stem(1:length(x1),x1);
title(' input signal x1')
subplot 212
stem(1:length(y2_1),y2_1);
title('Output signal y2_1')

%%
% System 2 X1
soundsc(x1,fs1)
soundsc(y2_1,fs1)
%%
% Observation:
% Signal x1 become smoother after passing system 2


%%
% System 2 Input X2
figure;
subplot 211
stem(1:length(x2),x2);
title(' input signal x2')
subplot 212
stem(1:length(y2_2),y2_2);
title('Output signal y2_2')
%%
% System 2 X2
soundsc(x2,fs2)
soundsc(y2_2,fs2)
%%
% Observation:
% Signal x2 become smoother after passing system 2

%%
% System 2 Input X3
figure;
subplot 211
stem(1:length(x3),x3);
title(' input signal x3')
subplot 212
stem(1:length(y2_3),y2_3);
title('Output signal y2_3')
%%
% System 2 X3
soundsc(x3,fs3)
soundsc(y2_3,fs3)
%%
% Observation: Signal x3 become smoother after passing system 2
%%
% System 2 Input X4
figure;
subplot 211
stem(1:length(x4),x4);
title(' input signal x4')
subplot 212
stem(1:length(y2_4),y2_4);
title('Output signal y2_4')
%%
% System 2 X4
soundsc(x4,fs4)
soundsc(y2_4,fs4)
%%
% Observation: Signal x4 become smoother after passing system 2
%%
% System 2 Input X5
figure;
subplot 211
stem(1:length(x5),x5);
title(' input signal x5')
subplot 212
stem(1:length(y2_5),y2_5);
title('Output signal y2_5')
%%
% System 2 X5
soundsc(x5,fs5)
soundsc(y2_5,fs5)
%%
% Observation: Signal x5 become smoother after passing system 2

%%
% System 3 Input X1
figure;
subplot 211
stem(1:length(x1),x1);
title(' input signal x1')
subplot 212
stem(1:length(y3_1),y3_1);
title('Output signal y3_1')
%%
% System 3 X1
soundsc(x1,fs1)
soundsc(y3_1,fs1)
%%
% Observation: Signal x1 is smoother after passing system 3 as compared to 
%system 2

%%
% System 3 Input X2
figure;
subplot 211
stem(1:length(x2),x2);
title(' input signal x2')
subplot 212
stem(1:length(y3_2),y3_2);
title('Output signal y3_2')
%%
% System 3 X2
soundsc(x2,fs2)
soundsc(y3_2,fs2)
%% 
%Observation: Signal x2 is smoother after passing system 3 as compared 
%to system 2
%%
% System 3 Input X3
figure;
subplot 211
stem(1:length(x3),x3);
title(' input signal x3')
subplot 212
stem(1:length(y3_3),y3_3);
title('Output signal y3_3')
%%
% System 3 X3
soundsc(x3,fs3)
soundsc(y3_3,fs3)
%%
% Observation: Signal x3 is smoother after passing system 3 as compared 
% to system 2

%%
% System 3 Input X4
figure;
subplot 211
stem(1:length(x4),x4);
title(' input signal x4')
subplot 212
stem(1:length(y3_4),y3_4);
title('Output signal y3_4')
%%
% System 3 X4
soundsc(x4,fs4)
soundsc(y3_4,fs4)
%%
% Observation: Signal x4 is smoother after passing system 3 as compared 
%to system 2
%%
% System 3 Input X5
figure;
subplot 211
stem(1:length(x5),x5);
title(' input signal x5')
subplot 212
stem(1:length(y3_5),y3_5);
title('Output signal y3_5')
%%
% System 3 X5
soundsc(x5,fs5)
soundsc(y3_5,fs5)
%%
% Observation: Signal x5 is smoother after passing system 3 as compared 
%to system 2

%%
% System 4 Input X1
figure;
subplot 411
stem(1:length(x1),x1);
title(' input signal x1')
subplot 412
stem(1:length(y4_1_L1),y4_1_L1);
title('L=50 Output signal')
subplot 413
stem(1:length(y4_1_L2),y4_1_L2);
title('L=100 Output signal')
subplot 414
stem(1:length(y4_1_L3),y4_1_L3);
title('L=400 Output signal')
%%
% System 4 X1
soundsc(x1,fs1)
soundsc(y4_1_L1,fs1)
soundsc(y4_1_L2,fs1) 
soundsc(y4_1_L3,fs1) 
%%
% As compared to the original signal, the signal x1 after passing system
% 4 appeared to be downsampled. Upon listening to the output of System 4 at
% L=50, the singer in the audio seems to be out of tempo with some echos in
% appearing in the background. Also the bell instrument become louder as
% compared to the singer. At L=100, the bell becomes more apparent and the
% delay in the singer's voice was longer. Finaly, at L=400, the audio cannot
% be heard clearly as there is a loud horn in the background.


%%
% System 4 Input X2
figure;
subplot 411
stem(1:length(x2),x2);
title(' input signal x2')
subplot 412
stem(1:length(y4_2_L1),y4_2_L1);
title('L=50 Output signal')
subplot 413
stem(1:length(y4_2_L2),y4_2_L2);
title('L=100 Output signal')
subplot 414
stem(1:length(y4_2_L3),y4_2_L3);
title('L=400 Output signal')
%%
% System 4 X2
soundsc(x2,fs2)
soundsc(y4_2_L1,fs2) % distorted voice
soundsc(y4_2_L2,fs2) % robotic voice
soundsc(y4_2_L3,fs2) % voice distortion is intensified with a feedback similar to 
%%
% At L=50, the voice becomes distorted with loud noise in the background, at 
% L=100, the voice's pitch becomes lower making it sound like a male robot.
% Finally, At L=400, the pitch of the input signal is at the lowest compared
% to the previous L values.


%%
% System 4 Input X3
figure;
subplot 411
stem(1:length(x3),x3);
title(' input signal x3')
subplot 412
stem(1:length(y4_3_L1),y4_3_L1);
title('L=50 Output signal')
subplot 413
stem(1:length(y4_3_L2),y4_3_L2);
title('L=100 Output signal')
subplot 414
stem(1:length(y4_3_L3),y4_3_L3);
title('L=400 Output signal')
%%
% System 4 X3
soundsc(x3,fs3)
soundsc(y4_3_L1,fs3) %high pitch like piano key 
soundsc(y4_3_L2,fs3) % guitar string
soundsc(y4_3_L3,fs3) %guitar string with lower chords
%%
% At L=50, the signal is heard as high pitch like piano key as compared to
% the input, at L=100, the sound is similar to a guitar string while at
% L=400, the sound becomes identical to the sound of base guitar string.


%%
%System 4 Input X4
figure;
subplot 411
stem(1:length(x4),x4);
title(' input signal x4')
subplot 412
stem(1:length(y4_4_L1),y4_4_L1);
title('L=50 Output signal')
subplot 413
stem(1:length(y4_4_L2),y4_4_L2);
title('L=100 Output signal')
subplot 414
stem(1:length(y4_4_L3),y4_4_L3);
title('L=400 Output signal')
%%
% System 4 X4
soundsc(x4,fs4)
soundsc(y4_4_L1,fs4)
soundsc(y4_4_L2,fs4)
soundsc(y4_4_L3,fs4)
%%
% At L=50, the signal is heard as high pitch like piano key as compared to
% the input, at L=100, the sound is similar to a guitar string while at
% L=400, the sound becomes identical to the sound of base guitar string.


%%
%System 4 Input X5
figure;
subplot 411
stem(1:length(x5),x5);
title(' input signal x5')
subplot 412
stem(1:length(y4_5_L1),y4_5_L1);
title('L=50 Output signal')
subplot 413
stem(1:length(y4_5_L2),y4_5_L2);
title('L=100 Output signal')
subplot 414
stem(1:length(y4_5_L3),y4_5_L3);
title('L=400 Output signal')
%%
% System 4 X5
soundsc(x5,fs5)
soundsc(y4_5_L1,fs5) % high pitch alarm
soundsc(y4_5_L2,fs5) %train horn
soundsc(y4_5_L3,fs5) %ship horn
%%
% The audio is transformed from a rain like sound to a high pitch alarm at
% L=50. At L=100, the sound is intensified making it comparable to a train
% horn, and finally at L=400, it seems like a loud ship horn.


%% 
% 4.Answer the following questions (Use MATLAB to support your answer):
% 
% a. Is the system BIBO stable? (Hint: use impz ())

%%
%# System 1
as1 = [1]; % output coefficient
bs1 = [0.5 0.5]; % system 1 input coef
figure(); 
zplane(bs1,as1); % system 1 Z pole
title('System 1 Z-Pole Graph')

s1N=1000;
s1n=0:s1N-1;
s1x = (s1n==0);
s1y=filter(bs1,as1,s1x);
figure(); 
stem(s1n,s1y);

%%
% System 1 is a BIBO stable system as poles are inside the boundary of the system as shown in the zplane.
% Also, impulse response graph of system 1 shows that it is finite thus, the values are bounded and does
% not change in value.

%%
%# System 2
as2 = [1 2 2]; % system 2 output coef
bs2 = [1]; % system 2 input coef
figure();
zplane(bs2,as2); % generate z-pole of system 2;
title('System 2 Z-Pole Graph')

s2N=1000;
s2n=0:s2N-1;
s2x = (s2n==0);
s2y=filter(bs2,as2,s2x);
figure();
stem(s2n,s2y);

%%
% System 2 is NOT a BIBO stable system with its poles are outside of the system's circle as shown in the
% zplane. This is supported by the impulse response graph with its amplitude changes by near end of time.

%%
%# System 3
as3 = [1 2 2]; % system 3 output coef
bs3 = [1.5 0.5]; % system 3 input coef
figure(); 
zplane(bs3,as3); % generate z-pole of system 3;
title('System 3 Z-Pole Graph')

s3N=1000;
s3n=0:s3N-1;
s3x = (s3n==0);
s3y=filter(bs3,as3,s3x);
figure();
stem(s3n,s3y);

%%
% System 3 is NOT a BIBO stable system as its poles are outside of the system's circle as shown in the zplane.
% Similar to system 2, its impulse reponse graph shows that its amplitude chanegs by near end of time

%%
%# System 4
L=50
a_l = zeros(1,L-1);  %
as4 = [1 a_l -0.5 -0.5] % output coefficient n, n-L, n-L-1
bs4 = [1];  % input coefficient n
figure(); 
zplane(bs4,as4);
title('System 4 Z-Pole Graph at L=50')

s4N=1000;
s4n=0:s4N-1;
s4x = (s4n==0);
s4y=filter(bs4,as4,s4x);
figure(); 
stem(s4n,s4y);
title('System 4 at L=50')

L=100
a_l = zeros(1,L-1);  %
as4 = [1 a_l -0.5 -0.5] % output coefficient n, n-L, n-L-1
bs4 = [1];  % input coefficient n
figure(); 
zplane(bs4,as4);
title('System 4 Z-Pole Graph at L=100')

s4N=1000;
s4n=0:s4N-1;
s4x = (s4n==0);
s4y=filter(bs4,as4,s4x);
figure(); 
stem(s4n,s4y);
title('System 4 at L=100')

L=400
a_l = zeros(1,L-1);  %
as4 = [1 a_l -0.5 -0.5] % output coefficient n, n-L, n-L-1
bs4 = [1];  % input coefficient n
figure(); 
zplane(bs4,as4);
title('System 4 Z-Pole Graph at L=400')

s4N=1000;
s4n=0:s4N-1;
s4x = (s4n==0);
s4y=filter(bs4,as4,s4x);
figure(); 
stem(s4n,s4y);
title('System 4 at L=400')

%%
% System 4 is a BIBO stable system as poles are inside the endge of the 
%boundary as shown in the zplane.

%%
% b. Is the system causal?

%%
% # System 1
% System 1 is causal as its output depends on the current and previous
% inputs. From the provided matlab implementation below, the values of y(n)
% depends on the past inputs  n-1.
%
%   function y=dt_1(x)
%   y= zeros(1,length(x)); % zero initialization
%   for n=1:length(x);
%       if n<2
%           y(n) = x(n);
%       else
%           y(n) = 0.5*x(n) + 0.5*x(n-1)
%       end
%   end
%   end

%%
% # System 2
% System 2 is a causal system as its output depends only on the past and
% present values of inputs and previous outputs as presented on the matlab
% implementation below.
% function y=dt_2(x)
%   y= zeros(1,length(x));
% 
%   for n=1:length(x);
%       if n==1
%           y(n) = x(n);
%       elseif n==2
%           y(n) = x(n) - 2*y(n-1);
%       else
%           y(n) = x(n) - 2*y(n-1) - 2*y(n-2)
%       end
%   end
%   end

%%
% # System 3
% System 3 is a causal system as its output depends on the past and present
% values of inputs and previous output( inputs n, n-1, outputs n-1, n-2).
% The matlab implementation of system 3 is presented below.
%   function y=dt_3(x)
%   y= zeros(1,length(x));
%   b = [1.5 0.5];  %input coefficients n , n-1
%   a = [1 2 2]; % output coefficients n, n-1, n-2
%   h = impz(b,a); %impulse response
%   y = conv(h,x);
%   end

%%
% System 4 is also a causal system with its input is dependent on current
% and past inputs and past outputs. Its matlab implementation is presented
% below.
%   function y=dt_4(x,L)
%   y= zeros(1,length(x));
%   b = [1];  % input coefficient n
%   a_l = zeros(1,L-1);  %
%   a = [1 a_l -0.5 -0.5] % output coefficient n, n-L, n-L-1
%   h = impz(b,a); %impulse response
%   y = conv(h,x);
%   end

%%
% # Is the system FIR or IIR?
%
% System 1: Since system 1 is causal and its pole is located at the
% boundary of the Z-pole graph(See the System 1 Z-Pole Graph) at X=-1, 
% therefore system 1 is IIR.
%
% System 2: Since system 2 is causal and its poles are beyond the
% boundary of the system nor the origin.
%
% System 3: Since system 3 is causal and its poles are outside the 
% boundary of the system nor the origin, it is an IIR.
%
% System 4: Since system 4 is causal and its poles are not in the origin
% but within the boundary of the system, it is an IIR.

%%
% # What does the system do? (Based from #2&#3, you may also use freqz(b,a)
% to have an insight on the filter / frequency response)

%%
% System 1
figure(); freqz(bs1,as1);

%%
% OBSERVATION: System 1 shifted the input signal from 0 to -50dB with a 
% phase of 0 to -90 degrees per cycle.


%%
% System 2
figure(); freqz(bs2,as2);

%%
% OBSERVATION: Sysem 2 shifted the magnitude of the input signal from -14 
% to 3dB with a phase of 0 to
% 350 degrees per cycle. This response makes the input audio smoother as
% compared to the input signal.
%%
% System 3
figure(); freqz(bs3,as3);

%%
% OBSERVATION: System 3 shifted the magnitude of the input signal from -8dB 
% to 5dB with a phase of 0 to 360 degrees per cycle. This shift makes the
% input audio smoother as compared to system 2.


%%
% System 4
L=50
a_l = zeros(1,L-1);  %
as4 = [1 a_l -0.5 -0.5] % output coefficient n, n-L, n-L-1
bs4 = [1];  % input coefficient n
figure(); freqz(bs4,as4);
%%
% OBSERVATION: System 4 has a fluctuating frequency response, which
% explains the inaudible noise heard at various L values of the system.
