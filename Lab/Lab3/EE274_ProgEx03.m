%% Paul Vincent S. Nonat  2018-21366   EE274_ProgEx03
%
%
%
%% A.1-2. The Bilateral Z-Transform
%% Sequence (a) $x(n) = (\frac{4}{3})^n u(1-n)$
% *Manual Solution*
%%
% $x(n) = (\frac{4}{3})^n u(-n+1)$
%%
% $X(z) = \sum_{n=-\infty}^{\infty} x(n)z^{-n}$
%%
% $X(z) = \sum_{n=-\infty}^{\infty} (\frac{4}{3})^n u(-n+1)z^{-n}$
%%
% *$Let \ k=-n+1 \ and \ n=1-k$*
%%
% Plugging in the value of $k$ and $n$, we have:
% $X(z) = \sum_{n=-\infty}^{\infty} (\frac{4}{3})^{1-k} u(k)z^{k-1}$
%%
% $X(z) = \sum_{n=0}^{\infty} (\frac{4}{3}) \cdot ((\frac{4}{3})^{-1})^{k} \cdot ((1/z)^{-1})^{k} \cdot z^{-1}$
%%
% $X(z) = (\frac{4z^{-1}}{3}) \ \sum_{n=0}^{\infty} (\frac{3}{4z^{-1}})^{k}$
%%
% $X(z) = {(\frac{4z^{-1}}{3}) \cdot (\frac{1}{1 - \frac{3}{4z^{-1}}}), \ 0 \ < \ \mid {z} \mid \ < \ \frac{4}{3}}$
%%
% $X(z) ={\frac{16z^{-2}}{-9+12z^{-1}}, \ 0 \ < \ \mid {z} \mid \ < \ \frac{4}{3}}$

%%
% *z-plane for 1.(a)*
a_a_coef=[-9, 12, 0];
a_b_coef=[0, 0, -16];
zplane(a_b_coef,a_a_coef);

%%
% *Verification of z-transform v. original sequence with first 8-coef.*
[delta_n,n]= impseq(0,0,7);
A_a_Xz=filter(a_b_coef,a_a_coef,delta_n) %A_a_Xz is z-transform sequence
A_a_Xn=[(4/3).^n].*stepseq(1,0,7) 
figure;
subplot 211
stem(n,A_a_Xz)
title('A.a X(n) First 8-coefficient')
xlabel('n')
subplot 212
stem(n,A_a_Xz)
title('A.a X(Z) First 8-coefficient')
xlabel('n')
%A_a_Xn is the original sequence.

%%
% *The coefficient values generated from $X(z)$ and $x(n)$ are the same.
% Therefore, the z-transform operation for sequence(a) is correct. *
% 

%% Sequence (b) $x(n) = 2^{- \mid {n} \mid} + (\frac{1}{3})^{\mid {n} \mid}$
%%
% $X(z) = \sum_{n=0}^{\infty} 2^{-n}z^{-n} + \sum_{n=0}^{\infty} (\frac{1}{3})^{n}z^{-n}$
%%
% $X(z) = \sum_{n=0}^{\infty} (\frac{z^{-1}}{2})^{n} + \sum_{n=0}^{\infty} (\frac{z^{-1}}{3})^{n}$
%%
% $X(z) = \frac{1}{1-\frac{z^{-1}}{2}}+\frac{1}{1-\frac{z^{-1}}{3}}$
%%
% $X(z) = \frac{2}{2-z^{-1}}+\frac{3}{3-z^{-1}}$
%%
% $X(z) = \frac{12-5z^{-1}}{(2-z^{-1})(3-z^{-1})},\ \mid z \mid \ > \ \frac{1}{3} \ \cap \ \mid z \mid \ > \ \frac{1}{2}$
%%
% $ X(z) = \frac{12-5z^{-1}}{6-5z^{-1}+z^{-2}},\ \mid z \mid \ > \ \frac{1}{3} \ \cap \ \mid z \mid \ > \ \frac{1}{2}$
%%
% *z-plane for 1.(b)*
b_a_coef=[6 -5 1];
b_b_coef=[12 -5 0];
zplane(b_b_coef,b_a_coef);

%%
% *Verification of z-transform v. original sequence with first 8-coef.*
[delta,n]= impseq(0,0,7);
A_b_Xz=filter(b_b_coef,b_a_coef,delta) %A_b_Xz is z-transform sequence
A_b_Xn=((2).^(-abs(n)))+((1/3).^(abs(n))) %A_b_Xn is the original sequence
figure;
subplot 211
stem(n,A_b_Xn)
title('A.b X(n) First 8-coefficient')
xlabel('n')
subplot 212
stem(n,A_b_Xz)
title('A.b X(Z) First 8-coefficient')
xlabel('n')
%%
% *The first 8-coefficient values from $X(z)$ and $x(n)$ are the same
% Therefore, the z-transform operation for sequence(b) is correct. *

%% A.3. $x(n)=(\frac{1}{3})^{n}u(n-2)+(0.9)^{n-3}u(n)$
%%
% $X(z)={\frac{3z^{-2}}{27-9z^{-1}}}+{\frac{1.3717}{1-0.9z^{-1}}}$
%%
% $X(z)={\frac{{37.036}-{12.34z^{-1}}+{3z^{-2}}-{2.7z^{-3}}}{{27}-{33.3z^{-1}}+{8.1z^{-2}}}} 
% \ \mid z \mid \ > \ \frac{1}{3} \ \cap \ \mid z \mid \ > \ {0.9}$

%%
% *z-plane for A.3*
A3_b_coef=[37.036, -12.34, 3, -2.7];
A3_a_coef=[27, -33.3, 8.1];
zplane(A3_b_coef,A3_a_coef);

%%
% *Verification of z-transform vs original sequence using the first 20-coef.*
[delta,n]= impseq(0,0,19);
A3_Xz=filter(A3_b_coef,A3_a_coef,delta) %A3_Xz is z-transform sequence
A3_Xn=(((1/3).^n).*(stepseq_n(2,0,19))+(((0.9).^(n-3)).*(stepseq_n(0,0,19)))) 
%A3_Xn is the original sequence, see stepseq_n.m
figure;
subplot 211
stem(n,A3_Xn)
title('X(n) First 20-coefficient')
xlabel('n')
subplot 212
stem(n,A3_Xz)
title('X(Z) First 20-coefficient')
xlabel('n')
%%
% *The first 20-coefficient values from $X(z)$ and $x(n)$ are the same
% Therefore, the z-transform operation for sequence(A 3) is correct. *
 
%% B.4. Inverse Z-Transform
%% Sequence(c) $X(z)={\frac{1-z^{-1}-4z^{-2}+4z^{-3}}{1-\frac{11}{4}z^{-1}+\frac{13}{8}z^{-2}-\frac{1}{4}z^{-3}}}$
B4_b_coef=[1, -1, -4, 4];
B4_a_coef=[1, (-11/4), (13/8), (-1/4)];
[B4_zeros, B4_poles, B4_C]=residuez(B4_b_coef,B4_a_coef);
%%
% $X(z)=\frac{0z}{z-2} - \frac{10z}{z-0.5} + \frac{27z}{z-0.25} - {16}$
%%
% $X(n)=u(-n)-(2^{-2n}(5 \times 2^{n+1}-27)(1-u(-n)))$
%%
% *Verification of z-transform v. ans sequence with first 8-coef.*
%%
% 
[delta,n]= impseq(0,0,8);
B4_Xz=filter(B4_b_coef,B4_a_coef,delta); %B4_Xz is z-transform sequence
%B4_Xn is inv. ztrans sequence
B4_Xn=-heaviside(-n)-((2.^(-2*n)).*(5.*(2.^(n+1))-27).*(1-heaviside(-n)));
B4_Xz(2:8)% First 8 coef of B4_Xz - Z-transf 
B4_Xn(2:8)% First 8 coef of B4_Xn - Inv. Z-transf
figure;
title('Inverse Z-Transform')
subplot 211
stem(n,B4_Xz)
title('X(Z) First 8-coefficient')
xlabel('n')
subplot 212
stem(n,B4_Xn)
title('X(n) First 8-coefficient')
xlabel('n')

%% C.5. Signal Generation
%%
% _Generate the periodic even symmetric square pulse x(n) from [0, 1].
% The period of the pulse is 1 second and a pulse with of 250ms with sampling freq.
% of 8KHz. Plot one period of x(n) and verify if you have the correct
% waveform._ 

period = 1
pulse_width=0.250
fs=8000
C5_time=0:period/fs:period; % time from frequenzy 8kHz
C5_x=square((2*pi*C5_time),(pulse_width/(period*2))*100); % generate x(n)

C5_x(end)=[]
C5_time(end)=[]
C5_x=(abs(C5_x)+C5_x)/2; % remove -1 samples to make x(n) [0,1]
C5_x = C5_x + flip(C5_x)

figure;
plot(C5_time, C5_x); 
title("x(n)"); % 2 periods w/ 250ms pw each 1,0.
xlabel("Time in second (s)");
ylabel("Amplitude");

%% 5.a. How many samples in one period?
%% 
sampled_period=(length(C5_x)) % samples in one period

%%
% *Since the Sampling frequency is 8kHz, therefore there are 8000 samples in
% one period.*

%% 5.b. How many samples with a value of 1?
% 
value1_period=(sum(C5_x(:)==1)) % in one period

%%
% *Since the generated signals has a period of 1 second and a pulse width
% of 0.250 seconds, thereforce 25% of the sampled datapoints contains a
% value of one which is equal to 2000.*

%% 5.c. How many zeros?
value0_period=sum(C5_x(:)==0) % in one period

%%
% *Since the generated signals has a period of 1 second and a pulse width
% of 0.250 seconds, thereforce 75% of the sampled datapoints contains a
% value of zero which is equal to 6000.*
figure;
bar(categorical({'samples'}),[value0_period;value1_period])
title("Samples with a value of 0")
ylabel('Number of Samples')
legend('Sampled Zeros','Sampled Ones')
%% C.6.Fourier Series Analysis Equation
%%
% _Using the analysis equation of the Fourier series, write a program
% that will compute the Fourier series coefficients of the periodic square
% pulse signal. Plot the magnitude and phase of the first 10 Fourier coef._

%% 
% plot of 100000 harmonic
C6_t = (0:1/8000:1);
C6_NHarmonics=20; 
cycles=1; 
C6_Nsamples=8000;
C6_y(1:C6_Nsamples)=pulse_width;
C6_j=1:C6_Nsamples;
for C6_k=1:C6_NHarmonics
 C6_x(C6_j)=(2*sin(0.25*pi*C6_k)/(pi*C6_k))*cos(C6_k*2*pi*cycles*C6_j/C6_Nsamples);
 C6_y=C6_y+C6_x;
end

C6_NHarmonics=10; 
cycles=1; 
C6_Nsamples=8000;
C6_y_harmonic(1:C6_Nsamples)=pulse_width;
C6_j=1:C6_Nsamples;
for C6_k=1:C6_NHarmonics
 C6_x(C6_j)=(2*sin(0.25*pi*C6_k)/(pi*C6_k))*cos(C6_k*2*pi*cycles*C6_j/C6_Nsamples);
 C6_y_harmonic=C6_y_harmonic+C6_x;
end

figure()
subplot 211
plot(C6_t(1:8000),C6_y);
title("Illustration of the Fourier synthesis with 20 harmonics");
xlabel("Time in second (s)"); 
ylabel("Amplitude");
subplot 212
plot(C6_t(1:8000),C6_y_harmonic)
title("Illustration of the Fourier synthesis with 10 harmonics");
xlabel("Time in second (s)");
ylabel("Amplitude");
% plot of 20 harmonic and 10 harmonic
%% 
% compute the magnitude and phase of 10 harmonics signal
C6y_fft = fft(C6_y_harmonic);
C6_magnitude = abs(C6y_fft);
C6_phaseangle = angle(C6y_fft);

figure()
subplot 211
stem(C6_magnitude(1:10)); 
title("Magnitude of first 10 Coef"); 
xlabel("coef"); 
ylabel("Magnitude");

subplot 212
stem(C6_phaseangle(1:10)); 
title("Phase of first 10 Coef"); 
xlabel("coef"); 
ylabel("Phase");

%% 6.a. What is the fundamental frequency of the square pulse?
%%
% *The fundamental frequency is defined by f=1/T. which indicates that
% one complete period is 1s with frequency =1Hz wherein 250ms is alloted for 
% on and 750ms for off.*

%% 6.b. Enumerate the Magnitude and Phase of first 10 coef.
disp("Magnitude");
C6_magnitude(1:10)

disp("Phase");
C6_phaseangle(1:10)

%% C.7. Fourier Series Synthesis Equation
%%
% _Using the synthesis equatin for the Foyrier series, synthesiez the
% original square pulse using the first 10 Fourier coefficients. Generate a
% plot of the original square pulse and the synthesized square pulse._ 

figure()
plot(C6_t(1:8000),C6_y_harmonic, 'color', 'r'); 
hold on; plot(C6_t(1:8000),C5_x(1:8000), 'color', 'b'); 
hold off; 
title("Comparison of Original and Harmonic Signal"); 
xlabel("Time in seconds (s)"); 
ylabel("Amplitude");
legend("10 Harmonic Square Signal","Original Signal")
%% 7.a. What is the average MSE of original square pulse vs synthesized pulse?
%
C7_mse_10harm = immse(C5_x(1:8000),C6_y_harmonic)
% *MSE is 0.92%*
%% 7.b. If you use 20 Fourier coef, what will be the MSE?
%
%
C6_NHarmonics=20; 
cycles=1; 
C6_Nsamples=8000;
C6_y_harmonic20(1:C6_Nsamples)=pulse_width;
C6_j=1:C6_Nsamples;
for C6_k=1:C6_NHarmonics
 C6_x(C6_j)=(2*sin(0.25*pi*C6_k)/(pi*C6_k))*cos(C6_k*2*pi*cycles*C6_j/C6_Nsamples);
 C6_y_harmonic20=C6_y_harmonic20+C6_x;
end

figure()
plot(C6_t(1:8000),C6_y_harmonic20, 'color', 'r'); 
hold on; plot(C6_t(1:8000),C5_x(1:8000), 'color', 'b'); 
hold off; 
title("Comparison of Original and Harmonic Signal"); 
xlabel("Time in seconds (s)"); 
ylabel("Amplitude");
legend("20 Harmonic Square Signal","Original Signal")

C7_mse_20harm = immse(C5_x(1:8000),C6_y_harmonic20)
% *MSE will be 0.51%*

%% 7.c. What is the effect on the fundamental freq if I increase the pulse width to 300ms? Explain.
% 
period = 1
pulse_width=0.3
fs=8000
C5_time=0:period/fs:period; % time from frequenzy 8kHz
C5_x=square((2*pi*C5_time),(pulse_width/(period*2))*100); % generate x(n)

C5_x(end)=[]
C5_time(end)=[]
C5_x=(abs(C5_x)+C5_x)/2; % remove -1 samples to make x(n) [0,1]
C5_x = C5_x + flip(C5_x)

figure;
plot(C5_time, C5_x); 
title("x(n) at PW=0.3s"); % 2 periods w/ 250ms pw each 1,0.
xlabel("Time in second (s)");
ylabel("Amplitude");
%%
% The fundamental frequency stays the same given that the increase in pulse
% width only increase the duty cycle of the signal, while the period stays
% the same.

%% 7.d. What is the effect on the Fourier coefficients if I change the pulse width
% *Depending on the change applied to the pulse width, the fourier
% coefficient will be shifted proportional to the change in the pulse width
% of the square wave.

%% 
%
%
C6_NHarmonics=10; 
cycles=1; 
C6_Nsamples=8000;
C6_y_harmonic_new_width(1:C6_Nsamples)=pulse_width;
C6_j=1:C6_Nsamples;
width=0.3
for C6_k=1:C6_NHarmonics
 C6_x(C6_j)=(2*sin(width*pi*C6_k)/(pi*C6_k))*cos(C6_k*2*pi*cycles*C6_j/C6_Nsamples);
 C6_y_harmonic_new_width=C6_y_harmonic_new_width+C6_x;
end

figure()
plot(C6_t(1:8000),C6_y_harmonic_new_width)
title("Illustration of the Fourier synthesis with 10 harmonics at Pw=0.3s");
xlabel("Time in second (s)");
ylabel("Amplitude");

%% 
% 
C6y_fft_newwidth = fft(C6_y_harmonic_new_width);
C6_magnitude_newwidth = abs(C6y_fft_newwidth);
C6_phaseangle_newwidth = angle(C6y_fft_newwidth);

figure()
subplot 411
stem(C6_magnitude(1:10)); 
title("Magnitude of first 10 Coef"); 
xlabel("coef"); 
ylabel("Magnitude");

subplot 412
stem(C6_phaseangle(1:10)); 
title("Phase of first 10 Coef"); 
xlabel("coef"); 
ylabel("Phase");
subplot 413
stem(C6_magnitude_newwidth(1:10)); 
title("Magnitude of first 10 Coef at PW=0.3s"); 
xlabel("coef"); 
ylabel("Magnitude");

subplot 414
stem(C6_phaseangle_newwidth(1:10)); 
title("Phase of first 10 Coef at PW=0.3s"); 
xlabel("coef"); 
ylabel("Phase");
%% 7.e. What is the effect on the Fourier coefficients if I change the period?
% *Depending on the change applied to the values of the period, the changes
% will alter the values of Fourier coefficient inversely.

%% 
%
%
C6_NHarmonics=10; 
cycles=2; 
C6_Nsamples=8000;
C6_y_harmonic_new_width(1:C6_Nsamples)=pulse_width;
C6_j=1:C6_Nsamples;
width=0.25
for C6_k=1:C6_NHarmonics
 C6_x(C6_j)=(2*sin(width*pi*C6_k)/(pi*C6_k))*cos(C6_k*2*pi*cycles*C6_j/C6_Nsamples);
 C6_y_harmonic_new_width=C6_y_harmonic_new_width+C6_x;
end

figure()
plot(C6_t(1:8000),C6_y_harmonic_new_width)
title("Illustration of the Fourier synthesis with 10 harmonics");
xlabel("Time in second (s)");
ylabel("Amplitude");

%% 
% 
C6y_fft_newwidth = fft(C6_y_harmonic_new_width);
C6_magnitude_newwidth = abs(C6y_fft_newwidth);
C6_phaseangle_newwidth = angle(C6y_fft_newwidth);

figure()
subplot 411
stem(C6_magnitude(1:10)); 
title("Magnitude of first 10 Coef"); 
xlabel("coef"); 
ylabel("Magnitude");

subplot 412
stem(C6_phaseangle(1:10)); 
title("Phase of first 10 Coef"); 
xlabel("coef"); 
ylabel("Phase");
subplot 413
stem(C6_magnitude_newwidth(1:10)); 
title("Magnitude of first 10 Coef at T=2"); 
xlabel("coef"); 
ylabel("Magnitude");

subplot 414
stem(C6_phaseangle_newwidth(1:10)); 
title("Phase of first 10 Coef at T=2"); 
xlabel("coef"); 
ylabel("Phase");
%% Functions used in the activity to generate sequence
% 

function [x,n]=impseq(n0,n1,n2)
    n = [n1:n2];
    x = [(n-n0) == 0];
end

%%
%
function [x,n]=stepseq(n0,n1,n2)
    n = [n1:n2];
    x = [(n0-n) < 0]; %change to less than to satisfy u(n-1)condition.
end

%%
%
function [x,n]=stepseq_n(n0,n1,n2)
    n = [n1:n2];
    x = [(n0-n) <= 0]; %change to less than to satisfy (n)condition.
end