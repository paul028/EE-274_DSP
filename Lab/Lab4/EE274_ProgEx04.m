%% Paul Vincent S. Nonat 2018-21366 EE274_ProgEx04
% # Part 1 - EE 274: Tutorials for Discrete Fourier Transform
% # Part 2 - EE 274 / CoE 197E Lab Exercise: Properties of the DFT
%

%% Part 1 - Tutorials for Discrete Fourier Transform
%
%
%% A. Resolving two signals of close frequency
% A.1-2
%
fs = 4000;
freq = [100 150];
amp = [1.2 0.8];
ph = [0 0];
t = [0:1/fs:5/min(freq)];
% generate the sinusoid
y =((amp'*ones(1,length(t))).*sin((2*pi*freq'*t)+ph'*ones(1,length(t))))';
figure(6); plot(sum(y'));

ty_sum = sum(y');
y_sum = ty_sum(1:100);
Nfft = 100;
y_fft = fft(y_sum, Nfft);
% generate the frequency axis
w = linspace(0,fs,Nfft);
% Plot the magnitude and the phase response
mag = abs(y_fft);
ph1 = angle(y_fft);
figure(7); subplot(2,1,1), stem(w, mag, '^');
axis ([0, 200, 0, 140]);
% use axis to zoom into plot to determine frequencies of sinusoids
title('The magnitude response from 0 to fs');
axis; % turn autoscaling back on
subplot(2,1,2), stem(w,ph1, '^');
title('The phase response from 0 to fs');

%%
% Question is pre-answered in the tutorial pdf file.
%
%% A.3.a-b.
%
fs = 4000;
freq = [400 150];
amp = [1.2 0.8];
ph = [0 0];
t = [0:1/fs:5/min(freq)];
% generate the sinusoid
y =((amp'*ones(1,length(t))).*sin((2*pi*freq'*t)+ph'*ones(1,length(t))))';
figure(); plot(sum(y'));

ty_sum = sum(y');
y_sum = ty_sum(1:100);
Nfft = 400;
y_fft = fft(y_sum, Nfft);
% generate the frequency axis
w = linspace(0,fs,Nfft);
% Plot the magnitude and the phase response
mag = abs(y_fft);
ph1 = angle(y_fft);
figure(); subplot(2,1,1), stem(w, mag, '^');
axis ([0, 200, 0, 140]);
% use axis to zoom into plot to determine frequencies of sinusoids
title('The magnitude response from 0 to fs');
axis; % turn autoscaling back on
subplot(2,1,2), stem(w,ph1, '^');
title('The phase response from 0 to fs');

%%
% # Can we now correctly determine the frequencies of the sinusoids?
% # What is the effect of increasing the number of frequency samples on the
% accuracy of the DFT?

%%
% *ANSWER*
% *Increasing the number of samples of the given sinusoids will allow us to
% correctly determine its frequencies with more certainty as it increases
% the number of phases and magnitudes of the DFT.*

%% A.4.
fs = 4000;
freq = [100 150];
amp = [1.2 0.8];
ph = [0 0];
t = [0:1/fs:2/min(freq)];
% generate the sinusoid
y =((amp'*ones(1,length(t))).*sin((2*pi*freq'*t)+ph'*ones(1,length(t))))';
figure(); plot(sum(y'));

ty_sum = sum(y');
y_sum = ty_sum();
Nfft = 100;
y_fft = fft(y_sum, Nfft);
% generate the frequency axis
w = linspace(0,fs,Nfft);

% Plot the magnitude and the phase response
mag = abs(y_fft);
ph1 = angle(y_fft);
figure(); subplot(2,1,1), stem(w, mag, '^');
axis ([0, 200, 0, 140]);
% use axis to zoom into plot to determine frequencies of sinusoids
title('The magnitude response from 0 to fs at N=100');
axis; % turn autoscaling back on
subplot(2,1,2), stem(w,ph1, '^');
title('The phase response from 0 to fs at N=100');

fs = 4000;
freq = [100 150];
amp = [1.2 0.8];
ph = [0 0];
t = [0:1/fs:2/min(freq)];
% generate the sinusoid
y =((amp'*ones(1,length(t))).*sin((2*pi*freq'*t)+ph'*ones(1,length(t))))';
figure(); plot(sum(y'));

ty_sum = sum(y');
y_sum = ty_sum();
Nfft = 400;
y_fft = fft(y_sum, Nfft);
% generate the frequency axis
w = linspace(0,fs,Nfft);

% Plot the magnitude and the phase response
mag = abs(y_fft);
ph1 = angle(y_fft);
figure(); subplot(2,1,1), stem(w, mag, '^');
axis ([0, 200, 0, 140]);
% use axis to zoom into plot to determine frequencies of sinusoids
title('The magnitude response from 0 to fs at N=400');
axis; % turn autoscaling back on
subplot(2,1,2), stem(w,ph1, '^');
title('The phase response from 0 to fs at N=400');

%%
% *As shown in the results above, the DFT can resolve the frequencies of
% the two signals. However with N=100, it cannot accurately resolve the
% frequencies due to the limited samples of magnitude and phase spectrum as
% compared to N=400, which contains more samples of phase and magnitude,
% thereby allowing us to resolve the frequency with better accuracy.

%% B. Effects of zero padding and signal length on the DFT
%% B.1.
fs = 4000;
freq = 1000;
amp = 1;
ph = 0;
t = [0:1/fs:2/freq];
% generate the sinusoid
y = amp*sin((2*pi*freq*t)+ph)';
y_sum = y;
Nfft = 100;
y_fft = fft(y_sum, Nfft);
% generate the frequency axis
w = linspace(0, fs/2, Nfft/2);
y_fft = y_fft(1:length(y_fft)/2);
% Plot the magnitude response,
mag = (y_fft);
figure(); plot(w, mag, '^');
title('The magnitude response from 0 to fs/2');

%% B.1.a.
% *The plot of the frequency response is similar to the sinc function
% because of the abs function applied to the sinusoid, thereby transporming
% its negative values into positive. As such, the sum values from 500:1500
% samples yields to a wave with a max magnitude of 4 as shown below
%

figure(); plot(w(1:25), y_fft(1:25));
figure(); plot(w(26:50), y_fft(26:50));
figure(); plot(w,abs(y_fft));

%% B.1.b.
% *Increasing the number of frequency samples is equivalent to having a more
% accurate magnitude response and its plot.*

%% B.1.c.
% *The same result as with B.1.b for increasing the number of samples.
%  Increasing the period means increasing the amplitude as show in the
%  generated figure.
fs = 4000;
freq = 1000;
amp = 1;
ph = 0;
t = [0:1/fs:5/freq]; % period inc to 5 from 2.
% generate the sinusoid
y = amp*sin((2*pi*freq*t)+ph)';
y_sum = y;
Nfft = 10000; %Increase samples to 10K
y_fft = fft(y_sum, Nfft);
% generate the frequency axis
w = linspace(0, fs/2, Nfft/2);
y_fft = y_fft(1:length(y_fft)/2);
% Plot the magnitude response,
mag = abs(y_fft);
figure(); plot(w, mag, '^');
title('The magnitude response from 0 to fs/2');

%% B.1.d.
% *Increasing the signal length means extending the range of lobes in the
% time domain. Therefore, if we will not modify the time range, the plot
% will not be seeon in the default range due to the increase in length.*

%% B.2.
fs = 4000;
freq = 1000;
amp = 1;
ph = 0;
t = [0:1/fs:2/freq];
% generate the sinusoid
y = amp*sin((2*pi*freq*t)+ph)';
y_sum = y;
y_fft = fft(y_sum);
y_fft = y_fft(1:length(y_fft)/2);
% Plot the magnitude response,
mag = abs(y_fft);
figure(); plot(mag);
title('The magnitude response from 0 to fs/2');
%% B.2.a.
% *As per the result of the code in this section and the fft() documatation,
% without specifying the sample number, function will return the fft as a
% vector of y_sum. As with the magnitude response, only the linear peak 
% values of magnitude of the original plot.*

%% B.2.b.
fs = 4000;
freq = 1000;
amp = 1;
ph = 0;
t = [0:1/fs:2/freq];
% generate the sinusoid
y = amp*sin((2*pi*freq*t)+ph)';
y(length(y)+1:100)=0; % 100 zero-padding
y_fft = fft(y);
y_fft = y_fft(1:length(y_fft)/2);
% Plot the magnitude response,
mag = abs(y_fft);
figure(); plot(mag);
title('The magnitude response from 0 to fs/2');
%%
%Based on the results show in the figure, padded zeros in the DFT acts
%similar to the frequency sample in the DFT thus, making the signal revert
%to the sinc like plot.

%% B.2.c.
fs = 4000;
freq = 1000;
amp = 1;
ph = 0;
t = [0:1/fs:5/freq]; %change period to 5
% generate the sinusoid
y = amp*sin((2*pi*freq*t)+ph)';
% y(length(y)+1:100)=0; % 100 zero-padding
y_fft = fft(y);
y_fft = y_fft(1:length(y_fft)/2);
% Plot the magnitude response,
mag = abs(y_fft);
figure(); plot(mag);
title('The magnitude response from 0 to fs/2');


%%
% *Increasing the signal period means increasing the magnitude of the peaks
% of the signal. Furthermore, the increase gradually eliminate the
% sidelobes in the magnitude spectrum.*


%% C 
% Write a Matlab code that will generate the presented figures Figure A and
% Figure B(see EE 274: Tutorials for Discrete Fourier Transform tutorial
% PDF file).

P=2; 
W=2; 
t=0:0.01:P;
R = rectpuls(t, W);
f=-5:0.01:5;
XR=exp(-j*pi*f).*sinc(f);
figure()
subplot 211
plot(t, R);
title("Time function: square pulse"); 
grid on;
subplot 212
plot(f,abs(XR));
title("Frequency spectrum: Sinc function");
grid on;

%%
%
Sf=8;
D2=1;
N2=8;

tn2=0:1/Sf:(N2-1)*1/Sf;
C2_xn=rectpuls(tn2, W);
ya2=fft(C2_xn, N2);
ya2=fftshift(ya2);
fo2=1/D2;
fa2=-(N2/2)*fo2:fo2:(N2/2-1)*fo2;


D3=2;
N3=16;
tn3=0:1/Sf:(N3-1)*1/Sf;
xn3=rectpuls(tn3, W);
ya3=fft(xn3, N3);
ya3=fftshift(ya3);
fo3=1/D3;
fa3=-(N3/2)*fo3:fo3:(N3/2-1)*fo3;
figure(); 
subplot 211; 
title('Sf=8Hz, D=1, N=8'); 
grid on; 
hold on;
plot(f, abs(XR)); 
stem(fa2,1/Sf*(abs(ya2)));
hold off;

subplot 212; 
title('Sf=8Hz, D=2, N=16'); 
grid on; 
hold on;
plot(f, abs(XR)); 
stem(fa3,1/Sf*(abs(ya3)));
hold off;
%%
%
Sf4=8;
D4=4;
N4=32;
tn4=0:1/Sf4:(N4-1)*1/Sf4;
xn4=rectpuls(tn4, W);
ya4=fft(xn4, N4);
ya4=fftshift(ya4);
fo4=1/D4;
fa4=-(N4/2)*fo4:fo4:(N4/2-1)*fo4;

f4=-8:0.01:8;
XR4=exp(-j*pi*f4).*sinc(f4);

Sf5=16;
D5=4;
N5=64;
tn5=0:1/Sf5:(N5-1)*1/Sf5;
xn5=rectpuls(tn5, W);
ya5=fft(xn5, N5);
ya5=fftshift(ya5);
fo5=1/D5;
fa5=-(N5/2)*fo5:fo5:(N5/2-1)*fo5;

figure();
subplot 211
plot(f, abs(XR)); 
title('Sf=8Hz, D=4, N=32'); 
grid on; 
hold on;
stem(fa4,1/Sf4*(abs(ya4)));
hold off;

subplot 212
plot(f4, abs(XR4));
title('Sf=16Hz, D=4, N=64'); 
grid on; 
hold on;
stem(fa5,1/Sf5*(abs(ya5)));
hold off;

%% Part 2 - Properties of DFT
%% A. Computation of DFT
%
clf;
w = -4*pi:8*pi/511:4*pi;
num = [2 1];den = [1 -0.6];
h = freqz(num, den, w);
x_axis = w/pi

subplot(4,1,1)
plot(x_axis,real(h));grid
title('Real part of H(e^{j\omega})')
xlabel('\omega /\pi');
ylabel('Amplitude');
subplot(4,1,2)
plot(x_axis,imag(h));grid
title('Imaginary part of H(e^{j\omega})')
xlabel('\omega /\pi');
ylabel('Amplitude');
subplot(4,1,3)
plot(x_axis,abs(h));grid
title('Magnitude Spectrum |H(e^{j\omega})|')
xlabel('\omega /\pi');
ylabel('Amplitude');
subplot(4,1,4)
plot(x_axis,angle(h));grid
title('Phase Spectrum arg[H(e^{j\omega})]')
xlabel('\omega /\pi');
ylabel('Phase in radians');

%% A.1. Compare and differentiate above code output and freqz(num,den).
% *The output waveform of freqz(num,den) with a w parameter plots every two
% interval in the x-space. While when we remove the w parameter in the
% freqz function, it maximize the signa plot on the entire -4:4 $w/\pi$
% space.

%% A.2. Symmetries in real and imag parts of DFT.
% *The two graphs are symmetric as they follow the same function, only that
% the real parts plot every 8*pi/511 in -4:4 $w/\pi$ space while the latter
% plots in the entire plane. DFT is a periodic function of omega with a
% period of 2*pi. All of the plots presented above are in the period of
% 2*pi, with the real part and magnitude are even symmetric while the
% imaginary part of phase are odd symmetric.*

%% A.3. Modify the above code and to satisfy the given function below.
%
% $U(e^{j \omega}) = \frac{0.7-0.5e^{j \omega}+0.3e^{-j2 \omega}+e^{-j3 \omega}}{1+0.3e^{-j \omega}-0.5^{-j2 \omega}+0.7e^{-j3 \omega}}$
N_A3 = 512;
g_n3 = [0.7 -0.5 0.3 1]; 
den_A3 = [1 0.3 -0.5 0.7];
[h_A_3, w_A3] = freqz(g_n3, den_A3, N_A3);
x_axis = w_A3/pi

figure()
subplot(4,1,1)
plot(x_axis,real(h_A_3));
grid on;
title('Real part of H(e^{j\omega})')
xlabel('\omega /\pi');
ylabel('Amplitude');
subplot(4,1,2)
plot(x_axis,imag(h_A_3));
grid on;
title('Imaginary part of H(e^{j\omega})')
xlabel('\omega /\pi');
ylabel('Amplitude');
subplot(4,1,3)
plot(x_axis,abs(h_A_3));
grid on;
title('Magnitude Spectrum |H(e^{j\omega})|')
xlabel('\omega /\pi');
ylabel('Amplitude');
subplot(4,1,4)
plot(x_axis,angle(h_A_3));
grid on
title('Phase Spectrum arg[H(e^{j\omega})]')
xlabel('\omega /\pi');
ylabel('Phase in radians');

%% A.4. Modify the code and evaluate the sequence. Use freqz and fft. Compare.
N_A4 = 512
w_A4 = -4*pi:8*pi/511:4*pi;
g_n = [1 3 5 7 9 11 13 15 17]; A4_den = 1;
h_freqz_A4 = freqz(g_n, A4_den, w_A4); % freqz function
hf_fft_A4 = fft(g_n, N_A4); % fft function
x_axis =w_A4/pi

%freqz
figure()
subplot(4,1,1)
plot(x_axis,real(h_freqz_A4));grid
title('Real part of H(e^{j\omega}) using freqz()')
xlabel('\omega /\pi');
ylabel('Amplitude');
subplot(4,1,2)
plot(x_axis,imag(h_freqz_A4));grid
title('Imaginary part of H(e^{j\omega}) using freqz()')
xlabel('\omega /\pi');
ylabel('Amplitude');
subplot(4,1,3)
plot(x_axis,abs(h_freqz_A4));grid
title('Magnitude Spectrum |H(e^{j\omega})| using freqz()')
xlabel('\omega /\pi');
ylabel('Amplitude');
subplot(4,1,4)
plot(x_axis,angle(h_freqz_A4));grid
title('Phase Spectrum arg[H(e^{j\omega})] using freqz()')
xlabel('\omega /\pi');
ylabel('Phase in radians');

% fft
figure()
subplot(4,1,1)
plot(x_axis,real(hf_fft_A4));grid
title('Real part of H(e^{j\omega}) using fft()')
xlabel('\omega /\pi');
ylabel('Amplitude');
subplot(4,1,2)
plot(x_axis,imag(hf_fft_A4));grid
title('Imaginary part of H(e^{j\omega}) using fft()')
xlabel('\omega /\pi');
ylabel('Amplitude');
subplot(4,1,3)
plot(x_axis,abs(hf_fft_A4));grid
title('Magnitude Spectrum |H(e^{j\omega})| using fft()')
xlabel('\omega /\pi');
ylabel('Amplitude');
subplot(4,1,4)
plot(x_axis,angle(hf_fft_A4));grid
title('Phase Spectrum arg[H(e^{j\omega})] using fft()')
xlabel('\omega /\pi');
ylabel('Phase in radians');
%%
% *freqz() and fft() function performs the same operation for DFT. However,
% freqz() compute the phase and magnitude periodically while fft() performs
% the same operation on one period.

%% B. Time shift property of the DFT
clf;
w_B = -pi:2*pi/255:pi; 
D_B = -10; %time shift
B_num = [1 2 3 4 5 6 7 8 9];
h1_B = freqz(B_num, 1, w_B);
h2_B = freqz([zeros(1,D_B) B_num], 1, w_B);

subplot(2,2,1)
plot(w_B/pi,abs(h1_B));grid
title('Magnitude of Original Sequence')
xlabel('\omega /\pi');
ylabel('Amplitude');
subplot(2,2,2)
plot(w_B/pi,abs(h2_B));grid
title('Magnitude of Time-Shifted Sequence')
xlabel('\omega /\pi');
ylabel('Amplitude');
subplot(2,2,3)
plot(w_B/pi,angle(h1_B));grid
title('Phase of Original Sequence')
xlabel('\omega /\pi');
ylabel('Phase in radians');
subplot(2,2,4)
plot(w_B/pi,angle(h2_B));grid
title('Phase of Time-Shifted Sequence')
xlabel('\omega /\pi');
ylabel('Phase in radians');
%%
% * As shown on the generated plots, the magnitude spectrum is not affected
% by the time shift, instead the time shift is applied to the phase
% spectrum. The variable that handles the time shift is D_B. Increasing the
% time shift makes the phase angle oscilllates faster in the the %w/pi% 
% domain.

%% C. Frequency-shift property of DFT
clf;
w = -pi:2*pi/255:pi; 
wo = -0.4*pi; %frequency shifted to the left 
num1 = [1 3 5 7 9 11 13 15 17];
L = length(num1);
h1 = freqz(num1, 1, w);
n = 0:L-1;
num2 = exp(wo*i*n).*num1;
h2 = freqz(num2, 1, w);
subplot(2,2,1)
plot(w/pi,abs(h1));grid
title('Magnitude of Original Sequence')
xlabel('\omega /\pi');
ylabel('Amplitude');
subplot(2,2,2)
plot(w/pi,abs(h2));grid
title('Magnitude of Frequency-Shifted Sequence')
xlabel('\omega /\pi');
ylabel('Amplitude');
subplot(2,2,3)
plot(w/pi,angle(h1));grid
title('Phase of Original Sequence')
xlabel('\omega /\pi');
ylabel('Phase in radians');
subplot(2,2,4)
plot(w/pi,angle(h2));grid
title('Phase of Frequency-Shifted Sequence')
xlabel('\omega /\pi');
ylabel('Phase in radians');
%%
% *The generated plot shows that both the phase and magnitude are shifted
% according to the magnitude of the parameter wo (negative shifts to the
% left, while positive shifts to the right).

%% D. Convolution property of DFT
clf;
w = -pi:2*pi/255:pi;
x1 = [1 3 5 7 9 11 13 15 17];
x2 = [1 -2 3 -2 1];
y = conv(x1,x2); %time domain convolution
h3 = freqz(y,1,w);

h1 = freqz(x1, 1, w); 
h2 = freqz(x2, 1, w); 
hp = h1.*h2; %DFT multiplication

subplot(2,2,1)
plot(w/pi,abs(hp));grid
title('Product of Magnitude Spectra')
xlabel('\omega /\pi');
ylabel('Amplitude');
subplot(2,2,2)
plot(w/pi,abs(h3));grid
title('Magnitude Spectrum of Convolved Sequence')
xlabel('\omega /\pi');
ylabel('Amplitude');
subplot(2,2,3)
plot(w/pi,angle(hp));grid
title('Sum of Phase Spectra')
xlabel('\omega /\pi');
ylabel('Phase in radians');
subplot(2,2,4)
plot(w/pi,angle(h3));grid
title('Phase Spectrum of Convolved Sequence')
xlabel('\omega /\pi');
ylabel('Phase in radians');

%%
% * As shown in the generated plots, the magnitude and phase from the
% multiplication of two DFT from the original sequences are the same to
% those obtained from the time-domain convolution of the said sequence.
% Hence, the convolution operation in time-domain is the same as the
% multiplication in DFT as presented in comparing the results obtained in
% the DFT version of y and the product of DFT multiplication hp.

%% E. Modulation property
clf;
w = -pi:2*pi/255:pi;
x1 = [1 3 5 7 9 11 13 15 17]; 
x2 = [1 -1 1 -1 1 -1 1 -1 1]; 
y = x1.*x2;

h1 = freqz(x1, 1, w); 
h2 = freqz(x2, 1, w); 
h3 = freqz(y,1,w); 

subplot(3,1,1)
plot(w/pi,abs(h1));grid
title('Magnitude Spectrum of First Sequence')
xlabel('\omega /\pi');
ylabel('Amplitude');
subplot(3,1,2)
plot(w/pi,abs(h2));grid
title('Magnitude Spectrum of Second Sequence')
xlabel('\omega /\pi');
ylabel('Amplitude');
subplot(3,1,3)
plot(w/pi,abs(h3));grid
title('Magnitude Spectrum of Product Sequence')
xlabel('\omega /\pi');
ylabel('Amplitude');

%%
% Based on the graph presented above, the DFT is modulated when you
% multiply the two sequence in the time domain as performed on the code
% above, where y is the product of x1 and x2 in the time domain. Hence, 
% the modulation property of the DFT says that the function is modulated 
% by another function if they are multiplied in the time domain.

%% F. Time-reversal property
clf;
w = -pi:2*pi/255:pi;
num = [1 2 3 4];
L = length(num)-1;
h1 = freqz(num, 1, w); 
h2 = freqz(fliplr(num), 1, w);
h3 = exp(w*L*i).*h2;
subplot(2,2,1)
plot(w/pi,abs(h1));grid
title('Magnitude Spectrum of Original Sequence')
xlabel('\omega /\pi');
ylabel('Amplitude');
subplot(2,2,2)
plot(w/pi,abs(h3));grid
title('Magnitude Spectrum of Time-Reversed Sequence')
xlabel('\omega /\pi');
ylabel('Amplitude');
subplot(2,2,3)
plot(w/pi,angle(h1));grid
title('Phase Spectrum of Original Sequence')
xlabel('\omega /\pi');
ylabel('Phase in radians');
subplot(2,2,4)
plot(w/pi,angle(h3));grid
title('Phase Spectrum of Time-Reversed Sequence')
xlabel('\omega /\pi');
ylabel('Phase in radians');

%%
% *As presented on the generated plot, the direction of phase spectrum is
% reversed as compared to the original DFT sequence. This time reversal is
% done by forming a new sequence using the function fliplr which reverse
% the samples contained on the original ramp sequence. The time reversal
% operation are as follows. First we call the freqz() function to set the
% h2 equal to DFT of the new sequence obtained from using fliplr() on the
% original ramp seuqnce. Then, h3 is set equal to the DFT of the time
% reversed ramp by multiplying h2 with a linear phase term to implement the
% shifting in the time domain.*
%

%% G.1 Circular-shifting property
% 
clf;
x = [0 2 4 6 8 10 12 14 16];
N = length(x)-1; 
n = 0:N;
y = circshift(x,5); %circshift is defined by the matlab function circshift in matlab2020a.
XF = fft(x);
YF = fft(y);

subplot(2,2,1);
stem(n,abs(XF));
grid;
title('Magnitude of DFT of Original Sequence');
xlabel('Frequency index k');
ylabel('|X[k]|');
subplot(2,2,2);
stem(n,abs(YF));
grid;
title('Magnitude of DFT of Circularly Shifted Sequence');
xlabel('Frequency index k');
ylabel('|Y[k]|');
subplot(2,2,3);
stem(n,angle(XF));grid;
title('Phase of DFT of Original Sequence');
xlabel('Frequency index k');
ylabel('arg(X[k])');
subplot(2,2,4);
stem(n,angle(YF));grid;
title('Phase of DFT of Circularly Shifted Sequence');
xlabel('Frequency index k');
ylabel('arg(Y[k])');
%%
% The above code demonstrate the circshift function in matlab for DFT.
% Circshift function operates as follows: the input sequence *x* is
% circularly shifted by *M* position. If *M > 0*, then circshift remove the
% leftmost M elements from the input function x and appends them on the
% right side of the remaining input sequence to obtain the circularly
% shifted sequence. In case when *M <0*, circshift first complements *M* by
% the length of input sequence *x*.

%%
% On the given example, we can see on its generated plots that the length
% of sequence *x* is *n=8* and the time shift is shifted to the left by five
% samples as shown in the phase term expression below.
% $W[N][kn0]=W[N][-k5]=exp(jk10pi/8)=exp(jk5pi/4)$ which increases the
% slope of phase.

%% G.2 Circular-convolution property

g1 = [1 2 3 4 5]; 
g2 = [2 2 0 1 1];
g1e = [g1 zeros(1,length(g2)-1)];
g2e = [g2 zeros(1,length(g1)-1)];
ylin = cconv(g1e,g2e);
disp('Linear convolution via circular convolution = ');
disp(ylin(1:9));
y = conv(g1, g2);
disp('Direct linear convolution = ');
disp(y);

%%
% The above code demonstrate the circular convolution property of DFT.
% circonv or cconv is a function that requires two inputs *g1e and g2e* of 
% equal length. Then, we let *g2e* be the infinite length periodic
% extension of *g2e*. Next, the elements from *1 to L* of the output vector
% *ylin* are obtained by taking the inner product between *g1e* and a
% length vector *L* which is obtained by circularly shifting to the right
% the time reversed vector *g2etr*. The output sample *y[n]* with *n*
% greater than or equal to *1* and less than or equal to *L*, the amounth
% of circular shift is *n-1* position. Based from the generated results,
% zero padding onto the match length is made possible for the
% implementation of linear convolution using circular convolution.
