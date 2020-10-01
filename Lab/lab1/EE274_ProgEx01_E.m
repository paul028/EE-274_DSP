%% EE 274 Digital Signal Processing 1 Lab Activity 1
% Name: Paul Vincent S. Nonat
%
%% E. Quantization
%Quantization is done by replacing each value of an analog signal $x(t)$ by
%the value of the nearest quantization level. To exemplify this oepration,
%let's simulate a unipolar ADC (Analog to Digital Converter) having the
%technical specifications: R= 10 Volts (full-scale range) and B = 3(number
%of bits).
% 
% # Write a MATLAB function y=adc_uni(x,R,B) where x and y are vectors
% containing the input signal and the quantized signal, respectively.
% # Test your function with an input ramp signal ranging from -5 to 15
% Volts (1 volt per step).
% # On the same graph, use the plot and stem functions to display the input
% signal and quantized signal respectively.

%%
% ADC_NU function test
R = 10;
B = 3;
x = -5:15;
y = adc_uni(x,R,B);
t = 0:length(x)-1;
figure(11)
plot(t,x,t,y)
plot(t,x,'g-*','LineWidth',2.2)
hold on
stem(t,y,'filled','LineWidth',2.2)
hold off
title('Ramp function unipolar quantization')
xlabel('Time in sec')
ylabel('Signal magnitude in volts')
axis([-0.1,20.1,-5.1,15.1]) 

%%
function y = adc_uni(x, R, B)
level = [0:R/(2^B):R-R/(2^B)];
temp = [-Inf,(level(2:end)-R/(2^(B+1))),Inf];
i=1
y=(x >= temp(i)).*(x < temp(i+1)).*level(i)
for i = 2:length(level)
    y = y + (x >= temp(i)).*(x < temp(i+1)).*level(i);
end
end