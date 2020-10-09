function y=dt_4(x,L)
y= zeros(1,length(x));
b = [1];  % input coefficient n
a_l = zeros(1,L-1);  %
a = [1 a_l -0.5 -0.5] % output coefficient n, n-L, n-L-1
h = impz(b,a); %impulse response
y = conv(h,x);