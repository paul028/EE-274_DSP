function y=dt_3(x)
y= zeros(1,length(x));
b = [1.5 0.5];  %input coefficients n , n-1
a = [1 2 2]; % output coefficients n, n-1, n-2
h = impz(b,a); %impulse response
y = conv(h,x);
end