function y=dt_2(x)
y= zeros(1,length(x)); % zero initialization

for n=1:length(x);
    if n==1
        y(n) = x(n);
    elseif n==2
        y(n) = x(n) - 2*y(n-1);
    else
        y(n) = x(n) - 2*y(n-1) - 2*y(n-2)
    end
end
end
