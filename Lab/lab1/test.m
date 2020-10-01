% This script creates a random signal, and then quantizes it.  The signal
% is oversampled and then decimated.
% * Oversampling is the process of taking in samples at a faster rate (in
% this script "os" times faster) than you need.
% * Decimation is the process Decimation is the process of taking only one
% of every "os" samples of an oversampled signal to get the final sampling
% rate.
%
% Processing of the oversampled signal can give some benefit, as you will
% see if you run the script.

fprintf('\n\nE71 Lab, Oversampling and Quantization\n');

b=3;                            % Number of bits.
N=100;                          % Number of samples in final signal.
maxOS = 4;                      % Max os rate = 2^maxOS.

% Choose the input type.
choice = questdlg('Choose input','Input',...
    'Sine','Sawtooth','Random','Random');

fprintf('Bits = %g, levels = %g, signal = %s.\n', b, 2^b, choice);
fprintf('SQNR measured in dB.\n\n');

%Generate set of points.
    os=2^maxOS;
    N_os=N*os;                   % Number of samples in oversampled signal
    n=0:(N_os-1);                % Index
    
    % Create the oversampled input data sequence.
    switch choice
        case 'Sine'
            xOrig=sin(2*pi*n/N_os);
        case 'Sawtooth'
             xOrig=sawtooth(2*pi*n/N_os);
        case 'Random'
             xOrig=randn(1,N_os);        % Random data
            % Smooth to begin to remove fast variations.
             xOrig=filter(ones(1,4*os)/4/os,1,xOrig);
             xOrig= xOrig/max(abs(xOrig));        % Scale to +/- 1
    end
    
% This large loop generates and analyzes data at several different
% oversampling rates (all powers of two).
for os_pow=0:maxOS,
    os=2^os_pow;                    % Oversampling rate.
    xds = downsample(xOrig,os);
   % x = filter(ones(1,4*os)/4/os,1,xOrig);
    x = filter(ones(1,4)/4,1,xds);
    x = x/max(abs(x));
    
    %Quantize the oversampled raw signal
    xq=floor((x+1)*2^(b-1));	%Signal is one of 2^b int values (0 to 2^b-1)
    xq=xq/(2^(b-1));            %Signal is from 0 to 2 (quantized)
    xq=xq-(2^(b)-1)/2^(b);      %Shift signal down (rounding)
    
    %Smooth (running average) the quantized oversampled signal
    x_qs=filter(ones(1,os)/os,1,xq);
    
    %Smooth the oversampled signal
    x_s=filter(ones(1,os)/os,1,x);
    
    %Quantize the oversampled smoothed signal
    x_sq=floor((x_s+1)*2^(b-1)); % Signal is one of 2^n int values (0 to 2^n-1)
    x_sq=x_sq/(2^(b-1));         % Signal is from 0 to 2 (quantized)
    x_sq=x_sq-(2^(b)-1)/2^(b);   % Shift signal down (rounding)
       
    xe_qs=x_s-x_qs;             % Error between smoothed and quantized/smoothed
    xe_sq=x_s-x_sq;             % Error between smoothed and smoothed/quantized
    

    subplot(3,1,1);
    stem(x,'m');
    hold on;
    stem(x_s,'b');
    hold off;
    set(gca,'YLim',[-1.1,1.1]);
    legend({'Original Signal','Smoothed Signal'});
    title('Original Signal (downsampled) and smoothed signal');
    
    subplot(3,1,2)
    stem(x_s,'b');
    hold on;
    stem(x_sq,'r');
    hold on;
    stem(xe_sq*2^(b-1),'g');
    legend('exact','quantized','error*2^{(b-1)}','Location','Southeast')
    set(gca,'YLim',[-1.1,1.1]);
    title(sprintf('Smoothed then quantized, %g bits, %g levels, os=%g ',b,2^b, os));
    hold off
    
    subplot(3,1,3)
    stem(x_s,'b');
    hold on;
    stem(x_qs,'r');
    hold on;
    stem(xe_qs*2^(b-1),'g');
    legend('exact','quantized','error*2^{(b-1)}','Location','Southeast')
    set(gca,'YLim',[-1.1,1.1]);
    title(sprintf('Quantized then smoothed, %g bits, %g levels, os=%g ',b,2^b, os));
    hold off
    
    sqnr_sq=10*log10(sum(x_s.^2)/sum(xe_sq.^2));
    sqnr_qs=10*log10(sum(x_s.^2)/sum(xe_qs.^2));
    
    fprintf('Oversampling = %g, ', os);
    fprintf('sqnr_sq = %g, sqnr_qs = %g, sqnr_qs-sqnr_sq = %g\n',...
        sqnr_sq, sqnr_qs, sqnr_qs-sqnr_sq);
    
    pause(1);
end