%% EE 274 Digital Signal Processing 1 Lab Activity 1
% Name: Paul Vincent S. Nonat
%

%% A. Signal Generation
% In this exercise, you will demonstrate your coding skills in MATLAB by
% generating the following signals:
%
% 
% # [y,n] = impseq(n0,a,b) 
% # [y,n] = stepseq(n0,a,b)
% # [y,n] = sigadd(x1,n1,x2,n2)
% # [y,n] = sigmult(x1,n1,x2,n2)
% # [y,n] = sigshift(x1,n1,n0)
% # [y,n] = sigfold(x1,n1)
% # [xe,xo,n] = evenodd(x1,n1)

%% 
% *[y,n] = impseq(n0,a,b)*
function [x,n]=impseq(n0,n1,n2)
% Generates x(n)= delta(n-n0); n1<=n<=n2
%--------------
%[x,n]=impseq(n0,n1,n2)
n=[n1:n2]; x=[(n-n0)==0]; 
end

%% 
% *[y,n] = stepseq(n0,a,b)*
function [x,n]=stepseq(n0,n1,n2)
% Generates x(n)= u(n-n0); n1<=n<=n2
%[x,n]=stepseq(n0,n1,n2)
n=[n1:n2]; x=[(n-n0)>=0]; 
end

%%
% *[y,n] = sigadd(x1,n1,x2,n2)*
function [y,n]=sigadd(x1,n1,x2,n2)
n=min(min(n1),min(n2)): max(max(n1),max(n2));
y1=zeros(1,length(n)); y2=y1;
y1(find((n>=min(n1))&(n<=max(n1))==1))=x1;
y2(find((n>=min(n2))&(n<=max(n2))==1))=x2;
y=y1+y2; 
end

%%
% *[y,n] = sigmult(x1,n1,x2,n2)*
function [y,n] = sigmult(x1,n1,x2,n2)
n=min(min(n1),min(n2)): max(max(n1),max(n2));
y1=zeros(1,length(n)); y2=y1;
y1(find((n>=min(n1))&(n<=max(n1))==1))=x1;
y2(find((n>=min(n2))&(n<=max(n2))==1))=x2;
y=y1*y2; 
end

%%
% *[y,n] = sigshift(x1,n1,n0)*
function [y,n] = sigshift(x,m,n0)
% implements y(n) = x(n-n0)
% -------------------------
% [y,n] = sigshift(x,m,n0)
%
n = m+n0; y = x;
end
%%
% *[y,n] = sigfold(x1,n1)*
function [y,n] = sigfold(x,n)
% implements y(n) = x(-n)
% -----------------------
% [y,n] = sigfold(x,n)
%
y = fliplr(x); n = -fliplr(n);
end
