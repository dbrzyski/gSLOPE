% The iterative version of group SLOPE 
function [gSLOPE, xg, SIGMA] = iGroupSLOPE(X, y, lambda, I, W, options) 

%% objects
epsilon = 1e-8;
[n,p]   = size(X);
indd    = [];
p1      = 1;
s1      = 0;
Xn      = zeros(n,1);
ebp     = zeros(1,n);

%% data structure
m     = max(I);
Lgths = zeros(m,1);
maxlength = 1;
for ii=1:m
    Lgths(ii) = sum(I==ii);
    if Lgths(ii)>maxlength
        maxlength = Lgths(ii);
    end
end
I2 = (p+1)*ones(maxlength,m);
for jj=1:m
    columnn = find(I==jj);
    columnn = columnn(:);
    I2(1:length(columnn),jj) = columnn;
end

%% Loop
while s1==0         
    sigma1  = sqrt(y'*(eye(n)-Xn*ebp)*y/(n-p1));
    [gSLOPE, xg] = GroupSLOPE(X, y, sigma1*lambda, I, W, options);
    indd1  = find(abs(xg)>epsilon);
    if length(indd)>n
        disp('solution not sparse')
    end
    if isequal(indd1,indd)
        s1 = 1;  
    else
        indd   = indd1;
        g_idxs = I2(:,indd);
        g_idxs = g_idxs(:);
        g_idxs = g_idxs(g_idxs<p+1);
        Xn     = X(:,g_idxs); 
        ebp    = (Xn'*Xn)\Xn';
        p1     = 1+size(Xn,2);
    end;
end; 
SIGMA = sigma1;
end