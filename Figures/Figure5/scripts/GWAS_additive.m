%% Addpath, file names
addpath('path to "\gSLOPE_matlab" ')
addpath('path to "\gSLOPE_matlab\SLOPE_code" ')
savepath = ''; % path for results
path = ''; %path to data

%% Load data
cd(path);
load('X.mat')

%% Adlas starting parameters
options = struct();
options.iterations = 1000;
options.verbosity = 0;
options.tolInfeas = 1e-6;
options.tolRelGap = 1e-6;

%% Objects
[n,p] = size(X);
XZ    = ExtendedMatrixNew(X);
pp    = size(XZ,2);
q     = 0.05; 
k1    = [80, 60, 40, 20, 1];
iter  = 100;
Z     = XZ(:,(p+1):end);
I     = [1:p,1:p];
L     = 2*ones(1,p);
W     = sqrt(L);
XZ2   = XZ;

fdr_gSLOPE       = zeros(1,length(k1));
sfdr_gSLOPE      = zeros(1,length(k1));
power_gSLOPE     = zeros(1,length(k1));
spower_gSLOPE    = zeros(1,length(k1));

fdr_SLOPE_X      = zeros(1,length(k1));
sfdr_SLOPE_X     = zeros(1,length(k1));
power_SLOPE_X    = zeros(1,length(k1));
spower_SLOPE_X   = zeros(1,length(k1));

fdr_SLOPE_XZ     = zeros(1,length(k1));
sfdr_SLOPE_XZ    = zeros(1,length(k1));
power_SLOPE_XZ   = zeros(1,length(k1));
spower_SLOPE_XZ  = zeros(1,length(k1));

Partial_FDP_gS      = zeros(length(k1), iter);
Partial_POW_gS      = zeros(length(k1), iter);
Partial_FDP_S_X     = zeros(length(k1), iter);
Partial_POW_S_X     = zeros(length(k1), iter);
Partial_FDP_S_XZ    = zeros(length(k1), iter);
Partial_POW_S_XZ    = zeros(length(k1), iter);

%% Orthogonalization of X matrix
for dd = 1:p
    Idd = find(I==dd);
    AA  = XZ2(:,Idd);
    [Q, R] = qr(AA,0);
    XZ2(:,Idd) = Q;        
end

%% Lambda

% lambda SLOPE_X
critical_pvalues = (1:p)*q/p;                       % critical p-values
lambda = icdf('normal',1-critical_pvalues/2,0,1);   % critical z-values 
lambda_OL = lambda;
for i=2:p
    lambda_OL(i)=  lambda(i)* sqrt(1 + sum(lambda_OL(1:(i-1)).^2)/(n-i));
end
u=lambda_OL(2:end)-lambda_OL(1:(end-1));
idx1=find(u>0,1,'first');
if (isempty(idx1))
     idx1 = p;
end
lambda_OL(idx1:end)=lambda_OL(idx1);
lambda_X = lambda_OL;

% lambda SLOPE_XZ                                          % nominal FDR control 
critical_pvalues = (1:pp)*q/pp;                       % critical p-values
lambda = icdf('normal',1-critical_pvalues/2,0,1);   % critical z-values 
lambda_OL = lambda;
for i=2:pp
    lambda_OL(i)=  lambda(i)* sqrt(1 + sum(lambda_OL(1:(i-1)).^2)/(n-i));
end
u=lambda_OL(2:end)-lambda_OL(1:(end-1));
idx1=find(u>0,1,'first');
if (isempty(idx1))
     idx1 = pp;
end
lambda_OL(idx1:end)=lambda_OL(idx1);
lambda_XZ = lambda_OL;

% lambda gSLOPE
lambda_gS = LambdaGENERAL(n, q, L, W);

%% Loop

% beta 
for i=1:length(k1)
    fprintf('k: %d\n', i) 
    k=k1(i);    
    fdp_gS        = zeros(1,iter); 
    stat_pow_gS   = zeros(1,iter);
    fdp_S_X       = zeros(1,iter); 
    stat_pow_S_X  = zeros(1,iter);
    fdp_S_XZ      = zeros(1,iter); 
    stat_pow_S_XZ = zeros(1,iter);

for b = 1:iter  
    fprintf('iter: %d\n', b)   
    ind1 = randperm(p,k);
    ind2 = p+ind1;
    II   = [ind1;ind2];
    II2  = [1:p;(p+1):pp];
      
    beta       = zeros(2*p,1); 
    beta(ind1) = (randsample([1,-1],k,1))'*5;
    beta(ind2) = (randsample([1,-1],k,1))'.*(2*rand(k,1)+3);    
    beta_II    = sqrt(sum(beta(II).^2));
    
    y  = XZ*beta + randn(n,1); 
    m  = mean(y);
    y  = y-m;    
    An = zeros(n,1); indd=[];
    
    %gSLOPE
    [gsl,GSL] = iGroupSLOPE(XZ2, y, lambda_gS, I, W, options);
    
    % SLOPE_X
    p1=1;
    s1=0; 
    ebp=zeros(1,n);    
    while s1==0 
        sigma1   = sqrt(y'*(eye(n)-An*ebp)*y/(n-p1));
        beta_hat_X = Adlas(X, y,sigma1*lambda_X,options);
        indd1    = find(abs(beta_hat_X)>0);
        if isequal(indd1,indd)
            s1=1;  
        else
            indd=indd1; An=X(:,indd);
            ebp=(An'*An)\(An');
            p1=1+length(indd);
        end
    end
    
    An = zeros(n,1); indd=[];
    
    % SLOPE_XZ
    p1=1;
    s1=0; 
    ebp=zeros(1,n);    
    while s1==0 
        sigma1   = sqrt(y'*(eye(n)-An*ebp)*y/(n-p1));
        beta_hat_XZ = Adlas(XZ,y,sigma1*lambda_XZ,options);
        indd1    = find(abs(beta_hat_XZ)>0);
        if isequal(indd1,indd)
            s1=1;  
        else
            indd=indd1; An=XZ(:,indd);
            ebp=(An'*An)\(An');
            p1=1+length(indd);
        end
    end
    beta_hat_II = sqrt(sum(beta_hat_XZ(II).^2));
    beta_hat_II2 = sqrt(sum(beta_hat_XZ(II2).^2));
    
    % Results for SLOPE_X
    SLOPE_X_supp        = find(abs(beta_hat_X)>0);
    n_discoveries       = length(SLOPE_X_supp);
    n_false_discoveries = length(setdiff(SLOPE_X_supp, ind1));
    n_true_discoveries  = n_discoveries - n_false_discoveries;     
    fdp_S_X(b)          = n_false_discoveries/max(n_discoveries,1);
    stat_pow_S_X(b)     = n_true_discoveries/k;

    % Results for SLOPE_XZ
    support             = (abs(beta_hat_II2)>0);
    n_discoveries       = sum(support>0);
    n_true_discoveries  = sum(abs(beta_hat_II)>0);
    n_false_discoveries = n_discoveries - n_true_discoveries;    
    fdp_S_XZ(b)         = n_false_discoveries/max(n_discoveries,1);
    stat_pow_S_XZ(b)    = n_true_discoveries/k;
    
    % Results for gSLOPE
    support                = (abs(GSL)>0);
    n_discoveries          = sum(support>0);
    n_true_discoveries     = sum(support(ind1) > 0);
    n_false_discoveries    = n_discoveries - n_true_discoveries;    
    fdp_gS(b)              = n_false_discoveries/max(n_discoveries,1);
    stat_pow_gS(b)         = n_true_discoveries/k;
    
    % Partial results
    Partial_FDP_gS(i,b)   = fdp_gS(b);
    Partial_POW_gS(i,b)   = stat_pow_gS(b);
    Partial_FDP_S_X(i,b)  = fdp_S_X(b);
    Partial_POW_S_X(i,b)  = stat_pow_S_X(b);
    Partial_FDP_S_XZ(i,b) = fdp_S_XZ(b);
    Partial_POW_S_XZ(i,b) = stat_pow_S_XZ(b);
    
    % Saving partial results
    cd(savepath)
    save Partial_FDP_gS.txt -ascii Partial_FDP_gS
    save Partial_POW_gS.txt -ascii Partial_POW_gS
    save Partial_FDP_S_X.txt -ascii Partial_FDP_S_X
    save Partial_POW_S_X.txt -ascii Partial_POW_S_X
    save Partial_FDP_S_XZ.txt -ascii Partial_FDP_S_XZ
    save Partial_POW_S_XZ.txt -ascii Partial_POW_S_XZ

end

fdr_gSLOPE(i)       = mean(fdp_gS);
sfdr_gSLOPE(i)      = std(fdp_gS);
power_gSLOPE(i)     = mean(stat_pow_gS);
spower_gSLOPE(i)    = std(stat_pow_gS);

fdr_SLOPE_X(i)      = mean(fdp_S_X);
sfdr_SLOPE_X(i)     = std(fdp_S_X);
power_SLOPE_X(i)    = mean(stat_pow_S_X);
spower_SLOPE_X(i)   = std(stat_pow_S_X);

fdr_SLOPE_XZ(i)     = mean(fdp_S_XZ);
sfdr_SLOPE_XZ(i)    = std(fdp_S_XZ);
power_SLOPE_XZ(i)   = mean(stat_pow_S_XZ);
spower_SLOPE_XZ(i)  = std(stat_pow_S_XZ);


cd(savepath)

save fdr_gSLOPE.txt -ascii fdr_gSLOPE
save sfdr_gSLOPE.txt -ascii sfdr_gSLOPE
save power_gSLOPE.txt -ascii power_gSLOPE
save spower_gSLOPE.txt -ascii spower_gSLOPE

save fdr_SLOPE_X.txt -ascii fdr_SLOPE_X
save sfdr_SLOPE_X.txt -ascii sfdr_SLOPE_X
save power_SLOPE_X.txt -ascii power_SLOPE_X
save spower_SLOPE_X.txt -ascii spower_SLOPE_X

save fdr_SLOPE_XZ.txt -ascii fdr_SLOPE_XZ
save sfdr_SLOPE_XZ.txt -ascii sfdr_SLOPE_XZ
save power_SLOPE_XZ.txt -ascii power_SLOPE_XZ
save spower_SLOPE_XZ.txt -ascii spower_SLOPE_XZ


end
