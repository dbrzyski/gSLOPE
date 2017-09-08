
function [clumps, clumps_n, clms_info, repr]  = geneSLOPE_XZ(Pheno, rho, q)

%This function calculates the geneSLOPE estimate for genetic matrix [X Z] 
%and phenotype data (model with dominance is assumed)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%INPUTS:
% Pheno - number of phenotype : 1 = TG, 2 = HDL, 3 = LDL, 4= TotalCHOL
% rho   - paramterer used in clustering
% q     - assumed FDR level

%% Addpath, file names
addpath('path to "\gSLOPE_matlab" ')
addpath('path to "\gSLOPE_matlab\SLOPE_code" ')
path = ''; %path to data

%% Loading data
cd(path)
Chr_sizes = importdata(''); % numbers of SNPs in each chromosome
Chr_sizes2 = [0; cumsum(Chr_sizes(1:(end-1)))];

%% Phenotypes
cd(path)
Phen        = importdata(''); %data with phenotype
y           = Phen(:,Pheno);
nanidx      = find(isnan(y));
y(nanidx,:) = [];
n           = length(y);

%% Objects
cut_level     = 0.05;
epsilon       = 1e-10;
p_gw          = sum(Chr_sizes);
Xcut          = zeros(n, 50000);
Xcut_info     = zeros(2, 50000);
Idx_list      = zeros(1,5000);
ind_b         = 1;
ind_e         = 0;

%% Adlas starting parameters
options           = struct();
options.verbosity = 0;
options.tolInfeas = 1e-7;
options.tolRelGap = 1e-7;

%% Loop (among all chromosomes)
cd(path)
for cc = 1:23
    pname = strcat( 'X_chr',int2str( cc ),'.mat' ); %data were saved separately for each chromosome
    X     = importdata(pname);
    X(abs(X)<1e-10) = NaN;
    X(nanidx,:) = [];   
    colMean     = nanmean(X);
    [~,col]   = find(isnan(X));
    X(isnan(X)) = colMean(col);    
    X           = bsxfun(@minus, X, mean(X));
    X(abs(X)<1e-10) = 0;
    X           = bsxfun( @rdivide, X, sqrt(sum(bsxfun( @times, X, X))));
    n           = size(X,1);    

    % ANOVA tests
    px = size(X,2);
    pv = zeros(1,px);
    for iii = 1:px
        a   = X(:,iii);
        a_rr = a(a~=0);
        y_rr = y(a~=0);
        pv(iii) = anova1(y_rr, a_rr,'off');
    end
           
    % information about pruned SNPs
    Cut_idxs = find(pv < cut_level);
    ind_e    = length(Cut_idxs)+ind_e;
    Xcut_info(1, ind_b:ind_e) = cc.*ones(1,length(Cut_idxs));
    Xcut_info(2, ind_b:ind_e) = Cut_idxs;
    Xcut(:, ind_b:ind_e)      = X(: , Cut_idxs);
    Idx_list(ind_b:ind_e)     = Cut_idxs + Chr_sizes2(cc);
    ind_b    = ind_e + 1;
end
Xcut_info = Xcut_info(:,1:ind_e);
Xcut      = Xcut(:,1:ind_e);
Idx_list  = Idx_list(:,1:ind_e);
X_corr    = Xcut'*Xcut;
   
% Clumps generating
[STRUCT, r_idxs] = StructNew(Xcut, y, rho, X_corr);
Xrps             = Xcut(:,r_idxs);
p                = size(Xrps, 2 );              

% Estendent matrix XZ
XZ = ExtendentMatrixNew(Xrps);
I  = [1:p,1:p];

% Lambda
critical_pvalues = (1:p_gw)*q/p_gw;                       
lambda           = icdf('normal',1-critical_pvalues/2,0,1);
lambda_OL        = lambda;
for i=2:500
    lambda_OL(i) =  lambda(i)* sqrt(1 + sum(lambda_OL(1:(i-1)).^2)/(n-i));
end
u    = lambda_OL(2:end)-lambda_OL(1:(end-1));
idx1 = find(u>0,1,'first');
if (isempty(idx1))
     idx1 = p_lambda;
end
lambda_OL(idx1:end)  = lambda_OL(idx1);
lambda               = lambda_OL;
lambda               = lambda(1: (2*p));

% SLOPE for representatives, iterative version      
indd = [];
p1   = 1;
s1   = 0;
Xn   = zeros(n,1);
ebp  = zeros(1,n);
while s1==0         
    sigma1 = sqrt(y'*(eye(n)-Xn*ebp)*y/(n-p1));
    SLOPE  = Adlas(XZ,y,sigma1*lambda,options);
    indd1  = find(abs(SLOPE)>epsilon);
    if length(indd)>n
        disp('solution not sparse')
    end
    if isequal(indd1,indd)
        s1 = 1;  
    else indd = indd1; Xn=XZ(:,indd); 
          ebp = (Xn'*Xn)\Xn';
          p1  = 1+length(indd);
    end
end            
II2          = [1:p;(p+1):(2*p)];
SLOPE_II2    = sqrt(sum(SLOPE(II2).^2));
supp         = find(SLOPE_II2>0);

% Chosen clumps
clumps   = [];
clumps_n = [];
repr     = [];
for c=1:length(supp)
    repr     = [repr, r_idxs(supp(c)) ];
    ChGroup   = find(STRUCT == STRUCT(r_idxs(supp(c))))';
    clumps    = [clumps; ChGroup];
    clumps_n  = [clumps_n; length(ChGroup)];
end
clms_info = Xcut_info(:,clumps);
clumps    = Idx_list(clumps);
repr      = Idx_list(repr);

end