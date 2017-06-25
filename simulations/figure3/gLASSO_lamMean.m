% This scripts estimates FDR and POWER for uncorected (O) and corrected (S)
% version of group SLOPE algorithm for case when groups sizes are equal

%% Addpath, file names
addpath('') %path to the directory "gSLOPE_matlab"
addpath('') %path to the directory "gSLOPE_code"
savepath = ''; % path for results

%% Objects
m         = 1000;
Lgths     = importdata('Lgths.txt');
I         = importdata('I.txt');
p         = length(I);
n         = 5000;
q1        = 0.1;
l_min     = min(Lgths);
l_max     = max(Lgths);
l_mean    = round(mean(Lgths));
W         = sqrt(Lgths);
K         = [3, 10, 20, 30, 40, 50, 60]; % set of true support sizes
iter      = 200; % number of iteration for each k from K
eps       = 1e-6; % estimated magnitudes smaller than eps are treated as zeros
k_max     = 50;
FDR_01        = zeros(1,length(K));
POWER_01      = zeros(1,length(K));
stdFDR_01     = zeros(1,length(K));
stdPOWER_01   = zeros(1,length(K));

%% Signal Strength
SigStrg = sqrt( 4*log(m)./(1-m.^(-2./Lgths)) - Lgths );
SigStrg = mean(SigStrg);
%% Starting parameters

%gSLOPE algorithm parameters
options            = struct();
options.verbosity  = 0;
options.tolInfeas  = 1e-7;
options.tolRelGap  = 1e-7;

%% Lambda and Lambda with correction
lambda1 = LambdaGENERAL(n, q1, Lgths, W);

%% Main loops
for jj = 1:length(K)
    fprintf('K:= %d\n', K(jj))
    k = K(jj);
    power_01    = zeros(1, iter);
    fdr_01      = zeros(1, iter);
    stdpower_01 = zeros(1, iter);
    stdfdr_01   = zeros(1, iter);
    meanLGTHs   = [];
    parfor ii = 1: iter
        fprintf('iter:= %d\n', ii)
        %% design matrix generating
        X = randn(n,p);
        X = bsxfun(@minus, X, mean(X)); % columns in X are centered
        X = bsxfun( @rdivide, X, sqrt(sum(bsxfun( @times, X, X)))); % columns in X are normalized
        for dd = 1:m
            Idd = find(I==dd);
            A   = X(:,Idd);
            [Q, R] = qr(A,0);
            X(:,Idd) = Q;        
        end % groups in X are orthogonalized

        % true beta generating
        beta       = zeros(p,1);
        group_supp = randsample(m,k);
        z          = randn(n,1);   
        for gg = 1:k
            ggroup       = find(I == group_supp(gg));
            ggsize       = length(ggroup);
            ggsignals    = abs(rand(ggsize,1)+0.1);
            beta(ggroup) = ggsignals*SigStrg/norm(ggsignals,2);
        end
        
        % observations generating
        y = X*beta + randn(n,1);
        
        % getting solutions
        [~, SL_01]   = iGroupSLOPE(X, y, ones(1,m)*lambda1(1), I, W, options);
        SL_01(abs(SL_01)<eps)   = 0;
        meanLGTHs = [meanLGTHs, Lgths(SL_01>0)];
        
        % getting FDR and POWER estimates
        discaveries_01        = find(abs(SL_01)>0);
        true_disc_01          = intersect(discaveries_01, group_supp);
        numb_disc_01          = length(discaveries_01);
        false_disc_01         = setdiff(discaveries_01, group_supp);      
        numb_false_disc_01    = length( setdiff(discaveries_01, group_supp) );
        fdr_01(ii)            = numb_false_disc_01/max(1,numb_disc_01);
        power_01(ii)          = (numb_disc_01-numb_false_disc_01)/k;   
    end
    FDR_01(jj)      = mean(fdr_01);
    POWER_01(jj)    = mean(power_01);
    stdFDR_01(jj)   = std(fdr_01);
    stdPOWER_01(jj) = std(power_01);
    meanLGTH(jj)    = mean(meanLGTHs);


 %% saving
 cd(savepath)
%    
 save gFDR_gLASSO_LM_01.txt -ascii FDR_01;
 save POWER_gLASSO_LM_01.txt -ascii POWER_01;
 save std_gFDR_gLASSO_LM_01.txt -ascii stdFDR_01;
 save std_POWER_gLASSO_LM_01.txt -ascii stdPOWER_01;
 save meanLGTH.txt -ascii meanLGTH;

 end
