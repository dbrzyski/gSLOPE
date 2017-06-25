% This scripts estimates FDR and POWER for uncorected (O) and corrected (S)
% version of group SLOPE algorithm for case when groups sizes are equal

%% Addpath, file names
addpath('') %path to the directory "gSLOPE_matlab"
addpath('') %path to the directory "gSLOPE_code"
savepath = ''; % path for results

%% Groups generating
I = [];
groups_sizes = 3:7;
groups_times = 200*ones(1,5);
groups_idx   = 1;
for ii=1:length(groups_sizes)
    current_size = groups_sizes(ii);
    for jj = 1:groups_times(ii)
        I = [I, groups_idx*ones(1, current_size)];
        groups_idx = groups_idx + 1;
    end
end

a = cumsum(groups_times);
b = zeros(1,a(end));
b(a - groups_times + 1) = 1;
Lgths = groups_sizes(cumsum(b));

%% Objects
p         = length(I);
m         = max(I);
n         = 5000;
q1        = 0.05;
q2        = 0.1;
l_min     = min(Lgths);
l_max     = max(Lgths);
l_mean    = round(mean(Lgths));
W         = sqrt(Lgths);
K         = [3, 25, 50, 75, 100, 125, 150, 175, 200, 225, 250]; % set of true support sizes
iter      = 300; % number of iteration for each k from K
eps       = 1e-8; % estimated magnitudes smaller than eps are treated as zeros
k_max     = 50;
FDR_005        = zeros(1,length(K));
POWER_005      = zeros(1,length(K));
FDR_01         = zeros(1,length(K));
POWER_01       = zeros(1,length(K));
stdFDR_005     = zeros(1,length(K));
stdPOWER_005   = zeros(1,length(K));
stdFDR_01      = zeros(1,length(K));
stdPOWER_01    = zeros(1,length(K));
PROP_TRUE_005  = zeros(length(K), length(groups_sizes));
PROP_TRUE_01   = zeros(length(K), length(groups_sizes));
PROP_FALSE_005 = zeros(length(K), length(groups_sizes));
PROP_FALSE_01  = zeros(length(K), length(groups_sizes));

%% Signal Strength
SigStrg = sqrt( 4*log(m)./(1-m.^(-2./Lgths)) - Lgths );
SigStrg = mean(SigStrg);
aa = (m/sum(sqrt(Lgths)))*SigStrg;
%BasicStrg = 1;
%% Starting parameters

%gSLOPE algorithm parameters
options            = struct();
options.verbosity  = 0;
options.tolInfeas  = 1e-7;
options.tolRelGap  = 1e-7;

%% Lambda and Lambda with correction
lambda1 = LambdaMAX(n, q1, Lgths, W, 'orthogonal');
lambda2 = LambdaMAX(n, q2, Lgths, W, 'orthogonal');

%% Main loops
for jj = 1:length(K)
    fprintf('K:= %d\n', K(jj))
    k = K(jj);
    power_005  = zeros(1, iter);
    fdr_005    = zeros(1, iter);
    power_01   = zeros(1, iter);
    fdr_01     = zeros(1, iter);
    stdpower_005 = zeros(1, iter);
    stdfdr_005   = zeros(1, iter);
    stdpower_01  = zeros(1, iter);
    stdfdr_01    = zeros(1, iter);
    propTRUE_005   = zeros(1,length(groups_sizes));
    propTRUE_01    = zeros(1,length(groups_sizes));
    propFALSE_005  = zeros(1,length(groups_sizes));
    propFALSE_01   = zeros(1,length(groups_sizes));
    for ii = 1: iter
        fprintf('iter:= %d\n', ii)
        
        % design matrix generating
        X = eye(p,p);
        
        % true beta generating
        beta       = zeros(p,1);
        group_supp = randsample(m,k);
        z          = randn(n,1);   
        for gg = 1:k
            ggroup       = find(I == group_supp(gg));
            ggsize       = length(ggroup);
            ggsignals    = abs(rand(ggsize,1)+0.1);
            beta(ggroup) = ggsignals*aa*sqrt(ggsize)/norm(ggsignals,2);
        end
        
        % observations generating
        y = X*beta + randn(n,1);
        
        % getting solutions
        [~,SL_005]   = GroupSLOPE(X, y, lambda1, I, W, options);
        [~,SL_01]    = GroupSLOPE(X, y, lambda2, I, W, options);
        SL_005(abs(SL_005)<1e-8)   = 0;
        SL_01(abs(SL_01)<1e-8)     = 0;
        
        % getting FDR and POWER estimates
        discaveries_005        = find(abs(SL_005)>0);
        discaveries_01         = find(abs(SL_01)>0);
        true_disc_005          = intersect(discaveries_005, group_supp);
        true_disc_01           = intersect(discaveries_01, group_supp); 
        numb_disc_005          = length(discaveries_005);
        numb_disc_01           = length(discaveries_01);
        false_disc_005         = setdiff(discaveries_005, group_supp);
        false_disc_01          = setdiff(discaveries_01, group_supp);       
        numb_false_disc_005    = length( setdiff(discaveries_005, group_supp) );
        numb_false_disc_01     = length( setdiff(discaveries_01, group_supp) );
        fdr_005(ii)            = numb_false_disc_005/max(1,numb_disc_005);
        fdr_01(ii)             = numb_false_disc_01/max(1,numb_disc_01);
        power_005(ii)          = (numb_disc_005-numb_false_disc_005)/k;
        power_01(ii)           = (numb_disc_01-numb_false_disc_01)/k;
        TrueChosenSizes_005    = Lgths(true_disc_005);
        TrueChosenSizes_01     = Lgths(true_disc_01);
        FalseChosenSizes_005   = Lgths(false_disc_005);
        FalseChosenSizes_01    = Lgths(false_disc_01);
        for lll = 1: length(groups_sizes)            
            propTRUE_005(lll) = propTRUE_005(lll) + length(find(TrueChosenSizes_005==groups_sizes(lll)));
            propTRUE_01(lll)  = propTRUE_01(lll) + length(find(TrueChosenSizes_01==groups_sizes(lll)));
            propFALSE_005(lll) = propFALSE_005(lll) + length(find(FalseChosenSizes_005==groups_sizes(lll)));
            propFALSE_01(lll) = propFALSE_01(lll) + length(find(FalseChosenSizes_01==groups_sizes(lll)));
        end       
    end
    FDR_005(jj)      = mean(fdr_005);
    FDR_01(jj)       = mean(fdr_01);
    POWER_005(jj)    = mean(power_005);
    POWER_01(jj)     = mean(power_01);
    stdFDR_005(jj)   = std(fdr_005);
    stdFDR_01(jj)    = std(fdr_01);
    stdPOWER_005(jj) = std(power_005);
    stdPOWER_01(jj)  = std(power_01);
    PROP_TRUE_005(jj,:)  = propTRUE_005./sum(propTRUE_005);
    PROP_TRUE_01(jj,:)   = propTRUE_01./sum(propTRUE_01);
    PROP_FALSE_005(jj,:) = propFALSE_005./sum(propFALSE_005);
    PROP_FALSE_01(jj,:)  = propFALSE_01./sum(propFALSE_01);
end

 %% saving
 cd(savepath)
% 
 save FDR_005_MAX_WEAK.txt -ascii FDR_005;
 save FDR_01_MAX_WEAK.txt -ascii FDR_01;
 save POWER_005_MAX_WEAK.txt -ascii POWER_005;
 save POWER_01_MAX_WEAK.txt -ascii POWER_01;
 
 save stdFDR_005_MAX_WEAK.txt -ascii stdFDR_005;
 save stdFDR_01_MAX_WEAK.txt -ascii stdFDR_01;
 save stdPOWER_005_MAX_WEAK.txt -ascii stdPOWER_005;
 save stdPOWER_01_MAX_WEAK.txt -ascii stdPOWER_01;
 
 save PROP_TRUE_005_MAX_WEAK.txt -ascii PROP_TRUE_005;
 save PROP_TRUE_01_MAX_WEAK.txt -ascii PROP_TRUE_01;
 save PROP_FALSE_005_MAX_WEAK.txt -ascii PROP_FALSE_005;
 save PROP_FALSE_01_MAX_WEAK.txt -ascii PROP_FALSE_01;
 
