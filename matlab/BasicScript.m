% This scripts estimates gFDR and POWER for uncorected (O) and corrected (S)
% version of group SLOPE algorithm for case when groups sizes are equal

%% Addpath, file names
addpath('') %path to the directory "gSLOPE_matlab"
addpath('') %path to the directory "gSLOPE_code"
savepath = ''; % path for results

%% Groups generating
I = [];
groups_sizes = [5,15];
groups_times = [250,250];

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
n         = p;
q         = 0.1;
l_min     = min(Lgths);
l_max     = max(Lgths);
W         = sqrt(Lgths);
K         = [10, 20]; % set of true support sizes
iter      = 50; % number of iteration for each k from K
eps       = 1e-8; % estimated magnitudes smaller than eps are treated as zeros
FDR_O     = zeros(1,length(K));
POWER_O   = zeros(1,length(K));
FDR_S     = zeros(1,length(K));
POWER_S   = zeros(1,length(K));
stdFDR_O     = zeros(1,length(K));
stdPOWER_O   = zeros(1,length(K));
stdFDR_S     = zeros(1,length(K));
stdPOWER_S   = zeros(1,length(K));

%% Signal Strength
%l_sig = round(mean(l_min,l_max));
l_sig = l_max;
SigStrg = 1.2*sqrt( 4*log(m)/(1-m^(-2/l_sig)) - l_sig );

%% Starting parameters

%gSLOPE algorithm parameters
options            = struct();
options.verbosity  = 0;
options.tolInfeas  = 1e-7;
options.tolRelGap  = 1e-7;

%% Lambda and Lambda with correction
l_lam = l_min;
w_lam = sqrt(l_lam);
lambda1 = LambdaGENERAL(n, q, l_lam.*ones(m,1), w_lam.*ones(m,1));
lambda2 = LambdaGENERAL(n, q, Lgths, W);


%% Main loops
for jj = 1:length(K)
    fprintf('K:= %d\n', K(jj))
    k = K(jj);
    power_O = zeros(1, iter);
    fdr_O   = zeros(1, iter);
    power_S = zeros(1, iter);
    fdr_S   = zeros(1, iter);
    stdpower_O = zeros(1, iter);
    stdfdr_O   = zeros(1, iter);
    stdpower_S = zeros(1, iter);
    stdfdr_S   = zeros(1, iter);
    for ii = 1: iter
        fprintf('iter:= %d\n', ii)
        
        % design matrix generating
        X = sqrt(1/n).*randn(n,p);
        X = bsxfun(@minus, X, mean(X)); % columns in X are centered
        X = bsxfun( @rdivide, X, sqrt(sum(bsxfun( @times, X, X)))); % columns in X are normalized
        
        % true beta generating
        beta1       = zeros(p,1);
        beta2       = zeros(p,1);
        group_supp1 = randsample((1:groups_times(2)) + groups_times(1),k);
        group_supp2 = group_supp1;
        z           = randn(n,1);   
        for gg = 1:k
            ggroup       = find(I == group_supp1(gg));
            ggsize       = length(ggroup);
            ggindex = randsample(ggroup,1);
            beta1(ggindex) = SigStrg;
            
            ggroup         = find(I == group_supp2(gg));
            ggsize         = length(ggroup);
            ggindex        = randsample(ggroup,1);
            beta2(ggindex) = SigStrg;
        end
        
        % observations generating
        z  = randn(n,1);
        y1 = X*beta1 + z;
        y2 = X*beta2 + z;
        
        % getting solutions
        [~,ODS]   = GroupSLOPE(X, y1, lambda1, I, W, options);
        [~,SLOPE]   = GroupSLOPE(X, y2, lambda2, I, W, options);
        ODS(abs(ODS)<1e-8)   = 0;
        SLOPE(abs(SLOPE)<1e-8)   = 0;
        
        % getting FDR and POWER estimates
        discaveries_O     = find(abs(ODS)>0);
        discaveries_S     = find(abs(SLOPE)>0);
        numb_disc_O       = length(discaveries_O);
        numb_disc_S       = length(discaveries_S);
        false_disc_O      = length( setdiff(discaveries_O, group_supp1) );
        false_disc_S      = length( setdiff(discaveries_S, group_supp2) );
        fdr_O(ii)         = false_disc_O/max(1,numb_disc_O);
        fdr_S(ii)         = false_disc_S/max(1,numb_disc_S);
        power_O(ii)       = (numb_disc_O-false_disc_O)/k;
        power_S(ii)       = (numb_disc_S-false_disc_S)/k;
    end
    FDR_O(jj)      = mean(fdr_O);
    FDR_S(jj)      = mean(fdr_S);
    POWER_O(jj)    = mean(power_O);
    POWER_S(jj)    = mean(power_S);
    stdFDR_O(jj)      = std(fdr_O);
    stdFDR_S(jj)      = std(fdr_S);
    stdPOWER_O(jj)    = std(power_O);
    stdPOWER_S(jj)    = std(power_S);
end

%% saving
cd(savepath)
save POWER_ss_ln_nonc.txt -ascii POWER_O;
save POWER_ss_ln_c.txt -ascii POWER_S;
save FDR_ss_ln_nonc.txt -ascii FDR_O;
save FDR_ss_ln_c.txt -ascii FDR_S;

save stdPOWER_ss_ln_nonc.txt -ascii stdPOWER_O;
save stdPOWER_ss_ln_c.txt -ascii stdPOWER_S;
save stdFDR_ss_ln_nonc.txt -ascii stdFDR_O;
save stdFDR_ss_ln_c.txt -ascii stdFDR_S;


