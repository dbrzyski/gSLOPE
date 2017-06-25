% To run this example please
% -----
% 1. make sure that suitable SLOPE software interface C code was compiled 
%    prior to use. To do this one should to localize directory 'SLOPE_code'
%    and  run the 'makemex' script from the Matlab command line;
% 2. provide the path to the directory "gSLOPE_matlab" in line 14
% 3. provide the path to the directory "gSLOPE_code" in line 15
% 4. (optional) you can provide the path in line 16 for the directory where
%    results should be stored and uncomment lines from 118 to 122 to save 
%    the results of this simulation.
% -----

%% Addpath, file names
addpath('D:\projects\gSLOPE\matlab') %path to the directory "gSLOPE_matlab"
addpath('D:\projects\gSLOPE\matlab\SLOPE_code') %path to the directory "gSLOPE_code"
%savepath = ''; % path for results

%% Setting
groupsSizes   = [3, 4, 5, 6, 7]; % considered sizes of groups
groupsNumbers = [200, 200, 200, 200, 200]; % number of groups of each size
n             = 2000; % number of rows in design matrix
q             = 0.05; % target gFDR level
K             = [3, 6, 9, 12, 15]; % set of true support sizes
iter          = 100; % number of iteration for each k from K
eps           = 1e-6; % estimated magnitudes smaller than eps are treated as zeros

%% Objects
m             = sum(groupsNumbers);
p             = groupsNumbers*groupsSizes';
Lgths         = [];
I             = [];
for hh = 1:length(groupsNumbers)
    Lgths = [Lgths,repmat(groupsSizes(hh), 1, groupsNumbers(hh))]; %#ok
end

for gg = 1:m
    I = [I,  gg*ones(1, Lgths(gg))]; %#ok
end
l_min          = min(Lgths);
l_max          = max(Lgths);
l_mean         = round(mean(Lgths));
W              = sqrt(Lgths);

gFDR       = zeros(1,length(K));
POWER      = zeros(1,length(K));
stdgFDR    = zeros(1,length(K));
stdPOWER   = zeros(1,length(K));

%% Signal Strength
SigStrg = sqrt( 4*log(m)./(1-m.^(-2./Lgths)) - Lgths );
SigStrg = mean(SigStrg);

%% Options for optimization software
options            = struct();
options.verbosity  = 0;
options.tolInfeas  = 1e-7;
options.tolRelGap  = 1e-7;

%% Lambda and Lambda with correction
lambda = LambdaGENERAL(n, q, Lgths, W);

%% Main loops
for jj = 1:length(K)
    fprintf('K:= %d\n', K(jj))
    k = K(jj);
    power     = zeros(1, iter);
    gfdr      = zeros(1, iter);
    stdpower  = zeros(1, iter);
    stdgfdr   = zeros(1, iter);
    for ii = 1: iter
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
        [~,gSL]            = iGroupSLOPE(X, y, lambda, I, W, options);
        gSL(abs(gSL)<eps)  = 0;
        
        % getting FDR and POWER estimates
        discaveries        = find(abs(gSL)>0);
        true_disc          = intersect(discaveries, group_supp); 
        numb_disc          = length(discaveries);
        false_disc         = setdiff(discaveries, group_supp);       
        numb_false_disc    = length( setdiff(discaveries, group_supp) );
        gfdr(ii)           = numb_false_disc/max(1, numb_disc);
        power(ii)          = (numb_disc - numb_false_disc)/k;    
    end
    gFDR(jj)     = mean(gfdr);
    POWER(jj)    = mean(power);
    stdgFDR(jj)   = std(gfdr);
    stdPOWER(jj) = std(power);
end

 %% saving

%  cd(savepath)
%  save gFDR.txt -ascii gFDR;
%  save POWER.txt -ascii POWER;
%  save stdgFDR.txt -ascii stdgFDR;
%  save stdPOWER.txt -ascii stdPOWER;


