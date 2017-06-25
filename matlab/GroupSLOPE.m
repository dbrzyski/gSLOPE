function [x, xg, info] = GroupSLOPE(A, b, lambda, I, W, options)
% MeanGroupAdlas solves group SLOPE problem
%
%       Minimize 1/2*||Ax-b||_2^2 + J_lambda(W*||b_I||_2),
%
% where W is diagonal matrix with weights of groups on diagonal and I is
% vector containing the structure of data in form: group jth is given by the
% set of indices I_j = find(I==j). By ||b_I||_2 we denote vector of l_2 norms 
% at groups level, i.e. ( ||b_I||_2 )_j: = || b(I_j) ||_2. 
%
% The options parameter is a structure with the following optional fields
% with [default value]: 
%
%    .iterations    Maximum number of iterations                  [10,000]
%    .verbosity     0 = nothing, 1 = major, 2 = every                  [1]
%    .fid           File identifier for output                [1 = stdout]
%    .optimIter     Iterations between optimality-condition checks     [1]
%    .gradIter      Iterations between full gradient computations     [20]
%    .tolInfeas     Maximum allowed dual infeasibility              [1e-6]
%    .tolRelGap     Stopping criterion for relative primal-dual gap [1e-6]
%    .xInit         Initial value of x                        [zeros(n,1)]
%
% The info output structure contains the following fields
%
%    .runtime       Runtime
%    .Aprods        Number of products with A
%    .ATprods       Number of products with A^T
%    .objPrimal     Primal objective
%    .objDual       Dual objective (possibly for infeasible dual point)
%    .infeas        Dual infeasibility
%    .status        Status: 1 = optimal, 2 = iterations
%

% -------------------------------------------------------------
% Start timer
% -------------------------------------------------------------
t0 = tic();

% Get problem dimension
[n,p]   = size(A);

% -------------------------------------------------------------
% Parse parameters
% -------------------------------------------------------------
if (nargin <  6), options = struct(); end;
if (nargin == 4), m = max(I); W = ones(1,m); end;
if (nargin <  4), I = 1:p; m = p; W = ones(1,m); end;

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
W = W(:);
W = (W.^(-1));
W = W(I);

iterations = getDefaultField(options,'iterations',10000);
verbosity  = getDefaultField(options,'verbosity',0);
fid        = getDefaultField(options,'fid',1);
optimIter  = getDefaultField(options,'optimIter',1);
gradIter   = getDefaultField(options,'gradIter',20);
tolInfeas  = getDefaultField(options,'tolInfeas',1e-6);
tolRelGap  = getDefaultField(options,'tolRelGap',1e-6);
xInit      = getDefaultField(options,'xInit',[]);

% Ensure that lambda is non-increasing
if ((length(lambda) > 1) && any(lambda(2:end) > lambda(1:end-1)))
   error('Lambda must be non-increasing.');
end
if (lambda(end) < 0)
   error('Lambda must be nonnegative');
elseif (lambda(1) == 0)
   error('Lambda must have at least one nonnegative entry.');
end

% Ensure that K have the same size as lambda
if length(lambda)~=m
    error( 'number of groups must be equal to dimension of lambda');
end

% Ensure that K discribes groups division properly
if length(I)~=p
    error('length of I has to be equal to p');
end

% -------------------------------------------------------------
% Initialize
% -------------------------------------------------------------

% Get initial lower bound on the Lipschitz constant
s = RandStream('mt19937ar','Seed',0);
x = randn(s,p,1); x = x / norm(x,2);
x = W.*(A'*(A*(W.*x)));            %x = A'*(A*x); 
L = norm(x,2);

% Constants for exit status
STATUS_RUNNING    = 0;
STATUS_OPTIMAL    = 1;
STATUS_ITERATIONS = 2;
STATUS_MSG = {'Optimal','Iteration limit reached'};

% Initialize parameters and iterates
if (isempty(xInit)), xInit = zeros(p,1); end;

t       = 1;
eta     = 2;
lambda  = lambda(:);
b       = b(:);
x       = xInit;
y       = x;
Ax      = A*(W.*x);      %Ax      = A*x;
%fPrev   = Inf;
iter    = 0;
status  = STATUS_RUNNING;
Aprods  = 2;
ATprods = 1;

% Deal with Lasso case
modeLasso = (numel(lambda) == 1);
if (modeLasso)
   proxFunction = @(v1,v2) proxL1(v1,v2);
else
   proxFunction = @(v1,v2,v3,v4) proxG(v1,v2,v3,v4);
end

if (verbosity > 0)
   fprintf(fid,'%5s  %9s   %9s  %9s  %9s\n','Iter','||r||_2','Gap','Infeas.','Rel. gap');
end


% -------------------------------------------------------------
% Main loop
% -------------------------------------------------------------
while (true)

   % Compute the gradient at f(y)
   if (mod(iter,gradIter) == 0) % Includes first iterations
      r = A*(W.*y)-b;      %      r = A*y-b;
      g = W.*(A'*r);       %      g = A'*r;
      f = (r'*r) / 2;
   else
      r = (Ax + ((tPrev - 1) / t) * (Ax - AxPrev)) - b;
      g = W.*(A'*r);            %      g = A'*r;
      f = (r'*r) / 2;
   end
   
   % Increment iteration count
   iter = iter + 1;

   % Check optimality conditions
   if ((mod(iter,optimIter) == 0))
      % Compute 'dual', check infeasibility and gap
      if (modeLasso)
         infeas = max(norm(g,inf)-lambda,0);

         objPrimal = f + lambda*norm(y,1);
         objDual   = -f - r'*b;
      else
         y2   = [y',0];
         y2   = y2.^2;
         y2   = y2(I2);
         y_I  = sqrt( sum(y2,1) )';         
         y_Is = sort(abs(y_I),'descend');
         g2   = [g',0];
         g2   = g2.^2;
         g2   = g2(I2);
         g_I  = sqrt( sum(g2,1) )';
         g_Is = sort(abs(g_I),'descend');
              
         % Compute primal and dual objective
         objPrimal =  f + lambda'*y_Is;
         objDual   = -f - r'*b;
         infeas    = max(max(cumsum(g_Is-lambda)),0);      
      end
      
      % Format string
      if (verbosity > 0)
         str = sprintf('   %9.2e  %9.2e  %9.2e',objPrimal - objDual, infeas/lambda(1), abs(objPrimal - objDual) / max(1,objPrimal));
      end
      
      % Check primal-dual gap
      if ((abs(objPrimal - objDual)/max(1,objPrimal) < tolRelGap) && ...
          (infeas < tolInfeas * lambda(1)))
         status = STATUS_OPTIMAL;
      end
   else
      str = '';
   end

   if (verbosity > 0)
      if ((verbosity == 2) || ...
         ((verbosity == 1) && (mod(iter,optimIter) == 0)))
      fprintf(fid,'%5d  %9.2e%s\n', iter,f,str);
      end
   end
   
   
   % Stopping criteria
   if (status == 0)
      if (iter >= iterations)
         status = STATUS_ITERATIONS;
      end
   end
   
   if (status ~= 0)
      if (verbosity > 0)
         fprintf(fid,'Exiting with status %d -- %s\n', status, STATUS_MSG{status});
      end
      break;
   end
   
   % Keep copies of previous values
   AxPrev = Ax;
   xPrev  = x;
   fPrev  = f;
   tPrev  = t;
   
   % Lipschitz search
   while (true)
      % Compute prox mapping
      x = proxFunction(y - (1/L)*g, lambda/L, I, I2);
      d = x - y;
      
      Ax = A*(W.*x);          %Ax = A*x;
      r  = Ax-b;
      f  = (r'*r)/2;
      q  = fPrev + d'*g + (L/2)*(d'*d);
      
      Aprods = Aprods + 1;
      
      if (q >= f*(1-1e-12))
         break;
      else
         L = L * eta;
      end
   end
   
   % Update
   t = (1 + sqrt(1 + 4*t^2)) / 2;
   y = x + ((tPrev - 1) / t) * (x - xPrev);
end

% Set solution
x = y;

% Information structure
info = struct();
if (nargout > 1)
   info.runtime   = toc(t0);
   info.Aprods    = Aprods + ceil(iter / gradIter);
   info.ATprods   = ATprods + iter;
   info.objPrimal = objPrimal;
   info.objDual   = objDual;
   info.infeas    = infeas;
   info.status    = status;
   info.L         = L;
end

x  = W.*x;
x2 = [x',0];
x2 = x2.^2;
x2 = x2(I2);
xg = sqrt( sum(x2,1) )';
%xg = sqrt(sum(bsxfun( @times, x2, x2)))';

end % Function MeanGroupAdlas


% ------------------------------------------------------------------------
function opt = getDefaultField(data,field,default)
% ------------------------------------------------------------------------
   if isfield(data,field)
      opt = data.(field);
   else
      opt = default;
   end
end


% ------------------------------------------------------------------------
function x = proxL1(y,lambda)
% ------------------------------------------------------------------------
   % Normalization
   y    = y(:);
   sgn  = sign(y);
   x    = sgn .* max(abs(y) - lambda,0);
end
