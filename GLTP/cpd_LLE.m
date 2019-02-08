%CPD_LLE The non-rigid CPD point-set registration. It is recommended to use
%   and umbrella function rcpd_register with an option opt.method='nonrigid'
%   instead of direct use of the current funciton.
%
%   [C, W, sigma2, iter, T] =cpd_MLLE(X, Y, beta, lambda, max_it, tol, viz, outliers, fgt, corresp, sigma2, knn, tau);
%
%   Input
%   ------------------
%   X, Y       real, double, full 2-D matrices of point-set locations. We want to
%              align Y onto X. [N,D]=size(X), where N number of points in X,
%              and D is the dimension of point-sets. Similarly [M,D]=size(Y).
%   beta       (try 1..5 ) std of Gaussian filter (Green's funciton) 
%   lambda     (try 1..5) regularization weight
%   max_it     (try 150) maximum number of iterations allowed
%   tol        (try 1e-5) tolerance criterium
%   viz=[0 or 1]    Visualize every iteration         
%   outliers=[0..1] The weight of noise and outliers, try 0.1
%   fgt=[0 or 1]    (default 0) Use a Fast Gauss transform (FGT). (use only for the large data problems)
%   corresp=[0 or 1](default 0) estimate the correspondence vector.
%
%
%   Output
%   ------------------
%   C      Correspondance vector, such that Y corresponds to X(C,:).
%   W      Non-rigid transformation cooeficients.
%   sigma2 Final sigma^2
%   iter   Final number or iterations
%   T      Registered Y point set
%


function  [C, W, sigma2, iter, T] =cpd_LLE(X, Y, beta, lambda0, max_it, tol, viz, outliers, fgt, corresp, sigma2, knn, tau0 , tau_annealing_handle, lambda_annealing_handle)

[N, D]=size(X); [M, D]=size(Y);

% Initialization
iter = 0;  ntol = tol + 10; W = zeros(M,D);
if ~exist('sigma2','var') || isempty(sigma2) || (sigma2==0)
    sigma2 = (M*trace(X'*X)+N*trace(Y'*Y)-2*sum(X)*sum(Y)')/(M*N*D);
end
sigma2_init = sigma2;
T = Y; 

% Construct affinity matrix G
G = cpd_G(Y, Y, beta);

% Local Linear Embedding
Phi = LLE(Y', knn, D);

iter=0; ntol=tol+10; L=1;
while (iter<max_it) && (ntol > tol) && (sigma2 > 1e-8) 
    
    % Plot the result on current iteration
    if viz, cpd_plot_iter(X, T); end;
    
    % Simulated Annealing
    %tau = tau0*(1-iter/max_it);
    tau_annealing_ratio = tau_annealing_handle(iter, max_it) ; %1-iter/max_it; %1- 1/(1 + exp(-20*(iter/max_it-0.9)));
    tau = tau0*tau_annealing_ratio;
    
    lambda_annealing_ratio = lambda_annealing_handle(iter, max_it) ; %1- 1/(1 + exp(-20*(iter/max_it-0.9)));
    lambda = lambda0*lambda_annealing_ratio;
    
    L_old = L;
    % Check wheather we want to use the Fast Gauss Transform
    if (fgt==0)  % no FGT
        [P1,Pt1, PX, L]=cpd_P(X,T, sigma2 ,outliers); st='';
    else         % FGT
        [P1, Pt1, PX, L, sigma2, st]=cpd_Pfast(X, T, sigma2, outliers, sigma2_init, fgt);
    end
    
    %L = L + lambda/2*trace(W'*G*W);
    L = L + lambda/2*trace(W'*G*W) + tau/2*trace(T'*Phi*T);
    ntol=abs((L-L_old)/L);
    
    if viz
        disp([' GLTP nonrigid ' st ' : dL= ' num2str(ntol) ', iter= ' num2str(iter) ' sigma2= ' num2str(sigma2)]);
    end
    
    % M-step. Solve linear system for W.
    dP1 = spdiags(P1,0,M,M); % precompute diag(P)
    %W=(dP*G+lambda*sigma2*eye(M))\(PX-dP*Y);
    W = (dP1*G + lambda*sigma2*eye(M) + tau*sigma2*Phi*G) \ (PX-dP1*Y - tau*sigma2*Phi*Y);
    
    % % same, but solve symmetric system, this can be a bit faster
    % % but can have roundoff errors on idP step. If you want to speed up
    % % use rather a lowrank version: opt.method='nonrigid_lowrank'.
    %
    % idP=spdiags(1./P1,0,M,M); 
    % W=(G+lambda*sigma2*idP)\(idP*PX-Y)

    % update Y postions
    T = Y + G*W;

    Np=sum(P1);sigma2save=sigma2;
    sigma2=abs((sum(sum(X.^2.*repmat(Pt1,1,D)))+sum(sum(T.^2.*repmat(P1,1,D))) -2*trace(PX'*T)) /(Np*D));
    
    iter=iter+1;
end

if viz, disp(['GLTP registration succesfully completed: iter ', num2str(iter),', ntol: ', num2str(ntol), ', sigma2: ',num2str(sigma2)] ); end;

%Find the correspondence, such that Y corresponds to X(C,:)
if corresp, C = cpd_Pcorrespondence(X,T,sigma2save,outliers); else C=0; end;
