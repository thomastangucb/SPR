%% nonlinear embedding preserving multiple local-linearities

%Input: the data set X, the neighborhood size k or the selected neighborhood NI, the reduced dimension d
%the neighborhood NI{i} include the indices of x_i and its neighbors
%Jing Wang and Zhenyue Zhang, Nonlinear Embedding Preserving Multiple Local-linearities, Patten Recognition, Vol.43, pp.1257¡ª1268, 2010

function [Phi, Y] = MLLE(X, k, d)

[~, N]=size(X);    % X:m*N, m dimension, N points
tol = 1.0e-3;     % the regularlization parameter

% Neighborhood selection
X2 = sum(X.^2,1);
D = repmat(X2,N,1)+repmat(X2',1,N)-2*X'*X;
[~, J] = sort(D,2);


% Local Information
for i=1:N
    xi = X(:,i);
    Ii = J(i,1:k+1); Ji = Ii(2:end);
    Gi = X(:,Ji)-repmat(xi,[1,k]);        % or [Ui,Si,Vi] = svd(Gi);
    [V,S] = schur(full(Gi'*Gi));
    [ei,JIi] = sort(diag(S),'descend');
    ratio(i) = sum(ei(d+1:k))/sum(ei(1:d));
    Theta{i} = V(:,JIi);                  % the local coordinates system
    Ev{i} = ei;
    C = Gi'*Gi;
    C = C + eye(k,k)*tol*trace(C);        % regularlization
    Cw = C\ones(k,1);                     % solve C*Cw=1
    Cw = Cw/sum(Cw);                      % enforce sum(Cw)=1
    Ow{i} = Cw;                           % considered as the optimal weight vector
end
temp = sort(ratio); eta = temp(ceil(N/2));

% Determine the number of weights
s = zeros(1,N);
for i=1:N
    ell = k-d;
    Lambda = Ev{i};
    while sum(Lambda(k-ell+1:end))/sum(Lambda(1:k-ell))>eta && ell>1
        ell = ell-1;
    end
    s(i) = ell;
end

Phi = sparse([],[],[],N,N,0);
for i = 1:N
    Ii = J(i,1:k+1);
    Vi = Theta{i};
    Ve = Vi(:,k-s(i)+1:end); ve = sum(Ve,1)';
    alpha = norm(ve)/sqrt(s(i));
    u = ve-alpha; normu = norm(u);
    if normu > 1.0e-5
        u = u/normu;
        Wi = (1-alpha)^2*Ow{i}*ones(1,s(i))+(2-alpha)*(Ve-(Ve*(2*u))*u');    % the multiple local weights
    else
        Wi = (1-alpha)^2*Ow{i}*ones(1,s(i))+(2-alpha)*Ve;                    % the multiple local weights
    end
    Phi(Ii,Ii) = Phi(Ii,Ii)+[-ones(1,s(i));Wi]*[-ones(1,s(i));Wi]';          % Construct the alignment matrix
end

%
if nargout > 1
    options.disp = 0; options.isreal = 1; options.issym = 1;
    [Q, Lamb] = eigs(sparse(Phi),d+1,0,options);  lambda = diag(Lamb);
    [lambda_s, II] = sort(abs(lambda));
    Y = Q(:,II(2:d+1))';
end
