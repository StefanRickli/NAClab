
function [lambda,X,C,S,res,lcond] = PEigRefine(A,lambda0,X0,C0,S0)
%
% Refine the pseudo-eigenvalue computed by PseudoEig.m
%  Input:       A ---  the matrix
%         lambda0 ---  the output lambda of PseudoEig.m
%              X0 ---  the output X of PseudoEig.m
%           C0,S0 ---  the output C and S of PseudoEig.m
%
% Output:  lambda ---  the approximate value of the defective eigenvalue
%               X ---  a set of eigenvector and generalized eigenvectors
%             C,S ---  matrix parameters in (A-lambda*I)X = X*S, C'*X = T
%             res ---  residual
%           lcond ---  mxk condition number
%             
   n = size(C0,1);   % get matrix size and geo. multiplicity
   k = size(S0,1);     % get the smallest Jordan block size
   
   C = C0; C(:,1) = X0(:,1)/norm(X0(:,1)); % unit eigenvector as C(:,1)
   
   Y = X0; Y(:,1) = C(:,1);  t = [0,Y(:,1)'*Y(:,2:k)]; % reset X0 in Y so that
   Y(:,2:k) = Y(:,2:k)-Y(:,1)*t(2:k);                  % C(:,1)'*Y(:,2:end) = 0
                                                       % so that C'*Y = T
   S = S0; S(1,2) = S(1,2)*norm(X0(:,1));  % update S so that
   for j = 3:k                             % (A-lambda0*I)*Y = Y*S
       S(1,j) = S(j-1,j)*t(j-1);
   end;
   
   % Orthonormalize Y as X and update S to S0
   [X,R] = qr(Y,0);
   if R(1,1) < 0, R(1,:) = -R(1,:); X(:,1) = -X(:,1); end; % so X(:,1) ~= -Y(:,1)
   S = R*S*inv(R);

   % use the generic Gauss-Newton iteration code
   [Z,res,lcond] = GaussNewton({@peigmap,{1,ones(n,k)},{A,C,S}},@peigjac,{lambda0,X},1);
   lambda = Z{1};  X=Z{2};
end
