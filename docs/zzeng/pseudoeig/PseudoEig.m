
function [lambda,X,C,S,res,lcond] = PseudoEig(A,lambda0,m,k,N)
%
% Compute an m x k pseudo-eigenvalue of a matrix A 
%  Input:       A ---  the matrix
%         lambda0 ---  initial approximation of the eigenvalue
%               m ---  geometric multiplicity
%               k ---  smallest Jordan block size
%               N ---  (optional) matrix whose columns span the same
%                       space spanned by the last m right singular 
%                       vectors of  A-lambda0*I
%
% Output:  lambda ---  the approximate value of the defective eigenvalue
%               X ---  a set of eigenvector and generalized eigenvectors
%             C,S ---  matrix parameters in (A-lambda*I)X = X*S, C'*X = T
%             res ---  residual
%           lcond ---  mxk condition number
%        
   n = size(A,1);   
   if nargin == 4
       C = Srand(n,m);          % set C at random
   else
       [Q,~] = qr(Srand(m,1));  % get a random orthogonal matrix
       C = N*Q;                 % set C at random with the same
                                % column space as N
   end
   
   % calculate initial X and S
   X0 = zeros(n,k);  S = zeros(k,k);
   X0(:,1) = [A-lambda0*eye(n,n); C']\[zeros(n,1); 1; zeros(m-1,1)]; % 1st column
   for j = 1:k-1
       x = [A-lambda0*eye(n,n); C']\[X0(:,j); zeros(m,1)];
       t = 1/norm(x,2);
       X0(:,j+1) = t*x;
       S(j,j+1) = t;
   end
   
   % use the generic Gauss-Newton iteration code
   [Z,res,lcond] = GaussNewton({@peigmap,{1,ones(n,k)},{A,C,S}},@peigjac,{lambda0,X0},0);
   lambda = Z{1};  X=Z{2};
end
