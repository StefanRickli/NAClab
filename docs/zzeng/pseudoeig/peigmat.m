
function M = peigmat(lambda,X,A,C,S)
%
% Compute an m x k pseudo-eigenvalue of a matrix A 

   % get the matrix for the linear transformation
   M = LinearTransformMatrix(@pemap,{ones(size(X)),1},{lambda,X,A,C,S}); 
end

function T = pemap(Y,lambda,lambda0,X0,A,C,S)
%
%  The Jacobian of the pseudo-eig. mapping
%    with components rearranged 
%
   X = Y(:,end:-1:1); Y0 = X0(:,end:-1:1);        % reverse the columns of X, X0
   
   U = -lambda*Y0 + A*X - lambda0*X - X*S;   % the Jacobian linear transformation
   V = C'*X;  
   
   T = [V(:,end:-1:1);U(:,end:-1:1)];  % reverse the columns of the image
end
