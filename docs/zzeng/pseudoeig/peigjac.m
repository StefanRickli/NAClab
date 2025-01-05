
function Z = peigjac(lambda,X,lambda0,X0,A,C,S)
%
% The Jacobian of the mapping 
%           (lambda,X) --> ((A-lambda*I)*X-X*S, C'*X-T)
%    at (lambda0,X0), where A, C, S, T are constant matrices.
% This Jacobian is a linear transformation that maps (lambda,X) to the cell
% 
    Z = {-lambda*X0+A*X-lambda0*X-X*S;  C'*X};
end
