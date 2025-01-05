
function Z = peigmap(lambda,X,A,C,S)
%  Function for the mapping 
%       (lambda,X) --> {(A-lambda*I)*X - X*S, C'*X-T}
%  where entries of T are all zeros except T(1,1) = 1
%     A, C, S are constant matrices
% 
    U = A*X-lambda*X-X*S;          % the first component
    V = C'*X; V(1,1) = V(1,1)-1;   % the second component
    Z = {U,V};                     % two components in a cell
end
