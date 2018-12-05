function [Q,s,t] = fit_isometry(X,Y)

% Fits a mapping with rotation (Q), stretch (s), and translation (t)
% between the n corresponding points defined by the matrices X and Y.
% Mapping structure is Y = s*Q*X + t.
% 
% Input:
% X: m x n matrix, where each column defines the location of a point
% Y: m x n matrix, where column i defines the location of the point
% corresponding to column i of X
% 
% Output:
% Q: m x m rotation matrix
% s: positive isotropic stretch factor
% t: m x 1 translation

if any(size(X) ~= size(Y))
    error('X and Y must have the same size.');
end

n = size(X,2);
Xhat = X - repmat(mean(X,2),1,n);
Yhat = Y - repmat(mean(Y,2),1,n);
[U,~,V] = svd(Xhat*Yhat');
Q = V*U'; %Rotation
if det(Q) < 0
    Q = Q - 2*U(:,end)*V(:,end)'; %If reflection, invert last dyad
end
s = trace(Xhat*Yhat'*Q)/norm(Xhat,'fro')^2; %Stretch
t = mean(Y,2) - s*Q*mean(X,2); %Translation

end