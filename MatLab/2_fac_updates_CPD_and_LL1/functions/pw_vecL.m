function out = pw_vecL(U,R,L,varargin)

% Fonction pour réaliser l'opération que je note \odot_{\text{vec}}
% entre les facteurs A et B avec des L_r différents
% pw_vec veut dire "partition-wise vectorization"

out = [];

for r=1:R
    % if r==1
        %out = X(:,1:L(1))*Y(:,1:L(1))'; out = out(:);
    %else
        %mat = X(:,sum(L(:,1:r-1))+1:sum(L(:,1:r)))*Y(:,sum(L(:,1:r-1))+1:sum(L(:,1:r)))';
        %out = [out mat(:)];
        %out = [out kr(X(:,(r-1)*L(r)+1:r*L(r)),Y(:,(r-1)*L(r)+1:r*L(r)))*ones(L(r),1)];
    % end
    % mat = U{r}{1}*U{r}{2}';
    % mat = mat';
    % out = [out mat(:)];
    out = [out kr(U{r}{2},U{r}{1})*ones(L(r),1)];
    
end

end