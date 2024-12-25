function out = pw_vecL(U,R,L,varargin)

% Fonction for computing [(A_1 \odot B_1)1_{L_1},...,(A_R \odot B_R)1_{L_R}]
% with A_r = U{r}{1}, B_r = U{r}{2} for r=1,...,R

out = [];

for r=1:R
    out = [out kr(U{r}{2},U{r}{1})*ones(L(r),1)];
    
end

end