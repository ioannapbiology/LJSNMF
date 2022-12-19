function [normalized_cut_matrix] = weight_to_norm_cut(weight_matrix)
%weight_to_norm_cut
%   from a weight/similarity matrix
%   we compute A = (D^-1/2)*W*(D^-1/2), the normalized cut matrix
%   (from LJ-SNMF)

D = diag(sum(weight_matrix,2));
s = size(D);
index = 1:s(1)+1:s(1)*s(2);  % Indices of the main diagonal
D_m_sqr = D; % initialize
D_m_sqr(index) = D(index) .^ (-1/2); % raise only diagonal elements to the power of (-1/2)
normalized_cut_matrix = D_m_sqr*weight_matrix*D_m_sqr; 


end

