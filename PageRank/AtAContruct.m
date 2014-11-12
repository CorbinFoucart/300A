% Function to construct (A^TA)^-1
% Corbin Foucart

% this function works recursively

% reads in a matrix with only the first diagonal entry complete
function rX = AtAContruct(R,X,col_num)

n = length(X(:,1));
R
X
col_num

% base case, we have completed the process, return X
if (col_num == 1)
    rX = X;    
else 
    % not the base case, we must build the last row and column
    % of sub-matrix C, which is a square matrix with the corner
    % filled in. we must find X_ij for i = c_cols - 1 : 1
    
    % build columns and rows
    j = col_num;
    for i = (j - 1):-1:1
        % compute x_ij
        total = 0;
        for alpha = (i+1):n
            total = total - R(i,alpha)*X(alpha,j);
        end
        total = total./(X(j,j));
        X(i,j) = total
        X(j,i) = total
    end

    % compute next diagonal entry c_ j-1 j-1
    k = col_num-1;
    total = 0;
    for alpha = (k + 1):n
        total = total - R(k, alpha)*X(alpha, k);
    end
    total = total./(R(k,k)).^2;
    X(k,k) = total;
    
    X = AtAContruct(R, X, (col_num-1));
    rX = X;
end

end