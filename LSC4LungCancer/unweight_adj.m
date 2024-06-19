function adj_matrix = unweight_adj(A)
% Extract row and column indices from the weighted adjacency matrix
[row, col] = find(A>0.3);

% Create an unweighted adjacency matrix using sparse
n = max(max(row), max(col));  % Number of nodes
adj_matrix = sparse(row, col, ones(size(row)), n, n);
end

