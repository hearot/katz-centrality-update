function [x, L] = katz_node(A, x, alpha, L_max, tol, w)
    % KATZ_NODE Gives an approximate update of the Katz centrality
    %     score after removing a node w.
    %
    % INPUTS:
    %     A      : adjacency matrix, in R^{n x n};
    %     x      : original Katz vector;
    %     alpha  : Katz parameter;
    %     L_max  : maximum number of iterations;
    %     w      : node to be removed.
    %
    % OUTPUTS:
    %     x      : updated Katz vector;
    %     L      : final number of iterations.

    L = 1;
    q = alpha * A(:, w);
    x_w = x(w);
    x = x - x_w * q;

    while x_w * norm(q) / norm(x) > tol && L < L_max
        q = alpha * A * q;
        q(w) = 0;
        x = x - x_w * q;
        L = L + 1;
    end

    x(w) = 1;
end