function [x, L] = katz_edge(A, x, alpha, L_max, tol, u, v)
    % KATZ_EDGE Gives an approximate update of the Katz centrality
    %     score after removing an edge e = {u, v}.
    %
    % INPUTS:
    %     A      : adjacency matrix, in R^{n x n};
    %     x      : original Katz vector;
    %     alpha  : Katz parameter;
    %     L_max  : maximum number of iterations;
    %     u      : source of the edge e;
    %     v      : target of the edge e.
    %
    % OUTPUTS:
    %     x      : updated Katz vector;
    %     L      : final number of iterations.

    x_u = x(u);
    x_v = x(v);
    x(u) = x_u - alpha * x_v;
    x(v) = x_v - alpha * x_u;

    L = 1;

    s = alpha^2 * A(:, u);
    s(v) = s(v) - alpha^2;

    t = alpha^2 * A(:, v);
    t(u) = t(u) - alpha^2;

    x = x - x_v * s - x_u * t;

    while norm(x_v * s + x_u * t) / norm(x) > tol && L < L_max
        s_u = s(u);
        s_v = s(v);
        s = alpha * A * s;
        s(v) = s(v) - alpha * s_u;
        s(u) = s(u) - alpha * s_v;

        t_u = t(u);
        t_v = t(v);
        t = alpha * A * t;
        t(v) = t(v) - alpha * t_u;
        t(u) = t(u) - alpha * t_v;

        x = x - x_v * s - x_u * t;
        L = L + 1;
    end
end