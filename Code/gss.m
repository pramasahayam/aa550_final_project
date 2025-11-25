function [opt_t_f, optval, state_history, thrust_history] = gss(lb, ub, tol, r0, rdot0, rf, flag)
gr = (sqrt(5) - 1)/2;
d = gr * (ub - lb);

x1 = lb + d;
x2 = ub - d;

while d > tol
    [x1_status, x1_optval, ~] = powered_descent(r0, rdot0, rf, x1, flag);
    [x2_status, x2_optval, ~] = powered_descent(r0, rdot0, rf, x2, flag);

    if ~strcmp(x1_status, 'Solved')
        if strcmp(x1_status, 'Failed') || strcmp(x1_status, 'Infeasible')
            x1_optval = 10^6;
        end
    end

    if ~strcmp(x2_status, 'Solved')
        if strcmp(x2_status, 'Failed') || strcmp(x2_status, 'Infeasible')
            x2_optval = 10^6;
        end
    end

    formatSpec = 'lb = %.2f, x1 = %.2f (%.4f), x2 = %.2f (%.4f), ub = %.2f\n';
    fprintf(formatSpec, lb, x1, x1_optval, x2, x2_optval, ub);

    if x1_optval < x2_optval
        lb = x2;
        d = gr * (ub - lb);
        x2 = x1;
        x1 = lb + d;
    elseif x1_optval > x2_optval
        ub = x1;
        d = gr * (ub - lb);
        x1 = x2;
        x2 = ub - d;
    end

    opt_t_f = (x1 + x2)/2;
end

fprintf('----------------------------------------------------------------------\nFinal Values\n');
formatSpec = 'lb = %.2f, x1 = %.2f, x2 = %.2f, ub = %.2f\n';
fprintf(formatSpec, lb, x1, x2, ub);

[~, optval, state_history, thrust_history] = powered_descent(r0, rdot0, rf, opt_t_f, flag);
fprintf('Optimal t_f: %.2f [s], Optimal Value: %.4f', opt_t_f, optval);
end

%[appendix]{"version":"1.0"}
%---
