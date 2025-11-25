function [status, optval, state_history, control_history] = powered_descent(r0, rdot0, rf, t_f, flag)
    g = [-3.7114 0 0]'; % [m/s^]
    
    m_wet = 1905; % [kg]
    m_dry = 1505;
    
    burn_rate = 4.53*10^-4; % [s/m]

    rho1 = 4972; % [N]
    rho2 = 13260; % [N]

    y0 = [r0; rdot0; log(m_wet)];
    
    glideslope_angle_bound = deg2rad(86);

    N = 100; % [steps]
    delta_t = t_f/N; % [s]

    A_c = [zeros(3,3) eye(3) zeros(3,1); zeros(4,7)];

    B_c = [zeros(3,4); eye(3), zeros(3,1); zeros(1,3) -burn_rate];

    [A_d, B_d] = c2d(A_c, B_c, delta_t);

    [xi_batched, Psi_batched] = build_batch_matrices(A_d, B_d, N, y0, g);
    
    S = [zeros(2, 1), eye(2)];
    c = [-tan(glideslope_angle_bound); zeros(2, 1)];

    omega = zeros(1, 4*N);
    for M = 4:4:4*N
        omega(1, M) = delta_t;
    end
    
    cvx_precision best
    cvx_begin quiet
        variable eta(4*(N+1), 1)

        k = (0:N)';

        z0 = log(m_wet - rho2 * burn_rate * delta_t * k);
        mu_1 = rho1 * exp(-z0);
        mu_2 = rho2 * exp(-z0);

        u_matrix = reshape(eta, 4, N+1);
        u_k = u_matrix(1:3, :);
        sigma_k = u_matrix(4, :)';

        y_k = [y0 reshape(xi_batched(:) + Psi_batched * eta, 7, N)];

        y_N = y_k(:, end);

        z_k = y_k(7, :)';

        if strcmp(flag, 'Fuel')
            minimize omega*eta(1:4*N, 1)
            subject to
                y_N(1:6) == [rf; zeros(3, 1)];
            
            pos_error = y_k(1:3, :) - repmat(rf, 1, N+1);

        elseif strcmp(flag, 'Landing Error')
            minimize norm(y_N(1:3))
            subject to
                y_N(1) == 0;
                y_N(4:6) == 0;

            pos_error = y_k(1:3, :) - repmat(y_N(1:3), 1, N+1);
        end
            
        subject to
            y_N(7) >= log(m_dry);

            norms(u_k, 2, 1)' <= sigma_k;

            mu_1 .* (1 - (z_k - z0) + 0.5 * (z_k - z0).^2) <= sigma_k <= mu_2 .* (1 - (z_k - z0));

            log(m_wet - rho2 * burn_rate * delta_t * k) <= z_k <= log(m_wet - rho1 * burn_rate * delta_t * k);
            
            norms(S * pos_error, 2, 1) + c' * pos_error <= 0;
    cvx_end
    
    status = cvx_status;
    optval = cvx_optval;
    control_history = reshape(eta, 4, N+1);
    state_history = [y0 reshape(xi_batched(:) + Psi_batched * eta, 7, N)];
end
%[text] 
function [xi_batched, Psi_batched] = build_batch_matrices(A, B, N, y0, g)
    xi_batched = zeros(7, N);
    Psi_batched = zeros(7*N, 4*(N+1));
    
    current_Lambda = B;
    
    for k = 1:N
        xi_batched(:, k) = A^k * y0 + current_Lambda * [g; 0];
    
        if k < N
            current_Lambda = A * current_Lambda + B;
        end
    
        if k == 1
            Psi_batched(1:7, 1:4) = B;
        else
            Psi_batched(7*k-6:7*k, :) = A * Psi_batched(7*(k-1)-6:7*(k-1), :);
            Psi_batched(7*k-6:7*k, 4*k-3:4*k) = B;
        end
    end
end

%[appendix]{"version":"1.0"}
%---
