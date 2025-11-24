function [status, optval, state_history, thrust_history] = powered_descent(r0, rdot0, rf, t_f, flag)
    g = [-3.7114 0 0]'; % [m/s^]
    
    m_wet = 1905; % [kg]
    m_dry = 1505;
    
    burn_rate = 4.53*10^-4; % [s/m]
    n_thrusters = 6;
    T1 = 0.3*3100; % [kN]
    T2 = 0.8*3100; % [kN]
    cant_angle = deg2rad(27); % [rad]
    
    % rho1 = n_thrusters * T1 * cos(cant_angle);
    rho1 = 4972; % [N]
    % rho2 = n_thrusters * T2 * cos(cant_angle);
    rho2 = 13260; % [N]

    y0 = [r0; rdot0; log(m_wet)];
    
    glideslope_angle_bound = deg2rad(86);

    N = 100; % [steps]
    delta_t = t_f/N; % [s]

    A_c = [zeros(3,3) eye(3) zeros(3,1);
           zeros(4,7)];

    B_c = [zeros(3,4);
           eye(3), zeros(3,1);
           zeros(1,3) -burn_rate];

    [A_d, B_d] = c2d(A_c, B_c, delta_t);

    [xi_batched, Psi_batched, Upsilon_batched] = build_batch_matrices(A_d, B_d, N, y0, g);

    [z0, mu_1, mu_2] = precompute_vars(N, delta_t, m_wet, burn_rate, rho1, rho2);
    
    S_gs = [zeros(2, 1), eye(2), zeros(2, 3)];
    c_gs = [-tan(glideslope_angle_bound); zeros(5, 1)];

    omega = zeros(1, 4*N);
    for M = 4:4:4*N
        omega(1, M) = delta_t;
    end

    if strcmp(flag, 'Fuel')
        cvx_precision best
        cvx_begin
            variable eta(4*(N+1), 1)
            minimize omega*eta(1:4*N, 1)
            subject to
                
            y_N = xi_batched(:, N) + Psi_batched(7*N-6:7*N, :) * eta;
            y_N(1:6) == [rf; zeros(3, 1)];
            y_N(7) >= log(m_dry);

            for k = 0:N
                u_k = Upsilon_batched(4*k+1:4*k+4, :) * eta;
                norm(u_k(1:3)) <= u_k(4);

                if k == 0
                    mu_1(1) * (1 - (y0(7) - z0(1)) + 0.5*(y0(7) - z0(1))^2) ...
                        <= u_k(4) <= mu_2(1) * (1 - (y0(7) - z0(1)));
                else
                    y_k = xi_batched(:, k) + Psi_batched(7*k-6:7*k, :) * eta;
    
                    mu_1(k+1) * (1 - (y_k(7) - z0(k+1)) + 0.5*(y_k(7) - z0(k+1))^2) ...
                        <= u_k(4) <= mu_2(k+1) * (1 - (y_k(7) - z0(k+1)));
    
                    log(m_wet - rho2 * burn_rate * delta_t * k) ...
                        <= y_k(7) <= log(m_wet - rho1 * burn_rate * delta_t * k);
    
                    norm(S_gs * y_k(1:6)) + c_gs' * y_k(1:6) <= 0;
                end
            end
        cvx_end

    elseif strcmp(flag, 'Landing Error')
        cvx_precision best
        cvx_begin
            variable eta(4*(N+1), 1)

            y_N = xi_batched(:, N) + Psi_batched(7*N-6:7*N, :) * eta;

            minimize norm(y_N(1:3))
            subject to

            y_N(1) == 0;
            y_N(4:6) == zeros(3, 1);
            y_N(7) >= log(m_dry);

            for k = 0:N
                u_k = Upsilon_batched(4*k+1:4*k+4, :) * eta;
                norm(u_k(1:3)) <= u_k(4);

                if k == 0
                    mu_1(1) * (1 - (y0(7) - z0(1)) + 0.5*(y0(7) - z0(1))^2) ...
                        <= u_k(4) <= mu_2(1) * (1 - (y0(7) - z0(1)));
                else
                    y_k = xi_batched(:, k) + Psi_batched(7*k-6:7*k, :) * eta;

                    mu_1(k+1) * (1 - (y_k(7) - z0(k+1)) + 0.5*(y_k(7) - z0(k+1))^2) ...
                        <= u_k(4) <= mu_2(k+1) * (1 - (y_k(7) - z0(k+1)));

                    log(m_wet - rho2 * burn_rate * delta_t * k) ...
                        <= y_k(7) <= log(m_wet - rho1 * burn_rate * delta_t * k);
                    
                    pos_error = [y_k(1:3) - y_N(1:3); zeros(3, 1)];

                    norm(S_gs * pos_error) + c_gs' * (pos_error) <= 0;
                end
            end
        cvx_end
    end
    
    status = cvx_status;
    optval = cvx_optval;

    state_history = zeros(7, N+1);
    thrust_history = zeros(3, N+1);
    % thrust_mag = zeros(1, N+1);

    for k = 0:N
        if k == 0
            state_history(:, 1) = y0;
            thrust_history(:, 1) = zeros(3, 1);
        else
            y_k = xi_batched(:, k) + Psi_batched(7*k-6:7*k, :) * eta;
            u_k = Upsilon_batched(4*k+1:4*k+4, :);

            state_history(:, k+1) = y_k;
            thrust_history(:, k+1) = [eye(3) zeros(3,1)] * u_k * eta;
            % thrust_mag(:, k+1) = norm(thrust_history(k+1), 2);
        end
    end
end
%[text] 
function [xi_batched, Psi_batched, Upsilon_batched] = build_batch_matrices(A, B, N, y0, g)
    xi_batched = zeros(7, N);
    Psi_batched = zeros(7*N, 4*(N+1));
    Upsilon_batched = zeros(4*N, 4*(N+1));
    
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
    
    for k = 1:N+1
        Upsilon_batched(4*k-3:4*k, 4*k-3:4*k) = eye(4);
    end
end
%[text] 
function [z0, mu_1, mu_2] = precompute_vars(N, delta_t, m_wet, burn_rate, rho1, rho2)
    z0 = zeros(N+1, 1);
    mu_1 = zeros(N+1, 1);
    mu_2 = zeros(N+1, 1);
    for k = 0:N
        z0(k+1) = log(m_wet - rho2 * burn_rate * delta_t * k);
        mu_1(k+1) = rho1 * exp(-z0(k+1));
        mu_2(k+1) = rho2 * exp(-z0(k+1));
    end
end

%[appendix]{"version":"1.0"}
%---
