function [status, optval, state_history, thrust_history] = solve_fuel_optimal(r0, rdot0, rf, t_f)
%[text] $A=\\left\\lbrack \\begin{array}{ccc}\n0 & I\_{3\\times 3}  & 0\\\\\n0 & 0 & 0\\\\\n0 & 0 & 0\n\\end{array}\\right\\rbrack \\in {\\mathbb{R}}^{7\\times 7}${"editStyle":"visual"}
%[text] $B=\\left\\lbrack \\begin{array}{cc}\n0 & 0\\\\\nI\_{3\\times 3}  & 0\\\\\n0 & -\\alpha \n\\end{array}\\right\\rbrack \\in {\\mathbb{R}}^{7\\times 4}${"editStyle":"visual"}
    g = [-3.7114 0 0]'; % [m/s^]
    
    m_wet = 1905; % [kg]
    
    burn_rate = 4.53*10^-4; % [s/m]
    n_thrusters = 6;
    T1 = 0.3*3100; % [kN]
    T2 = 0.8*3100; % [kN]
    thrust_angle = deg2rad(4);
    rho1 = n_thrusters * T1 * cos(thrust_angle);
    rho2 = n_thrusters * T2 * cos(thrust_angle);

    y0 = [r0; rdot0; log(m_wet)];

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

    e_sigma = [0 0 0 1]';
    E = [eye(6), zeros(6,1)];
    F = [zeros(1,6), 1];
    E_u = [eye(3), zeros(3, 1)];

    omega = zeros(1, 4*N);
    for M = 4:4:4*N
        omega(1, M) = delta_t;
    end

    cvx_solver sdpt3
    cvx_begin quiet
        variable eta(4*(N+1), 1)
        minimize omega*eta(1:4*N, 1)
        subject to
        for k = 0:N
            if k == 0
                Upsilon_k = Upsilon_batched(4*k+1:4*k+4, :);
                xi_Psi_eta = xi_batched(:, N) + Psi_batched(7*N-6:7*N, :) * eta;
    
                mu_1(1) * (1 - (F*y0 - z0(1)) + 0.5*(F*y0 - z0(1))^2) ...
                    <= e_sigma' * Upsilon_k * eta <= mu_2(1) * (1 - (F*y0 - z0(1)))
    
                norm(E_u * Upsilon_k * eta, 2) <= e_sigma' * Upsilon_k * eta
    
                [eye(6) zeros(6, 1); zeros(1, 7)] * xi_Psi_eta == [rf ;zeros(4, 1)]
            else
                Upsilon_k = Upsilon_batched(4*k+1:4*k+4, :);
                Fxi_Psi_eta = F * (xi_batched(:, k) + Psi_batched(7*k-6:7*k, :) * eta);
    
                norm(E_u * Upsilon_k * eta, 2) <= e_sigma' * Upsilon_k * eta
    
                mu_1(k+1) * (1 - (Fxi_Psi_eta - z0(k+1)) + 0.5*(Fxi_Psi_eta - z0(k+1))^2) ...
                    <= e_sigma' * Upsilon_k * eta <= mu_2(k+1) * (1 - (Fxi_Psi_eta - z0(k+1)))
    
                log(m_wet - rho2 * burn_rate * delta_t * k) ...
                    <= Fxi_Psi_eta <= log(m_wet - rho1 * burn_rate * delta_t * k)
            end
        end
    cvx_end

    optval = cvx_optval;
    status = cvx_status;

    state_history = zeros(7, N+1);
    thrust_history = zeros(3, N+1);
    % thrust_mag = zeros(1, N+1);

    for k = 0:N
        if k == 0
            state_history(:, 1) = y0;
            thrust_history(:, 1) = zeros(3, 1);
        else
            xi_Psi_eta = xi_batched(:, k) + Psi_batched(7*k-6:7*k, :) * eta;
            Upsilon_k = Upsilon_batched(4*k+1:4*k+4, :);

            state_history(:, k+1) = xi_Psi_eta;
            thrust_history(:, k+1) = [eye(3) zeros(3,1)] * Upsilon_k * eta;
            % thrust_mag(:, k+1) = norm(thrust_history(k+1), 2);
        end
    end
end
%[text] $\\Phi\_k =A^k${"editStyle":"visual"}
%[text] $\\Lambda\_k =B+\\textrm{AB}+\\ldotp \\ldotp \\ldotp +A^{k-1} B=A\\left(\\Lambda\_{k-1} \\right)+B${"editStyle":"visual"}
%[text] $\\xi\_k =\\Phi\_k y\_0 +\\Lambda\_k \\left\\lbrack \\begin{array}{c}\ng\\\\\n0\n\\end{array}\\right\\rbrack${"editStyle":"visual"}
%[text] $\\Psi\_k =\\textrm{this}\\;\\textrm{shit}\\;\\textrm{tuff}\\;\\textrm{bruh}\\;\\textrm{idk}${"editStyle":"visual"}
%[text] $\\Upsilon\_k =\\left\\lbrack \\begin{array}{ccccccc}\n0\_{4\\textrm{x4}}  & \\ldotp \\ldotp \\ldotp  & 0\_{4\\textrm{x4}}  & I\_{4\\textrm{x4}}  & 0\_{4\\textrm{x4}}  & \\ldotp \\ldotp \\ldotp  & 0\_{4\\textrm{x4}} \n\\end{array}\\right\\rbrack${"editStyle":"visual"}
%[text] The batched $\\xi${"editStyle":"visual"}, $\\Psi${"editStyle":"visual"}, and $\\Upsilon${"editStyle":"visual"} matrices over all temporal nodes are found ahead of time, as required by CVX.
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
