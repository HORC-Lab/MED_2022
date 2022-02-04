function [Jx,Ju] = jac_comp(x_0_star,u_0_star)
% function [Jx,Ju] = jac_comp(x_0_star,u_0_star)

% Description:
% This function approximates the Jacobians Jx and Ju using finite
% difference, evaluated at the optimal pair (x_0_star,u_0_star). Initially,
% the system is simulated for one step, using the optimal MS state and control 
% input vectors, and the subsequent MS state is derived. Then, each one of
% the optimal MS states is separately perturbed by a really small value,
% and the system is again simulated for one step to derive the subsequent
% MS state. After that, each one of the optimal control inputs is separately 
% perturbed by a really small value, and the system is again simulated for 
% one step to derive the subsequent MS state. Finally, the Jacobians are
% calculated using finite difference, and after being multiplied by the
% matrix A, they are returned by the function back to main program. For more
% information, the reader is referred to the work by Liu et al., described
% in the publication "Dynamic walking...", DOI: 10.1109/ICRA.2015.7139999,
% and the corresponding dissertation by Yiping Liu.

% INPUTS:
% x_0_star: optimal initial MS state for the specific forward speed as derived from the nonlinear optimizer
% u_0_star: optimal control input vector for the specific forward speed as derived from the nonlinear optimizer

% OUTPUTS:
% Jx: scaled approximated Jacobian of return map with respect to the state vector
% Ju: scaled approximated Jacobian of return map with respect to the control input vector

% NOTE:
% No need to declare global variables here. The internal functions
% "sim_ms_to_ms_compliant" and "sim_td_ms_compliant" can access the global 
% variables directly.
%================================================================================================%
%% Initialization
diff_jac = sqrt(eps);                                                           % Extremely small value utilized for the finite difference
A = diag([1 -1 1 1 -1]);                                                        % Matrix required to enforce 2-step periodic locomotion

%% Derive subsequent MS state after one step for optimal MS state and control input pair
des_steps = 2;
% Convert x-y-z slice states to the corresponding initial conditions for the Dual-SLIP model
[x2,x3,x4,x7,x8,x9] = xyz_to_ICs(real(x_0_star(1)),real(x_0_star(2)),real(x_0_star(3)),real(x_0_star(4)),real(x_0_star(5)),0); 
% Assume that the feet is at zero height with zero velocity at t=0
init_conds = [0,x2,x3,x4,0,0,x7,x8,x9,0];                                       % Augment the initial condition variables with zeros for the 3D case of the Model
x_slice_x0_u0_star = sim_ms_to_ms_compliant(init_conds,u_0_star,des_steps);     % Simulate system for one step for the given optimal pair 
%================================================================================================%
%% Apply the small perturbation to the optimal MS state and control input vector
u_0_star_pert = u_0_star+diff_jac;
x_0_star_pert = x_0_star+diff_jac;
%================================================================================================%
%% Calculate the Jacobian with respect to the optimal MS state
df_dx = [];                                                                     % Initialize an empty vector to store the elements of the Jacobian
for i = 1:5                                                                     % Repeat the same process for all MS states
    temp = x_0_star;                                                            % Temporary variable that stores the optimal MS state
    temp(i) = x_0_star_pert(i);                                                 % Apply the perturbation only to one of the optimal MS states
    des_steps = 2;
    % Convert x-y-z slice states to the corresponding initial conditions for the Dual-SLIP model
    [x2,x3,x4,x7,x8,x9] = xyz_to_ICs(real(temp(1)),real(temp(2)),real(temp(3)),real(temp(4)),real(temp(5)),0); 
    % Assume that the feet is at zero height with zero velocity at t=0
    init_conds = [0,x2,x3,x4,0,0,x7,x8,x9,0];                                   % Augment the initial condition variables with zeros for the 3D case of the Model

    x_slice_Jx = sim_ms_to_ms_compliant(init_conds,u_0_star,des_steps);         % Simulate system for one step for the perturbed optimal pair
    df_dx = [df_dx (x_slice_Jx - x_slice_x0_u0_star)'/diff_jac];                % Derive and store finite difference for this perturbed MS state
end
%================================================================================================%
%% Calculate the Jacobian with respect to the optimal control input
du_dx = [];                                                                     % Initialize an empty vector to store the elements of the Jacobian
for i = 1:3                                                                     % Repeat the same process for all control inputs
    temp = u_0_star;                                                            % Temporary variable that stores the optimal control input
    temp(i) = u_0_star_pert(i);                                                 % Apply the perturbation only to one of the optimal control inputs
    des_steps = 2;
    % Convert x-y-z slice states to the corresponding initial conditions for the Dual-SLIP model
    [x2,x3,x4,x7,x8,x9] = xyz_to_ICs(real(x_0_star(1)),real(x_0_star(2)),real(x_0_star(3)),real(x_0_star(4)),real(x_0_star(5)),0); 
    % Assume that the feet is at zero height with zero velocity at t=0
    init_conds = [0,x2,x3,x4,0,0,x7,x8,x9,0];                                   % Augment the initial condition variables with zeros for the 3D case of the Model

    x_slice_Ju = sim_ms_to_ms_compliant(init_conds,temp,des_steps);             % Simulate system for one step for the perturbed optimal pair
    du_dx = [du_dx (x_slice_Ju - x_slice_x0_u0_star)'/diff_jac];                % Derive and store finite difference for this perturbed control input
end
%================================================================================================%
%% Derive the desired quantities for the implementation of the LQR controller
Jx = A*df_dx;
Ju = A*du_dx;
end