%% "three_D_Dual_SLIP_compliant_sim_MED_2022.m"
% Author: Chrysostomos Karakasis

% Description:
% In this code, the 3D Dual-SLIP model is simulated for a set of optimal initial
% conditions found using the nonlinear optimizer by Liu et al. (2015), but for compliant terrain.
% The system is initiated using an optimal set of initial conditions found
% for the ideal rigid terrain, using the nonlinear optimizer by Liu et al. (2015).
% Initially, the ground stiffness was set to the rigid value of 50 MN/m.
% Then, after 10 steps the system experienced a one-step unilateral lower
% stiffness perturbation, after which the ground stiffness was set back to rigid.
% Specifically, at the 10th Touchdown event the ground stiffness under the leg
% about to land was lowered to a specific value and was kept constant
% throughout the whole step, while the ground stiffness of the leg in support
% remained to rigid. Next, after the completion of that step, the ground stiffness
% was set to rigid for both legs. In order to handle that perturbation, a
% biomechanics-inspired controller is proposed. Specifically, at the 10th
% Touchdown, the stiffness of the leg in support is preadjusted to a higher
% value, while the stiffness of the leg experiencing the perturbation is
% also increased throughout the step. At the 11th Midstance step, the
% stiffness of the leg about to land on rigid terrain is again set to
% default.

% NOTE:
% For more information see report "Report for 3D Dual SLIP Model Analysis"

% Last Update: 1/30/2022

warning("off")
clear all;
close all;

% Definition of global variables that will be shared among functions
global m l_0 g m_f_A m_f_B k_g_A k_g_B b_g_A c_a b_g_B n_pow k_A k_B thld_vel thld_pos thld_l opts max_step_size;
%================================================================================================%
%%  Sec.(1) - The following values were utilized as proposed in "Dynamic walking..." DOI: 10.1109/ICRA.2015.7139999
m = 80;                                                         % (kg) The mass of the CoM
l_0 = 1;                                                        % (m) The rest length of both legs
g = 9.81;                                                       %( m/s) Free fall acceleration
%================================================================================================%
%% Sec. (2) - Define parameters for compliant terrain
m_f_A = 1;                                                      % (kg) The mass of the left foot (Set to zero for the rigid case)
m_f_B = 1;                                                      % (kg) The mass of the right foot (Set to zero for the rigid case)
k_g_A = 0.5*10^8;                                               % (N/m) The stiffness of the ground that the left foot feels (Set to 0.5*10^8 to simulate rigid case)
k_g_B = 0.5*10^8;                                               % (N/m) The stiffness of the ground that the right foot feels (Set to 0.5*10^8 to simulate rigid case)
c_a = 0.2;                                                      % (scalar) Usually between 0.01-0.5 depending on the materials and impact velocity (Check Vasilopoulos et al. (2014) DOI: 10.1109/IROS.2014.6943251)
b_g_A = 1.5*c_a*k_g_A;                                          % (N*s/m) The damping of the ground that the left foot feels (Function of the stiffness - Eq.(5) at Vasilopoulos et al. (2014))
b_g_B = 1.5*c_a*k_g_B;                                          % (N*s/m) The damping of the ground that the left foot feels (Function of the stiffness - Eq.(5) at Vasilopoulos et al. (2014))
n_pow = 1.5;                                                    % (scalar) In the case of Hertzian non-adhesive contact is equal to 1.5 (Vasilopoulos et al. (2014))
%------------------------------------------------------------------------------------------------%
% Assign a random Initial Foot Position of Leg A (leg in support at t=0)
x_A_f = 2;                                                      % (m) Initial position of supporting leg's foot on the x-axis
y_A_f = 1;                                                      % (m) Initial position of supporting leg's foot on the y-axis
z_A_f = 0;                                                      % (m) Initial vertical position of supporting leg's foot
foot_pos_A_B = [x_A_f y_A_f];                                   % Store initial foot position of Leg A (leg in support at t=0)
%------------------------------------------------------------------------------------------------%
des_steps = 100;                                                % Set desired number of steps for the simulation of the system
step_disturb = 10;                                              % (scalar) Step at which the system will experience the stiffness perturbation
stiff_disturb = 25000;                                          % (N/m) The value at which the ground stiffness will be set to during the perturbation
stiff_default = k_g_A;                                          % The default value at which the ground stiffness will be set to after the perturbation
%================================================================================================%
%% Sec.(3) - Defining Threshold values for Equality Conditions for Gait Event Detection
thld_vel = 0.01;                                                % Threshold for the detection of Midstance (MS) and Lowest Height (LH) gait events
thld_pos = 0.0005;                                              % Threshold for the detection of Touchdown (TD) gait events
thld_l = 0.007;                                                 % Threshold for the detection of Lift Off (LO) gait events
%================================================================================================%
%% Initialization of Parameters required for Simulation
opts = odeset('RelTol',1e-6,'AbsTol',1e-6);                     % Absolute and Relative error tolerances for ode113 (Based on values mentioned in Sec.4.2.4 in Hartmut Geyer PhD Thesis)
max_step_size = 0.01;                                           % Time duration of each step simulation
next = 0;                                                       %% Boolean variable used to indicate that the system failed before reaching the desired number of steps
%================================================================================================%
%% Optimal Pair of Initial MS State and Control Input Vector for the System
% Definition of Optimal Initial MS State
% For the implementation by Liu et al. (2015), we will use the slice MS state x = [x_rel y_rel z_c x_dot_rel y_dot_rel]
% x0 = [x_0_d y_0_d z_0 x_dot_0_d y_dot_0_d]
%-------------------------------------------------------------------------%
% Fixed parameters
x_0_d = 0;                                                      % (m) Desired relative position in the x-axis of the point mass wrt to foot in support at MS
y_0_d = 0.05;                                                   % (m) Desired relative position in the y-axis of the point mass wrt to foot in support at MS
y_dot_0_d = 0;                                                  % (m/s) Desired relative velocity in the y-axis at MS (Zero value is needed to satisfy the periodic gait conditions (Liu et al. (2015))
%-------------------------------------------------------------------------%
% Definition of optimal initial conditions as they were derived using the nonlinear optimizer by Liu et al. (2015)
% to initiate a periodic gait with a desired forward velocity at MS for the original 3D Dual-SLIP Model
x_dot_0_d = 1;                                                  % (m/s) Initial forward velocity x_dot_{0,d} approximately between 0.7 m/s and 1.3 m/s.
z_0 = 0.987530775175059;                                        % (m) Initial vertical position of the point mass
x_0_star = [x_0_d y_0_d z_0 x_dot_0_d y_dot_0_d]';              % Derive optimal initial MS state
z_dot_0 = 0;                                                    % System is always initiazed at MS
%------------------------------------------------------------------------------------------------%
% Definition of Optimal Control Input
a_o_a = 72.7411042049505;                                       % (deg) Angle of Attack associated with the forward touchdown angle
theta = deg2rad(180-a_o_a);                                     % (rad) Forward touchdown angle
phi = 0.190971419146386;                                        % (rad) Lateral touchdown angle
k_A = 14163.5345312576;                                         % (N/m) Spring stiffness of both legs
k_B=k_A;                                                        % Set the same stiffness to both legs
u_0_star = [phi theta k_A]';                                    % Derive optimal control input vector
z_th = l_0*sin(theta);                                          % (m) Threshold height at which touchdown takes place (determined by angle of attack) %See Eq.(163)
%================================================================================================%
%% Convert slice MS state to the state proposed for our implementation
%Function that converts x-y-z slice states to the corresponding initial conditions for the Dual-SLIP model
[x2,x3,x4,x6,x7,x8] = xyz_to_ICs(x_0_star(1),x_0_star(2),x_0_star(3),x_0_star(4),x_0_star(5)+0.0,z_dot_0);
% This is where the user can also apply disturbance to any of the states at t=0

% Assume that the foot is at zero height with zero velocity at t=0
init_conds = [0,x2,x3,x4,0,0,x6,x7,x8,0];                          % Augment the initial condition variables with zeros for the 3D case of the Model
%-------------------------------------------------------------------------%
%The following equations are a reminder of the physical meaning of the state variables x1-x10
%See Section 2 in report "Report for 3D Dual SLIP Model Analysis"
% x1 = z_f_A                                                    % Vertical position of the foot for the leg in support
% x2 = r_A                                                      % Lenght of the foot for leg A
% x3 = theta1                                                   % Angle between the ground projection of the point mass and the x-axis
% x4 = theta2                                                   % Angle between the spring leg in support and the the xy plane
% x5 = z_f_B                                                    % Vertical position of the foot for the leg not in support
% x6 = dot{z}_f_A                                               % Vertical velocity of the foot for the leg in support
% x7 = dot{r}_A                                                 % Rate of change of the length of the foot for the in support
% x8 = dot{theta1}                                              % Angular velocity between the ground projection of the point mass and the x-axis
% x9 = dot{theta2}                                              % Angular velocity between the spring leg in support and the the xy plane
% x10 = dot{z}_f_B                                              % Vertical velocity of the foot for the leg not in support
%End of reminder
%================================================================================================%
%% Defition and Initialization regarding the Energy of the System
% Initial Kinetic Energy of the System (Eq. 124-125)
T_initial =  0.5*m*( (init_conds(7)+init_conds(6)*sin(init_conds(4)))^2 + (init_conds(2)*init_conds(8)*cos(init_conds(4)))^2 + (init_conds(2)*init_conds(9)+init_conds(6)*cos(init_conds(4)))^2 ) + 0.5*m_f_A*( (init_conds(6)*sin(init_conds(4)))^2 + (init_conds(6)*cos(init_conds(4)))^2 );
% Initial Potential Energy of the System  (Eq. 126)
V_initial = 0.5*k_A*(l_0 - init_conds(2))^2 + m*g*(init_conds(2)*sin(init_conds(4)) + init_conds(1)) + m_f_A*g*init_conds(1);
% Initial Total Energy of the System
E_initial = T_initial + V_initial;
%================================================================================================%
%% Approximate Jacobians evaluated at the optimal pair (x_0_star,u_0_star)
[Jx,Ju] = jac_comp(x_0_star,u_0_star);

%% Initializations required during the simulation
steps=0;                                                        % Number of successful physical steps using the model
Dxn_vector = [];                                                % Initialization of vector that will store the state errors across all steps
Dun_vector = [];                                                % Initialization of vector that will store the control input errors across all steps

%================================================================================================%
%================================================================================================%
%% Initiate simulation of the system using the initial conditions defined above
%-------------------------------------------------------------------------%
% First, use the following function to search for the first Midstance gait
% event, using the initial conditions and the optimal control inputs for
% the system. Avoiding this step causes the system to act weird.
x_slice = sim_ms_to_ms_compliant(init_conds,u_0_star,1)';
%-------------------------------------------------------------------------%
while(steps < des_steps)                                        % While loop executed for every step. "des_steps" specifies the desired number of steps
    
    steps = steps + 1;                                          % Increase the number of successful physical steps
    %================================================================================================%
    %% Apply LQR Controller
    q1 = 1; q2 = 1; q3 = 1; q4 = 1; q5 = 1;                     % Select the diagonal elements for the Q matrix of the LQR controller
    r1 = 1; r2 = 1; r3 = 1;                                     % Select the diagonal elements for the R matrix of the LQR controller
    % Define the diagonal Q-R matrices required for the LQR controller
    Q_dlqr = [q1 0 0 0 0; 0 q2 0 0 0; 0 0 q3 0 0; 0 0 0 q4 0; 0 0 0 0 q5];
    R_dlqr=[r1 0 0 ; 0 r2 0; 0 0 r3];
    %-------------------------------------------------------------------------%
    A = diag([1 -1 1 1 -1]);                                    % Matrix required to enforce 2-step periodic locomotion
    B = diag([-1 1 1]);                                         % Matrix required to account for the sign-alternating forward touchdown angle at each step
    x_slice_star_n = A^(steps-1)*x_0_star;                      % Calculate desired-optimal MS state at step n
    u_star_n = B^(steps-1)*u_0_star;                            % Calculate desired-optimal control input at step n
    Dxn = x_slice - x_slice_star_n;                             % Calculate error between actual and desired MS state at step n
    Dxn_vector = [Dxn_vector Dxn];                              % Store the MS state errors at every step for plotting purposes
    %-------------------------------------------------------------------------%
    A_dlqr = Jx;                                                % Utilize Jacobians evaluated at the optimal pair (x_0_star,u_0_star)
    B_dlqr = Ju;                                                % to linearize the return ma
    K = -dlqr(A_dlqr,B_dlqr,Q_dlqr,R_dlqr);                     % Derive the time-invariant, stable feedback gain from the LQR controller
    %-------------------------------------------------------------------------%
    eig(A_dlqr+B_dlqr*K);                                       % Verify that the feedback regulated system is stable
    rank(ctrb(A_dlqr,B_dlqr));                                  % Verify that the (Jx, Ju) pair is controllable
    %-------------------------------------------------------------------------%
    un = u_star_n + (B^(steps-1))*K*(A^(steps-1))*Dxn;          % Calculate the feedback law at the step n
    phi = un(1);                                                % Update the lateral touchdown angle
    theta = un(2);                                              % Update the forward touchdown angle
    z_th = l_0*sin(theta);                                      % Update Threshold height at which Touchdown takes place (determined by angle of attack)
    Dun_vector = [Dun_vector un-u_star_n];                      % Store the control input errors at every step for plotting purposes
    
    %================================================================================================%
    %% Apply second part of the biomechanics-inspired modified controller
    % The first part of the proposed controller is applied at Touchdown (Line )
    % Based on the proposed controller, we wish to maintain a high stiffness
    % gain for the leg experiencing the low stiffness during the step after
    % the perturbation.
    if steps == step_disturb + 1                                % Step where the leg on soft terrain in now in support
        %k_A now corresponds to the leg stiffness of the leg in support which is in on the soft terrain
        k_A = un(3)*11;                                         % Apply high stiffness gain to the latest output of the LQR gain
        k_B = un(3);                                            % Apply default stiffness for the second leg about to step on the rigid terrain
    else
        k_A = un(3);                                            % Change leg stiffness for both legs according to the LQR
        k_B = k_A;
    end
    %================================================================================================%
    %% Plotting Section
    % In this section, the MS state and control input errors are plotted
    % to observe the system response to the perturbation and to verify
    % the validity of the proposed controller.
    figure(1)                                                   % Plot the MS state error with respect to steps
    subplot(5,1,1)
    plot(Dxn_vector(1,:),'o-')                                  % Error for the first MS state with respect to steps
    grid on;
    ylabel("$\Delta x$ [m]",'interpreter','latex')
    subplot(5,1,2)
    plot(Dxn_vector(2,:),'o-')                                  % Error for the second MS state with respect to steps
    grid on;
    ylabel("$\Delta y$ [m]",'interpreter','latex')
    subplot(5,1,3)
    plot(Dxn_vector(3,:),'o-')                                  % Error for the third MS state with respect to steps
    grid on;
    ylabel("$\Delta z$ [m]",'interpreter','latex')
    subplot(5,1,4)
    plot(Dxn_vector(4,:),'o-')                                  % Error for the fourth MS state with respect to steps
    grid on;
    ylabel("$\Delta \dot{x}$ [m/s]",'interpreter','latex')
    subplot(5,1,5)
    plot(Dxn_vector(5,:),'o-')                                  % Error for the fifth MS state with respect to steps
    grid on;
    ylabel("$\Delta \dot{y}$ [m/s]",'interpreter','latex')
    xlabel("Step (Mid-stance) count")                           % X-axis represents the number of successful steps and the corresponding MS gait events
    sgtitle({['State Error at Midstance Events'],['Steps: ',num2str(steps),' - Max error: ',num2str(max(max(abs(Dxn_vector))))],['Feet Mass: ',num2str(m_f_A),' kg - Ground Stiffness: ',num2str(k_g_A/1000),' kN/m - Ground Damping: ',num2str(b_g_A/1000),' kNs/m'],['Stiffness Disturbance: ',num2str(stiff_disturb),' N/m - Step: ',num2str(step_disturb)]},'interpreter','latex','fontsize',18)
    %-------------------------------------------------------------------------%
    figure(2)                                                   % Plot the MS state error with respect to steps
    subplot(3,1,1)
    plot(Dun_vector(1,:),'o-')                                  % Error for the first control input with respect to steps
    grid on;
    ylabel("$\Delta \phi$ [rad]",'interpreter','latex')
    subplot(3,1,2)
    plot(Dun_vector(2,:),'o-')                                  % Error for the second control input with respect to steps
    grid on;
    ylabel("$\Delta \theta$ [rad]",'interpreter','latex')
    subplot(3,1,3)
    plot(Dun_vector(3,:),'o-')                                  % Error for the third control input with respect to steps
    grid on;
    ylabel("$\Delta k$ [N/m]",'interpreter','latex')
    xlabel("Step (Mid-stance) count")                           % X-axis represents the number of successful steps and the corresponding MS gait events
    sgtitle({['Inputs Error at Midstance Events'],['Steps: ',num2str(steps),' - Max error: ',num2str(max(max(abs(Dun_vector))))],['Feet Mass: ',num2str(m_f_A),' kg - Ground Stiffness: ',num2str(k_g_A/1000),' kN/m - Ground Damping: ',num2str(b_g_A/1000),' kNs/m'],['Stiffness Disturbance: ',num2str(stiff_disturb),' N/m - Step: ',num2str(step_disturb)]},'interpreter','latex','fontsize',18)
    %================================================================================================%
    %% Stiffness Perturbation
    % This will applied at Touchdown, but it being coded here to avoid
    % passing arguments to the function.
    if steps == step_disturb                                    % Apply stiffness perturbation ONLY on that step
        k_g_B = stiff_disturb;                                  % Set ground stiffness below the leg about to land to lower value specified in the beginning of the code
        b_g_B = 1.5*c_a*k_g_B;                                  % Update the damping of the surface, as it is a function of the stiffness
    else
        k_g_B = stiff_default;                                  % Reset ground stiffness back to default (rigid)
        b_g_B = 1.5*c_a*k_g_B;                                  % Update the damping of the surface, as it is a function of the stiffness
    end
    %================================================================================================%
    %% Dynamic Simulation of the system from Touchdown to Midstance
    % In this section, the system is simulated using a function to identify
    % the gait events Touchdown, Lowest Height, Lift Off and Midstance.
    [x_slice,next,init_conds,x_A_f,y_A_f] = sim_td_ms_compliant(z_th,theta,phi,k_A,k_B,init_conds,x_A_f,y_A_f,steps,step_disturb);
    % In case of an error the variable "next" will be set to "1" and the
    % simulation must stop.
    if next == 1
        break;
    end
end

