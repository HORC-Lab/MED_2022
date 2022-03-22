function [x_slice,next,init_conds,x_A_f,y_A_f] = sim_td_ms_compliant(z_th,theta,phi,k_A,k_B,init_conds,x_A_f,y_A_f,steps,step_disturb)
% function [x_slice,next,init_conds,x_A_f,y_A_f] = sim_td_ms_compliant(z_th,theta,phi,k_A,k_B,init_conds,x_A_f,y_A_f,steps,step_disturb)

% Description:
% This function simulates the system from TD to MS, using the given
% configuration and control parameters. Initially, the system is being
% simulating using Single Support dynamics until the TD gait event takes
% place. Then, at TD, the first part of the proposed biomechanic-inspired
% controller is implemented. After that, Double Support dynamics are
% utilized until the LH and LO take place. Lastly, Single Support dynamics
% are adopted and the system is simulated until the MS gait event takes
% place. At that point, the MS state is evaluated and the function returns
% the value back to the main function, where the LQR and the second part of
% the proposed controller are being implemented.
% For information, the reader is referred to the technical report by
% Karakasis: "Report for 3D Dual SLIP Model Analysis".

% INPUTS:
% z_th:         (m) Threshold height at which touchdown takes place (determined by angle of attack) %See Eq.(163)
% theta:        (rad) Forward touchdown angle
% phi:          (rad) Lateral touchdown angle
% k_A:          (N/m) Spring stiffness of leg A
% k_B:          (N/m) Spring stiffness of leg B
% init_conds:   Initial conditions to be utilized for the next simulation
% x_A_f:        (m) Current position of supporting leg's foot on the x-axis
% y_A_f:        (m) Current position of supporting leg's foot on the y-axis
% steps:        Number of successful physical steps using the model
% step_disturb: Step at which the system will experience the stiffness perturbation

% OUTPUTS:
% x_slice:      The MS state derived at the last step of the simulation
% next:         Boolean variable used to indicate that the system failed before reaching the desired number of steps
% init_conds:   Initial conditions to be utilized for the next simulation
% x_A_f:        (m) Current position of supporting leg's foot on the x-axis
% y_A_f:        (m) Current position of supporting leg's foot on the y-axis

% NOTE:
%
%================================================================================================%
% The global variables defined in the main function
% "three_D_Dual_SLIP_compliant_sim_MED_2022" have to be declared here, so
% that they can be accessed in this function's workspace.
global m l_0 g m_f_A m_f_B k_g_A k_g_B b_g_A c_a b_g_B n_pow k_A k_B thld_vel thld_pos thld_l opts max_step_size;
%-------------------------------------------------------------------------%
%% Initializations
t_total=[];                                                         % In this vector we will store all time variables that the ode113 solver returns
sol=[];                                                             % In this vector we will store all state variables that the ode113 solver returns
i=1;                                                                % Iteration counter for the outer loop
i_i = 1;                                                            % Iteration counter for the inner loop
i_th = 1100;                                                        % Maximum number of inner loop iterations
size_thld = 300;                                                    % Maximum allowed size of time points in the 0.01s interval during one simulation (ode113)
%================================================================================================%
%% Analysis of the system equations that will be simulated using the "ode" step integrator during ISS for compliant terrain
syms x(t) [10 1] real;                                              % Define the full state of the system
%----- Equations for the Single Support Phase -----%
F_g_A = k_g_A*(-x1(t))^n_pow - b_g_A*x6(t)*(-x1(t))^n_pow;          % Interaction force to capture the compliance of the surface based on the Hunt-Crossley (HC) model
%Eq.137 in Report by Karakasis
eqn1 = diff(x1(t),t) == x6(t);
%Eq.138 in Report by Karakasis
eqn2 = diff(x2(t),t) == x7(t);
%Eq.139 in Report by Karakasis
eqn3 = diff(x3(t),t) == x8(t);
%Eq.140 in Report by Karakasis
eqn4 = diff(x4(t),t) == x9(t);
%Eq.141 in Report by Karakasis
eqn5 = diff(x5(t),t) == x10(t);
%Eq.142 in Report by Karakasis
eqn6 = (m+m_f_A)*diff(x6(t),t) == -(m+m_f_A)*g +F_g_A - m*diff(x7(t),t)*sin(x4(t)) + m*x2(t)*sin(x4(t))*(x9(t))^2 - m*x2(t)*diff(x9(t),t)*cos(x4(t)) - 2*m*x7(t)*x9(t)*cos(x4(t));
%Eq.143 in Report by Karakasis
eqn7 = m*diff(x7(t),t) + m*diff(x6(t),t)*sin(x4(t)) == k_A*(l_0 - x2(t)) - m*g*sin(x4(t)) + m*x2(t)*(x9(t))^2 + m*x2(t)*( x8(t)*cos(x4(t)) )^2;
%Eq.144 in Report by Karakasis
eqn8 = m*diff(x8(t),t)*(cos(x4(t))*x2(t))^2 == 2*m*(x2(t)^2)*x8(t)*x9(t)*cos(x4(t))*sin(x4(t)) - 2*m*x2(t)*x7(t)*x8(t)*(cos(x4(t)))^2;
%Eq.145 in Report by Karakasis
eqn9 = m*diff(x9(t),t)*(x2(t)^2) + m*x2(t)*diff(x6(t),t)*cos(x4(t)) == -m*g*x2(t)*cos(x4(t)) - 2*m*x2(t)*x7(t)*x9(t) - m*((x2(t)*x8(t))^2)*sin(x4(t))*cos(x4(t));
%Eq.146 in Report by Karakasis
eqn10 = diff(x10(t),t) == 0;
%----- Equations for the Single Support Phase -----%
%-------------------------------------------------------------------------%
eqns = [eqn1, eqn2, eqn3, eqn4, eqn5, eqn6, eqn7, eqn8, eqn9, eqn10];
vars = [x1(t) x2(t) x3(t) x4(t) x5(t) x6(t) x7(t) x8(t) x9(t) x10(t)];
[M,F] = massMatrixForm(eqns,vars);                                  % Extract mass matrix and right side of semilinear system of differential algebraic equations                      
f = M\F;                                                            % System of algebraic expressions
odefun = odeFunction(f,vars);                                       % Convert symbolic expressions to function handle for ODE solvers
%================================================================================================%
%================================================================================================%
TD_frames = [];                                                     %In this vector we will store all time instances of every Touchdown event encounter during the following iteration
i_i=1;                                                              % Initialization of the inner iteration counter
while (1)                                                           % Inner loop for ISS phase until TD
    [t_i,sol_i] = ode113(odefun,[0 max_step_size],init_conds,opts); % Solve the system's equations for the ISS Phase using the specified ICs for 0.01 secs with the error tolerances "opts" specified earlier
    if i == 1                                                       % For the first iteration of the loop
        t_total=[t_total;t_i];                                      % Adjust the local time vector, so that it corresponds to the total time elapsed and append it to the total time vector
        sol=[sol;real(sol_i)];                                      % Append the local solutions for the state to the total state solutions vector
    else
        t_total=[t_total(1:end-1);t_i+t_total(end)];                % Adjust the local time vector, so that it corresponds to the total time elapsed and append it to the total time vector
        sol=[sol(1:end-1,:);real(sol_i)];                           % Append the local solutions for the state to the total state solutions vector
    end
    init_conds=[sol(end,1); sol(end,2); sol(end,3); sol(end,4); sol(end,5); sol(end,6); sol(end,7); sol(end,8); sol(end,9); sol(end,10)];    %Set the last state of the system for this iteration, as ICs for the next one
    
    %Compute the vertical position and velocity of the CoM
    z_c = sol_i(:,1) + sol_i(:,2).*sin(sol_i(:,4));
    z_c_dot = sol_i(:,6) + sol_i(:,7).*sin(sol_i(:,4))+sol_i(:,2).*sol_i(:,9).*cos(sol_i(:,4));
    
    TD_frame=find( (z_c_dot<0) & (abs(z_c-z_th)<thld_pos) & (sol_i(:,2)<l_0+thld_l) );  % Check whether the condition for the Touchdown event is satisfied
    TD_frames= [TD_frames; TD_frame+length(sol)-length(sol_i)];     % Store all encountered frames of events wrt the total time vector, in order to find afterwards the most accurate

    i=i+1;                                                          % Increase outer loop counter
    i_i = i_i+1;                                                    % Increase inner loop counter
        
    % Check whether some some frames that satisfy the Touchdown gait
    % event condition have been found. Then, check whether any
    % candidate TD frames were found in this iteration. As long as more
    % candidates are found, the iteration will keep going on. The
    % iteartion, will stop when none TD candidates where found
    % in this iteration, but some were found in the previous ones.
    % Moreover, the iteration will stop immediately if the length of
    % the leg in support exceeds the natural length beyond the allowed
    % threshold, or if the maxinum number of iterations was reached.
    % Moreover, if the size of the vector t_i is greater than the threshold
    % value "size_thld", the simulation should stop, as this indicates that
    % the ode solver takes too much time to run.
    if ( (~isempty(TD_frames) && isempty(TD_frame) ) || ~isempty((find(sol(:,2)>l_0+thld_l))) || i_i>i_th || size(t_i,1)>size_thld)
        break;
    end
end
%-------------------------------------------------------------------------%
if(isempty(TD_frames))                                              % Check whether the condition for the Touchdown event was satisfied for any time instant during the iteration
    disp('No Touchdown event occured during the simulation. Please pick different initial conditions.')
    next = 1;                                                       % Update boolean variable to indicate that the simulation should stop
    x_slice = [];                                                   % Return an empty array to indicate that the system failed
    init_conds = [];                                                % Return an empty array to indicate that the system failed
    x_B_f = 0;                                                      % Return a zero value to indicate that the system failed
    y_B_f = 0;                                                      % Return a zero value to indicate that the system failed
    return;                                                         % In case no Touchdown events took place the simulation is terminated
end

%Compute the vertical position and velocity of the CoM during the total
%duration of the simulation so far, to identify the most precise frames of the following events
z_c = sol(:,1) + sol(:,2).*sin(sol(:,4));
z_c_dot = sol(:,6) + sol(:,7).*sin(sol(:,4))+sol(:,2).*sol(:,9).*cos(sol(:,4));
%The above variables were determined for the whole duration of the simulation
%in order to find directly the time instant of the most accurate event wrt to the total time
TD_frame = TD_frames(find(abs(z_c(TD_frames)-z_th)==min(abs(z_c(TD_frames)-z_th)),1,'first'));        % Determine which frame amongs the Touchdown frames is the most accurate and define it as the Touchdown frame
%-------------------------------------------------------------------------%
if ~isempty(find( (sol(1:TD_frame,2)-l_0) >thld_l ))                % Check if the length of the leg in support exceeds the natural length beyond the allowed threshold,
    next = 1;                                                       % Update boolean variable to indicate that the simulation should stop
    x_slice = [];                                                   % Return an empty array to indicate that the system failed
    init_conds = [];                                                % Return an empty array to indicate that the system failed
    x_B_f = 0;                                                      % Return a zero value to indicate that the system failed
    y_B_f = 0;                                                      % Return a zero value to indicate that the system failed
    return;                                                         % The simulation is terminated
end
%-------------------------------------------------------------------------%
%Keep only data up until the Touchdown Event and disregard any data acquired
%using the SS dynamics after Touchdown.
sol=sol(1:TD_frame,:);
t_total=t_total(1:TD_frame,:);
%Set as initial conditions the state at Touchdown, so that the Double
%Support dynamics can begin their implementation exactly at Touchdown
init_conds=[sol(TD_frame,1); sol(TD_frame,2); sol(TD_frame,3); sol(TD_frame,4); sol(TD_frame,5); sol(TD_frame,6); sol(TD_frame,7); sol(TD_frame,8); sol(TD_frame,9); sol(TD_frame,10)];
z_c_td = z_c(TD_frame);                                             % Height of the COM when the TD took place
%-------------------------------------------------------------------------%
%Determine the position of the point mass at Touchdown
x_c_TD = x_A_f+sol(TD_frame,2)*cos( sol(TD_frame,3) )*cos( sol(TD_frame,4) );
y_c_TD = y_A_f+sol(TD_frame,2)*sin( sol(TD_frame,3) )*cos( sol(TD_frame,4) );

% Determine the length of the second leg at Touchdown, as it cannot be
% assumed that it is equal to the rest length, due to precision errors
initial_leg_B_length = ( sol(TD_frame,1) + sol(TD_frame,2)*sin(sol(TD_frame,4)) - sol(TD_frame,5) )/sin(theta);
% Determine the x-y location of the second leg that just landed on the
% ground, using the provided control variables
x_B_f = x_c_TD-initial_leg_B_length*cos(theta)*cos(phi);
y_B_f = y_c_TD-initial_leg_B_length*cos(theta)*sin(phi);
%================================================================================================%
%================================================================================================%
%% Apply the first part of the proposed biomechanics-inspired controller
% If this is the Touchdown event of the step when the perturbation should 
% take place, then the stiffness of the leg in support is preadjusted to a
% higher value, while the stiffness of the leg about to step on the soft
% surface is also increased. If this is not the perturbation step, then the
% stiffness of both legs maintains that same value (k_A), as it was determined by
% the LQR controller at Midstance.
if steps == step_disturb                                            % Step where the system steps on the soft terrain
    k_B = k_A*7.452;                                                   % Stiffness of the second leg about to step on the soft terrain
    k_A = k_A*3.105;                                                  % Stiffness of the first leg in support currently walkin on hard terrain
end
%================================================================================================%
%================================================================================================%
%% Analysis of the system equations that will be simulated using the "ode" step integrator during DS for compliant terrain
%----- Equations for the Double Support Phase -----%
syms x(t) [10 1] real;                                              % Define the full state of the system
% Derive the length of the second leg as a function of the states
r_B = sqrt((x_A_f-x_B_f+x2(t)*cos(x3(t))*cos(x4(t)))^2 + (y_A_f-y_B_f+x2(t)*sin(x3(t))*cos(x4(t)))^2 + ( x1(t)-x5(t)+x2(t)*sin(x4(t)) )^2 );
F_g_A = k_g_A*(-x1(t))^n_pow - b_g_A*x6(t)*(-x1(t))^n_pow;          % Interaction force to capture the compliance of the surface under the leg A based on the Hunt-Crossley (HC) model
F_g_B = k_g_B*(-x5(t))^n_pow - b_g_B*x10(t)*(-x5(t))^n_pow;         % Interaction force to capture the compliance of the surface under the leg B based on the Hunt-Crossley (HC) model
%Eq.180 in Report by Karakasis
eqn1_ds = diff(x1(t),t) == x6;
%Eq.181 in Report by Karakasis
eqn2_ds = diff(x2(t),t) == x7;
%Eq.182 in Report by Karakasis
eqn3_ds = diff(x3(t),t) == x8;
%Eq.183 in Report by Karakasis
eqn4_ds = diff(x4(t),t) == x9;
%Eq.184 in Report by Karakasis
eqn5_ds = diff(x5(t),t) == x10;
%Eq.185 in Report by Karakasis
eqn6_ds = diff(x6(t),t)*( m+m_f_A ) + m*x2(t)*diff(x9(t),t)*cos(x4(t)) + m*diff(x7(t),t)*sin(x4(t)) == -( m+m_f_A )*g + F_g_A + (k_B*(l_0 - r_B)*( x1(t) - x5(t) + x2(t)*sin(x4(t)) )/r_B) + m*x2(t)*((x9(t))^2)*sin(x4(t)) - 2*m*x7(t)*x9(t)*cos(x4(t));
%Eq.186 in Report by Karakasis
eqn7_ds = m*( diff(x6(t),t)*sin(x4(t)) + diff(x7(t),t) ) == k_A*(l_0 - x2(t)) - m*g*sin(x4(t)) + m*x2(t)*( ((x9(t))^2) +((x8*cos(x4(t)))^2) ) + k_B*( (l_0-r_B)/r_B )*( (x_A_f-x_B_f)*cos(x3(t))*cos(x4(t)) + (y_A_f-y_B_f)*sin(x3(t))*cos(x4(t)) + (x1(t)-x5(t))*sin(x4(t)) + x2(t) );
%Eq.187 in Report by Karakasis
eqn8_ds = m*diff(x8(t),t)*(cos(x4(t))*x2(t))^2 == 2*m*(x2(t)^2)*x8(t)*x9(t)*cos(x4(t))*sin(x4(t)) - 2*m*x2(t)*x7(t)*x8(t)*(cos(x4(t)))^2 + k_B*((l_0-r_B)/r_B)*x2(t)*( -(x_A_f-x_B_f)*sin(x3(t))*cos(x4(t)) + (y_A_f-y_B_f)*cos(x3(t))*cos(x4(t)));
%Eq.188 in Report by Karakasis
eqn9_ds = m*x2(t)*( x2(t)*diff(x9(t),t) + x7(t)*x9(t) + diff(x6(t),t)*cos(x4(t))) == -m*x2(t)*x7(t)*x9(t) - m*g*x2(t)*cos(x4(t)) - m*((x2(t)*x8(t))^2)*cos(x4(t))*sin(x4(t)) - k_B*((l_0-r_B)/r_B)*x2(t)*( (x_A_f-x_B_f)*cos(x3(t))*sin(x4(t)) + (y_A_f-y_B_f)*sin(x3(t))*sin(x4(t)) - (x1(t)-x5(t))*cos(x4(t)) );
%Eq.189 in Report by Karakasis
eqn10_ds = m_f_B*diff(x10(t),t) == -m_f_B*g + F_g_B - k_B*( (l_0-r_B)/r_B )*( x1(t)-x5(t)+x2(t)*sin(x4(t)) );
%----- Equations for the Double Support Phase -----%
%-------------------------------------------------------------------------%
eqns_ds = [eqn1_ds, eqn2_ds, eqn3_ds, eqn4_ds, eqn5_ds, eqn6_ds, eqn7_ds, eqn8_ds, eqn9_ds, eqn10_ds];
vars_ds = [x1(t) x2(t) x3(t) x4(t) x5(t) x6(t) x7(t) x8(t) x9(t) x10(t)];
[M_ds,F_ds] = massMatrixForm(eqns_ds,vars_ds);                      % Extract mass matrix and right side of semilinear system of differential algebraic equations
f_ds = M_ds\F_ds;                                                   % System of algebraic expressions           
odefun = odeFunction(f_ds,vars_ds);                                 % Convert symbolic expressions to function handle for ODE solvers
%================================================================================================%
LH_frames=[];                                                       % In this vector we will store all time instances of every Lowest Height event encounter during the following iteration
i_i = 1;                                                            % Initialization of the inner iteration counter
while (1)                                                           % Inner loop for DS phase until LH
    [t_i,sol_i] = ode113(odefun,[0 max_step_size],init_conds,opts); % Solve the system's equations for the ISS Phase using the specified ICs for 0.001 secs with the error tolerances "opts" specified earlier
    t_total=[t_total(1:end-1);t_i+t_total(end)];                    % Adjust the local time vector, so that it corresponds to the total time elapsed and append it to the total time vector
    sol=[sol(1:end-1,:);real(sol_i)];                               % Append the local solutions for the state to the total state solutions vector
    init_conds=[sol(end,1); sol(end,2); sol(end,3); sol(end,4); sol(end,5); sol(end,6); sol(end,7); sol(end,8); sol(end,9); sol(end,10)];    %Set the last state of the system for this iteration, as ICs for the next one
        
    %Compute the vertical position and velocity of the CoM
    z_c = sol_i(:,1) + sol_i(:,2).*sin(sol_i(:,4));
    z_c_dot = sol_i(:,6) + sol_i(:,7).*sin(sol_i(:,4))+sol_i(:,2).*sol_i(:,9).*cos(sol_i(:,4));
    
    % Compute the lengths of the two legs
    lA = sol_i(:,2);
    lB = sqrt((x_A_f-x_B_f+sol_i(:,2).*cos(sol_i(:,3)).*cos(sol_i(:,4))).^2 + (y_A_f-y_B_f+sol_i(:,2).*sin(sol_i(:,3)).*cos(sol_i(:,4))).^2 + ( sol_i(:,1) - sol_i(:,5) + sol_i(:,2).*sin(sol_i(:,4)) ).^2 );
    
    % Check whether the condition for the Lowest Height event is satisfied
    LH_frame = find( (abs(z_c_dot)<thld_vel) & (z_c < z_c_td) & (lA < l_0+thld_l) & (lB < l_0+thld_l) );
    LH_frames = [LH_frames; LH_frame+length(sol)-length(sol_i)];    % Store all encountered frames of events wrt the total time vector, in order to find afterwards the most accurate

    i=i+1;                                                          % Increase outer loop counter
    i_i = i_i + 1;                                                  % Increase inner loop counter
        
    % Check whether some some frames that satisfy the Lowest Height gait
    % event condition have been found. Then, check whether any
    % candidate LH frames were found in this iteration. As long as more
    % candidates are found, the iteration will keep going on. The
    % iteartion, will stop when none LH candidates where found
    % in this iteration, but some were found in the previous ones.
    % Moreover, the iteration will stop immediately if the length of
    % the legs in support exceeds the natural length beyond the allowed
    % threshold, or if the maxinum number of iterations was reached.
    % Moreover, if the size of the vector t_i is greater than the threshold
    % value "size_thld", the simulation should stop, as this indicates that
    % the ode solver takes too much time to run.
    if ((~isempty(LH_frames) && isempty(LH_frame) ) || ~isempty((find(lA>l_0+thld_l))) || ~isempty((find(lB>l_0+thld_l))) || i_i>i_th || size(t_i,1)>size_thld)
        break;
    end
end
%-------------------------------------------------------------------------%
if(isempty(LH_frames))          % Check whether the condition for the Lowest Height event was satisfied for any time instant during the iteration
    disp('No Lowest Height event occured during the simulation. Please pick different initial conditions.');
    next = 1;                                                       % Update boolean variable to indicate that the simulation should stop
    x_slice = [];                                                   % Return an empty array to indicate that the system failed
    init_conds = [];                                                % Return an empty array to indicate that the system failed
    x_B_f = 0;                                                      % Return a zero value to indicate that the system failed
    y_B_f = 0;                                                      % Return a zero value to indicate that the system failed
    return;                                                         % In case no Lowest Height events took place the simulation is terminated
end

%Compute the vertical position and velocity of the CoM during the total
%duration of the simulation so far, to identify the most precise frames of the following events
z_c = sol(:,1) + sol(:,2).*sin(sol(:,4));
z_c_dot = sol(:,6) + sol(:,7).*sin(sol(:,4))+sol(:,2).*sol(:,9).*cos(sol(:,4));
%The above variables were determined for the whole duration of the simulation
%in order to find directly the time instant of the most accurate event wrt to the total time
LH_frame = LH_frames( find( (abs(z_c_dot(LH_frames))) == min((abs(z_c_dot(LH_frames)))),1,'first'));    %Determine which frame amongst the Lowest Height frames is the most accurate and define it as the Lowest Height frame
%-------------------------------------------------------------------------%
if ~isempty(find( (sol(1:LH_frame,2)-l_0) >thld_l ))            % Check if the length of the leg in support exceeds the natural length beyond the allowed threshold,
    next = 1;                                                       % Update boolean variable to indicate that the simulation should stop
    x_slice = [];                                                   % Return an empty array to indicate that the system failed
    init_conds = [];                                                % Return an empty array to indicate that the system failed
    x_B_f = 0;                                                      % Return a zero value to indicate that the system failed
    y_B_f = 0;                                                      % Return a zero value to indicate that the system failed
    return;                                                         % The simulation is terminated
end
%-------------------------------------------------------------------------%
%Keep only data up until the Lowest Height Event and disregard any data acquired
%using the DS dynamics after Lowest Height.
sol=sol(1:LH_frame,:);
t_total=t_total(1:LH_frame,:);
init_conds=[sol(LH_frame,1); sol(LH_frame,2); sol(LH_frame,3); sol(LH_frame,4); sol(LH_frame,5); sol(LH_frame,6); sol(LH_frame,7); sol(LH_frame,8); sol(LH_frame,9); sol(LH_frame,10)];    %Set the last state of the system for this iteration, as ICs for the next one
%================================================================================================%
%================================================================================================%
%% Analysis of the system equations that will be simulated using the "ode" step integrator during DS for compliant terrain
% Althought the same equations are utilized for both LH and LO events, they
% have to redefined to incorporate any changes to the utilized parameters,
% such as the stiffness of the legs. This is for future use.
%----- Equations for the Double Support Phase -----%
syms x(t) [10 1] real;                                              % Define the full state of the system
% Derive the length of the second leg as a function of the states
r_B = sqrt((x_A_f-x_B_f+x2(t)*cos(x3(t))*cos(x4(t)))^2 + (y_A_f-y_B_f+x2(t)*sin(x3(t))*cos(x4(t)))^2 + ( x1(t)-x5(t)+x2(t)*sin(x4(t)) )^2 );
F_g_A = k_g_A*(-x1(t))^n_pow - b_g_A*x6(t)*(-x1(t))^n_pow;          % Interaction force to capture the compliance of the surface under the leg A based on the Hunt-Crossley (HC) model
F_g_B = k_g_B*(-x5(t))^n_pow - b_g_B*x10(t)*(-x5(t))^n_pow;         % Interaction force to capture the compliance of the surface under the leg B based on the Hunt-Crossley (HC) model
%Eq.180 in Report by Karakasis
eqn1_ds = diff(x1(t),t) == x6;
%Eq.181 in Report by Karakasis
eqn2_ds = diff(x2(t),t) == x7;
%Eq.182 in Report by Karakasis
eqn3_ds = diff(x3(t),t) == x8;
%Eq.183 in Report by Karakasis
eqn4_ds = diff(x4(t),t) == x9;
%Eq.184 in Report by Karakasis
eqn5_ds = diff(x5(t),t) == x10;
%Eq.185 in Report by Karakasis
eqn6_ds = diff(x6(t),t)*( m+m_f_A ) + m*x2(t)*diff(x9(t),t)*cos(x4(t)) + m*diff(x7(t),t)*sin(x4(t)) == -( m+m_f_A )*g + F_g_A + (k_B*(l_0 - r_B)*( x1(t) - x5(t) + x2(t)*sin(x4(t)) )/r_B) + m*x2(t)*((x9(t))^2)*sin(x4(t)) - 2*m*x7(t)*x9(t)*cos(x4(t));
%Eq.186 in Report by Karakasis
eqn7_ds = m*( diff(x6(t),t)*sin(x4(t)) + diff(x7(t),t) ) == k_A*(l_0 - x2(t)) - m*g*sin(x4(t)) + m*x2(t)*( ((x9(t))^2) +((x8*cos(x4(t)))^2) ) + k_B*( (l_0-r_B)/r_B )*( (x_A_f-x_B_f)*cos(x3(t))*cos(x4(t)) + (y_A_f-y_B_f)*sin(x3(t))*cos(x4(t)) + (x1(t)-x5(t))*sin(x4(t)) + x2(t) );
%Eq.187 in Report by Karakasis
eqn8_ds = m*diff(x8(t),t)*(cos(x4(t))*x2(t))^2 == 2*m*(x2(t)^2)*x8(t)*x9(t)*cos(x4(t))*sin(x4(t)) - 2*m*x2(t)*x7(t)*x8(t)*(cos(x4(t)))^2 + k_B*((l_0-r_B)/r_B)*x2(t)*( -(x_A_f-x_B_f)*sin(x3(t))*cos(x4(t)) + (y_A_f-y_B_f)*cos(x3(t))*cos(x4(t)));
%Eq.188 in Report by Karakasis
eqn9_ds = m*x2(t)*( x2(t)*diff(x9(t),t) + x7(t)*x9(t) + diff(x6(t),t)*cos(x4(t))) == -m*x2(t)*x7(t)*x9(t) - m*g*x2(t)*cos(x4(t)) - m*((x2(t)*x8(t))^2)*cos(x4(t))*sin(x4(t)) - k_B*((l_0-r_B)/r_B)*x2(t)*( (x_A_f-x_B_f)*cos(x3(t))*sin(x4(t)) + (y_A_f-y_B_f)*sin(x3(t))*sin(x4(t)) - (x1(t)-x5(t))*cos(x4(t)) );
%Eq.189 in Report by Karakasis
eqn10_ds = m_f_B*diff(x10(t),t) == -m_f_B*g + F_g_B - k_B*( (l_0-r_B)/r_B )*( x1(t)-x5(t)+x2(t)*sin(x4(t)) );
%----- Equations for the Double Support Phase -----%
%-------------------------------------------------------------------------%
eqns_ds = [eqn1_ds, eqn2_ds, eqn3_ds, eqn4_ds, eqn5_ds, eqn6_ds, eqn7_ds, eqn8_ds, eqn9_ds, eqn10_ds];
vars_ds = [x1(t) x2(t) x3(t) x4(t) x5(t) x6(t) x7(t) x8(t) x9(t) x10(t)];
[M_ds,F_ds] = massMatrixForm(eqns_ds,vars_ds);                      % Extract mass matrix and right side of semilinear system of differential algebraic equations
f_ds = M_ds\F_ds;                                                   % System of algebraic expressions           
odefun = odeFunction(f_ds,vars_ds);                                 % Convert symbolic expressions to function handle for ODE solvers
%================================================================================================%
LO_frames=[];                                                       % In this vector we will store all time instances of every Lift Off event encounter during the following iteration
i_i = 1;                                                            % Initialization of the inner iteration counter
while (1)                                                           % Inner loop for DS phase until LO
    [t_i,sol_i] = ode113(odefun,[0 max_step_size],init_conds,opts);                                  %Solve the system's equations for the ISS Phase using the specified ICs for 0.001 secs with the error tolerances "opts" specified earlier
    t_total=[t_total(1:end-1);t_i+t_total(end)];                                                      %Adjust the local time vector, so that it corresponds to the total time elapsed and append it to the total time vector
    sol=[sol(1:end-1,:);real(sol_i)];                                                                        %Append the local solutions for the state to the total state solutions vector
    init_conds=[sol(end,1); sol(end,2); sol(end,3); sol(end,4); sol(end,5); sol(end,6); sol(end,7); sol(end,8); sol(end,9); sol(end,10)];    %Set the last state of the system for this iteration, as ICs for the next one
    
    %Compute the vertical position and velocity of the CoM
    z_c = sol_i(:,1) + sol_i(:,2).*sin(sol_i(:,4));
    z_c_dot = sol_i(:,6) + sol_i(:,7).*sin(sol_i(:,4))+sol_i(:,2).*sol_i(:,9).*cos(sol_i(:,4));

    % Compute the lengths of the two legs
    lA = sol_i(:,2);
    lB = sqrt((x_A_f-x_B_f+sol_i(:,2).*cos(sol_i(:,3)).*cos(sol_i(:,4))).^2 + (y_A_f-y_B_f+sol_i(:,2).*sin(sol_i(:,3)).*cos(sol_i(:,4))).^2 + ( sol_i(:,1) - sol_i(:,5) + sol_i(:,2).*sin(sol_i(:,4)) ).^2 );
    
    % Check whether the condition for the Lift Off event is satisfied
    LO_frame = find( (z_c_dot>0) & (abs(lA-l_0)<thld_l));
    LO_frames = [LO_frames; LO_frame+length(sol)-length(sol_i)];    % Store all encountered frames of events wrt the total time vector, in order to find afterwards the most accurate
    
    i=i+1;                                                          % Increase outer loop counter
    i_i = i_i + 1;                                                  % Increase inner loop counter

    % Check whether some some frames that satisfy the Lift Off gait
    % event condition have been found. Then, check whether any
    % candidate LO frames were found in this iteration. As long as more
    % candidates are found, the iteration will keep going on. The
    % iteartion, will stop when none LO candidates where found
    % in this iteration, but some were found in the previous ones.
    % Furthermore, the iteration will stop immediately if the length of
    % the legs in support exceeds the natural length beyond the allowed
    % threshold, or if the maxinum number of iterations was reached.
    % Moreover, if the size of the vector t_i is greater than the threshold
    % value "size_thld", the simulation should stop, as this indicates that
    % the ode solver takes too much time to run.
    % Finally, check whether the height of the foot in support has any 
    % positive values, meaning that it is floating on the air.
    if ((~isempty(LO_frames) && isempty(LO_frame)) || ~isempty((find(lA>l_0+thld_l))) || ~isempty((find(lB>l_0+thld_l))) || i_i>i_th || size(t_i,1)>size_thld || ~isempty(find(sol_i(:,1)>0)))
        break;
    end
end
%-------------------------------------------------------------------------%
% Keep only the LO candidate frames that correspond to an non-positive
% foot height.
LO_frames = intersect(LO_frames,find(sol(:,1)<=0));

if(isempty(LO_frames))                                              % Check whether the condition for the Lift Off event was satisfied for any time instant during the iteration
    disp('No Lift Off event occured during the simulation. Please pick different initial conditions.')
    next = 1;                                                       % Update boolean variable to indicate that the simulation should stop
    x_slice = [];                                                   % Return an empty array to indicate that the system failed
    init_conds = [];                                                % Return an empty array to indicate that the system failed
    x_B_f = 0;                                                      % Return a zero value to indicate that the system failed
    y_B_f = 0;                                                      % Return a zero value to indicate that the system failed
    return;                                                         % In case no Lift Off events took place the simulation is terminated
end

%Determine which frame amongst the Lift Off frames is the most accurate and define it as the Lift Off frame
LO_frame = LO_frames(find(abs(sol(LO_frames,2)-l_0)==min(abs(sol(LO_frames,2)-l_0)),1,'first'));
%-------------------------------------------------------------------------%
if ~isempty(find( (sol(1:LO_frame,2)-l_0) >thld_l ))                % Check if the length of the leg in support exceeds the natural length beyond the allowed threshold,
    next = 1;                                                       % Update boolean variable to indicate that the simulation should stop
    x_slice = [];                                                   % Return an empty array to indicate that the system failed
    init_conds = [];                                                % Return an empty array to indicate that the system failed
    x_B_f = 0;                                                      % Return a zero value to indicate that the system failed
    y_B_f = 0;                                                      % Return a zero value to indicate that the system failed
    return;                                                         % The simulation is terminated
end
%-------------------------------------------------------------------------%
%Keep only data up until the Lift Off Event and disregard any data acquired
%using the DS dynamics after Lift Off.
sol=sol(1:LO_frame,:);
t_total=t_total(1:LO_frame,:);
%================================================================================================%
%================================================================================================%
%% Modifications for SSS phase
%Switch back to Single Support Phase but using the other leg.
%Therefore, the state variables have to be redefined and then we need
%to express them as functions of the old state variables in order to
%derive new initial conditions.
%-------------------------------------------------------------------------%
%Intermediate variables
x_c = x_A_f + sol(:,2).*cos(sol(:,3)).*cos(sol(:,4));
x_c_dot = (sol(:,7).*cos(sol(:,3)).*cos(sol(:,4))-sol(:,2).*sol(:,8).*sin(sol(:,3)).*cos(sol(:,4))-sol(:,2).*sol(:,9).*cos(sol(:,3)).*sin(sol(:,4)) );

y_c = y_A_f + sol(:,2).*sin(sol(:,3)).*cos(sol(:,4));
y_c_dot = (sol(:,7).*sin(sol(:,3)).*cos(sol(:,4))+sol(:,2).*sol(:,8).*cos(sol(:,3)).*cos(sol(:,4))-sol(:,2).*sol(:,9).*sin(sol(:,3)).*sin(sol(:,4)) );

z_c = sol(:,2).*sin(sol(:,4))+sol(:,1);
z_c_dot = sol(:,6) + sol(:,7).*sin(sol(:,4))+sol(:,2).*sol(:,9).*cos(sol(:,4));
%x'_1 Eq.200 in Report by Karakasis
z_f_B = sol(:,5);
%x'_2 Eq.201 in Report by Karakasis
r_B_actual = sqrt((x_A_f-x_B_f+sol(:,2).*cos(sol(:,3)).*cos(sol(:,4))).^2 + (y_A_f-y_B_f+sol(:,2).*sin(sol(:,3)).*cos(sol(:,4))).^2 + ( sol(:,1)-sol(:,5)+sol(:,2).*sin(sol(:,4)) ).^2);
%x'_3 Eq.202 in Report by Karakasis
theta_1_B = atan2( (y_B_f-y_c),(x_B_f-x_c) );
%x'_4 Eq.203 in Report by Karakasis
theta_2_B = pi - asin( (z_c - sol(:,5))./( r_B_actual ) );
%x'_5 Eq.204 in Report by Karakasis
z_f_A = sol(:,1);
%x'_6 Eq.205 in Report by Karakasis
z_f_B_dot = sol(:,10);
%x'_7 Eq.206 in Report by Karakasis
rdba_1 = x_A_f-x_B_f+sol(:,2).*cos(sol(:,3)).*cos(sol(:,4));
rdba_2 = sol(:,2).*sol(:,8).*sin(sol(:,3)).*cos(sol(:,4)) - sol(:,7).*cos(sol(:,3)).*cos(sol(:,4)) + sol(:,2).*sol(:,9).*cos(sol(:,3)).*sin(sol(:,4));
rdba_3 = y_A_f-y_B_f+sol(:,2).*sin(sol(:,3)).*cos(sol(:,4));
rdba_4 = sol(:,7).*sin(sol(:,3)).*cos(sol(:,4)) + sol(:,2).*sol(:,8).*cos(sol(:,3)).*cos(sol(:,4)) - sol(:,2).*sol(:,9).*sin(sol(:,3)).*sin(sol(:,4));
rdba_5 = sol(:,1) - sol(:,5) + sol(:,2).*sin(sol(:,4));
rdba_6 = sol(:,7).*sin(sol(:,4)) + sol(:,6) - sol(:,10) + sol(:,2).*sol(:,9).*cos(sol(:,4));
r_dot_B_actual = ( (rdba_3.*rdba_4) - (rdba_1.*rdba_2) + (rdba_5.*rdba_6) )./r_B_actual;
%-------------------------------------------------------------------------%
%x'_8 Eq.212 in Report by Karakasis
theta_1_B_dot = ( -y_c_dot.*(x_B_f-x_c)+x_c_dot.*(y_B_f-y_c) ).*( (cos(theta_1_B)).^2 )./( (x_B_f-x_c).^2 );
%x'_9 Eq.214 in Report by Karakasis
theta_2_B_dot = -( (sol(:,6) + (sol(:,7).*sin(sol(:,4))) + (sol(:,2).*sol(:,9).*cos(sol(:,4))) - sol(:,10) ).*r_B_actual - ( sol(:,1) + sol(:,2).*sin(sol(:,4)) - sol(:,5) ).*r_dot_B_actual )./( (r_B_actual.^2).*cos(pi - theta_2_B) );
%x'_10 Eq.215 in Report by Karakasis
z_f_A_dot = sol(:,6);
%-------------------------------------------------------------------------%
%Set as initial conditions the state at Lift Off, so that the Secondary
%Initial Support dynamics can begin their %implementation exactly at Lift
%Off.

%Set the vertical position and velocity of the first leg (that now
%becomes second) to zero
init_conds=[z_f_B(LO_frame); r_B_actual(LO_frame); theta_1_B(LO_frame) ; theta_2_B(LO_frame); 0; z_f_B_dot(LO_frame); r_dot_B_actual(LO_frame); theta_1_B_dot(LO_frame); theta_2_B_dot(LO_frame); 0];
%================================================================================================%
%================================================================================================%
syms x(t) [10 1] real;                                              % Define the full state of the system
%----- Equations for the Single Support Phase -----%
F_g_B = k_g_B*(-x1(t))^n_pow - b_g_B*x6(t)*(-x1(t))^n_pow;          % Interaction force to capture the compliance of the surface based on the Hunt-Crossley (HC) model
%Eq.216 in Report by Karakasis
eqn1 = diff(x1(t),t) == x6(t);
%Eq.217 in Report by Karakasis
eqn2 = diff(x2(t),t) == x7(t);
%Eq.218 in Report by Karakasis
eqn3 = diff(x3(t),t) == x8(t);
%Eq.219 in Report by Karakasis
eqn4 = diff(x4(t),t) == x9(t);
%Eq.220 in Report by Karakasis
eqn5 = diff(x5(t),t) == x10(t);
%Eq.221 in Report by Karakasis
eqn6 = (m+m_f_B)*diff(x6(t),t) == -(m+m_f_B)*g +F_g_B - m*diff(x7(t),t)*sin(x4(t)) + m*x2(t)*sin(x4(t))*(x9(t))^2 - m*x2(t)*diff(x9(t),t)*cos(x4(t)) - 2*m*x7(t)*x9(t)*cos(x4(t));
%Eq.222 in Report by Karakasis
eqn7 = m*diff(x7(t),t) + m*diff(x6(t),t)*sin(x4(t)) == k_B*(l_0 - x2(t)) - m*g*sin(x4(t)) + m*x2(t)*(x9(t))^2 + m*x2(t)*( x8(t)*cos(x4(t)) )^2;
%Eq.223 in Report by Karakasis
eqn8 = m*diff(x8(t),t)*(cos(x4(t))*x2(t))^2 == 2*m*(x2(t)^2)*x8(t)*x9(t)*cos(x4(t))*sin(x4(t)) - 2*m*x2(t)*x7(t)*x8(t)*(cos(x4(t)))^2;
%Eq.224 in Report by Karakasis
eqn9 = m*diff(x9(t),t)*(x2(t)^2) + m*x2(t)*diff(x6(t),t)*cos(x4(t)) == -m*g*x2(t)*cos(x4(t)) - 2*m*x2(t)*x7(t)*x9(t) - m*((x2(t)*x8(t))^2)*sin(x4(t))*cos(x4(t));
%Eq.225 in Report by Karakasis
eqn10 = diff(x10(t),t) == 0;
%----- Equations for the Single Support Phase -----%
%-------------------------------------------------------------------------%
eqns = [eqn1, eqn2, eqn3, eqn4, eqn5, eqn6, eqn7, eqn8, eqn9, eqn10];
vars = [x1(t) x2(t) x3(t) x4(t) x5(t) x6(t) x7(t) x8(t) x9(t) x10(t)];
[M,F] = massMatrixForm(eqns,vars);                                  % Extract mass matrix and right side of semilinear system of differential algebraic equations
f = M\F;                                                            % System of algebraic expressions
odefun = odeFunction(f,vars);                                       % Convert symbolic expressions to function handle for ODE solvers
%-------------------------------------------------------------------------%
% Now the second leg is considered as the leg in support. Switch the
% foot position of the leg A
x_A_f = x_B_f;
y_A_f = y_B_f;
%-------------------------------------------------------------------------%
% Switch ground stiffness and damping between feet.
% Now the second leg will be in support. Although the equations from LO
% to MS are defined above using the details of the second leg (subscript b)
% , after the MS, the equations utilize the details of the first leg
% (subscript a). Hence, it needs to be updated.
k_g_A = k_g_B;
b_g_A = 1.5*c_a*k_g_A;
k_A = k_B;
%-------------------------------------------------------------------------%
MS_frames = [];                                                     % In this vector we will store all time instances of every Midstance event encounter during the following iteration
i_i=1;                                                              % Initialization of the inner iteration counter
while (1)                                                           % Inner loop for ISS phase until MS
    [t_i,sol_i] = ode113(odefun,[0 max_step_size],init_conds,opts); % Solve the system's equations for the ISS Phase using the specified ICs for 0.001 secs with the error tolerances "opts" specified earlier
    t_total=[t_total(1:end-1);t_i+t_total(end)];                    % Adjust the local time vector, so that it corresponds to the total time elapsed and append it to the total time vector
    sol=[sol(1:end-1,:);real(sol_i)];                               % Append the local solutions for the state to the total state solutions vector
    init_conds=[sol(end,1); sol(end,2); sol(end,3); sol(end,4); sol(end,5); sol(end,6); sol(end,7); sol(end,8); sol(end,9); sol(end,10)];    %Set the last state of the system for this iteration, as ICs for the next one
    
    % Compute the vertical position and velocity of the CoM
    z_c = sol_i(:,1) + sol_i(:,2).*sin(sol_i(:,4));
    z_c_dot = sol_i(:,6) + sol_i(:,7).*sin(sol_i(:,4))+sol_i(:,2).*sol_i(:,9).*cos(sol_i(:,4));
    
    % Check whether the condition for the Midstance event is satisfied
    MS_frame=find( (abs(z_c_dot)<thld_vel) & (z_c > z_th) & (sol_i(:,2)<l_0+thld_l) );%Default
    MS_frames= [MS_frames; MS_frame+length(sol)-length(sol_i)];     % Store all encountered frames of events wrt the total time vector, in order to find afterwards the most accurate
    
    i=i+1;                                                          % Increase outer loop counter
    i_i = i_i+1;                                                    % Increase inner loop counter

    % Check whether some some frames that satisfy the Midstance gait
    % event condition have been found. Then, check whether any
    % candidate MS frames were found in this iteration. As long as more
    % candidates are found, the iteration will keep going on. The
    % iteartion, will stop when none MS candidates where found
    % in this iteration, but some were found in the previous ones.
    % Moreover, the iteration will stop immediately if the length of
    % the leg in support exceeds the natural length beyond the allowed
    % threshold, or if the maxinum number of iterations was reached.
    % Moreover, if the size of the vector t_i is greater than the threshold
    % value "size_thld", the simulation should stop, as this indicates that
    % the ode solver takes too much time to run.
    if ( (~isempty(MS_frames) && isempty(MS_frame) ) || ~isempty((find(sol(:,2)>l_0+thld_l))) || ~isempty((find(sol(:,1)>0))) || i_i>i_th || size(t_i,1)>size_thld)
        break;
    end
end
MS_frames = intersect(MS_frames,find(sol(:,1)<=0)); 
%-------------------------------------------------------------------------%
if(isempty(MS_frames))                                              % Check whether the condition for the Midstance event was satisfied for any time instant during the iteration
    disp('No Midstance event occured during the simulation. Please pick different initial conditions.');
    next = 1;                                                       % Update boolean variable to indicate that the simulation should stop
    x_slice = [];                                                   % Return an empty array to indicate that the system failed
    init_conds = [];                                                % Return an empty array to indicate that the system failed
    x_B_f = 0;                                                      % Return a zero value to indicate that the system failed
    y_B_f = 0;                                                      % Return a zero value to indicate that the system failed
    return;                                                         % In case no Midstance events took place the simulation is terminated
end

%Compute the vertical position and velocity of the CoM during the total
%duration of the simulation so far, to identify the most precise frames of the following events
z_c = sol(:,1) + sol(:,2).*sin(sol(:,4));
z_c_dot = sol(:,6) + sol(:,7).*sin(sol(:,4))+sol(:,2).*sol(:,9).*cos(sol(:,4));
%The above variables were determined for the whole duration of the simulation
%in order to find directly the time instant of the most accurate event wrt to the total time
MS_frame = MS_frames(find(abs(z_c_dot(MS_frames))==min(abs(z_c_dot(MS_frames))),1,'first'));        %Determine which frame amongst the Midstance frames is the most accurate and define it as the Midstance frame
%-------------------------------------------------------------------------%
if ~isempty(find( (sol(1:MS_frame,2)-l_0) >thld_l ))    % Check if the length of the leg in support exceeds the natural length beyond the allowed threshold,
    next = 1;                                                       % Update boolean variable to indicate that the simulation should stop
    x_slice = [];                                                   % Return an empty array to indicate that the system failed
    init_conds = [];                                                % Return an empty array to indicate that the system failed
    x_B_f = 0;                                                      % Return a zero value to indicate that the system failed
    y_B_f = 0;                                                      % Return a zero value to indicate that the system failed
    return;                                                         % The simulation is terminated
end
%-------------------------------------------------------------------------%
%Keep only data up until the Midstance Event and disregard any data acquired
%after Midstance.
sol=sol(1:MS_frame,:);
t_total=t_total(1:MS_frame,:);
init_conds=[sol(MS_frame,1); sol(MS_frame,2); sol(MS_frame,3); sol(MS_frame,4); sol(MS_frame,5); sol(MS_frame,6); sol(MS_frame,7); sol(MS_frame,8); sol(MS_frame,9); sol(MS_frame,10)];
%-------------------------------------------------------------------------%
% Evaluate the slice state at MS
x_c_MS = x_A_f+sol(MS_frame,2)*cos( sol(MS_frame,3) )*cos( sol(MS_frame,4) );
y_c_MS = y_A_f+sol(MS_frame,2)*sin( sol(MS_frame,3) )*cos( sol(MS_frame,4) );
x_c_dot_MS = sol(MS_frame,7)*cos( sol(MS_frame,3) )*cos( sol(MS_frame,4) ) +...
    (-1)*sol(MS_frame,2)*sin( sol(MS_frame,3) )*sol(MS_frame,8)*cos( sol(MS_frame,4) ) +...
    (-1)*sol(MS_frame,2)*cos( sol(MS_frame,3) )*sin( sol(MS_frame,4) )*sol(MS_frame,9);
y_c_dot_MS = sol(MS_frame,7)*sin( sol(MS_frame,3) )*cos( sol(MS_frame,4) ) +...
    sol(MS_frame,2)*cos( sol(MS_frame,3) )*sol(MS_frame,8)*cos( sol(MS_frame,4) ) +...
    (-1)*sol(MS_frame,2)*sin( sol(MS_frame,3) )*sin( sol(MS_frame,4) )*sol(MS_frame,9);
x_slice = [x_c_MS-x_A_f y_c_MS-y_A_f z_c(MS_frame) x_c_dot_MS y_c_dot_MS]';
next = 0;                                                           % Indicate that the simulatin from TD to MS was completed successfully
end