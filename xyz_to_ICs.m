function [x2,x3,x4,x7,x8,x9] = xyz_to_ICs(x_0_d,y_0_d,z_0,x_dot_0_d,y_dot_0_d,z_dot_0)
% function [x2,x3,x4,x7,x8,x9] = xyz_to_ICs(x_0_d,y_0_d,z_0,x_dot_0_d,y_dot_0_d,z_dot_0)

% Description:
% This function converts the MS state to the full state that is being
% utilized for the simulation of the system.

% INPUTS:
% x_0_d:            (m) Desired relative position in the x-axis of the point mass wrt to foot in support at MS
% y_0_d:            (m) Desired relative position in the y-axis of the point mass wrt to foot in support at MS
% z_0:              (m) Initial vertical position of the point mass
% x_dot_0_d:        (m/s) Initial forward velocity x_dot_{0,d} approximately between 0.7 m/s and 1.3 m/s.
% y_dot_0_d:        (m/s) Desired relative velocity in the y-axis at MS (Zero value is needed to satisfy the periodic gait conditions (Liu et al. (2015))
% z_dot_0           (m/s) Initial vertical velocity of the point mass

% OUTPUTS:
% x2:               Length of the foot for leg A
% x3:               Angle between the ground projection of the point mass and the x-axis
% x4:               Angle between the spring leg in support and the the xy plane
% x7:               Rate of change of the length of the foot for the in support
% x8:               Angular velocity between the ground projection of the point mass and the x-axis
% x9:               Angular velocity between the spring leg in support and the the xy plane

% Note:
% z_dot_0 is always 0.
%================================================================================================%

%Eq.83 in Report by Karakasis
x2 = sqrt( x_0_d^2 + y_0_d^2 + z_0^2 ); 
%Eq.84 in Report by Karakasis
x3 = atan2( y_0_d,x_0_d ); 
%Eq.85 in Report by Karakasis
x4 = atan2( z_0,sqrt(x_0_d^2 + y_0_d^2) );
%Eq.86 in Report by Karakasis
x7 = (x_0_d*x_dot_0_d + y_0_d*y_dot_0_d + z_0*z_dot_0) / (x2);
%Eq.87 in Report by Karakasis
x9 = (z_dot_0 - x7*sin(x4)) / (x2*cos(x4));
%Eq.88 in Report by Karakasis
x8 = (  -x_dot_0_d + x7*cos(x3)*cos(x4) -x2*cos(x3)*sin(x4)*x9 ) / (x2*sin(x3)*cos(x4));
end