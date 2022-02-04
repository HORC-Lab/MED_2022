In this project, the 3D Dual-SLIP model is simulated for a set of optimal initial conditions found using the nonlinear optimizer by Liu et al. (2015), but for compliant terrain.
The system is initiated using an optimal set of initial conditions found for the ideal rigid terrain, using the nonlinear optimizer by Liu et al. (2015).

Initially, the ground stiffness was set to the rigid value of 50 MN/m. Then, after 10 steps the system experienced a one-step unilateral lower stiffness perturbation, after which the ground stiffness was set back to rigid.
Specifically, at the 10th Touchdown event the ground stiffness under the leg about to land was lowered to a specific value and was kept constant throughout the whole step, while the ground stiffness of the leg in support remained to rigid. Next, after the completion of that step, the ground stiffness was set to rigid for both legs. In order to handle that perturbation, a
biomechanics-inspired controller is proposed. Specifically, at the 10th Touchdown, the stiffness of the leg in support is preadjusted to a higher
value, while the stiffness of the leg experiencing the perturbation is
also increased throughout the step. At the 11th Midstance step, the
stiffness of the leg about to land on rigid terrain is again set to
default.