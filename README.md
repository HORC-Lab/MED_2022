 Abstract
Bipedal walking is one of the most important hallmarks of human that robots have been trying to mimic for many decades. Although previous control methodologies have achieved robot walking on some terrains, there is a need for a framework allowing stable and robust locomotion over a wide range of compliant surfaces. This work proposes a novel biomechanics-inspired controller that adjusts the stiffness of the legs in support for robust and dynamic bipedal locomotion over compliant terrains. First, the 3D Dual-SLIP model is extended to support for the first time locomotion over compliant surfaces with variable stiffness and damping parameters. Then, the proposed controller is compared to a Linear-Quadratic Regulator (LQR) controller, in terms of robustness on stepping on soft terrain. The LQR controller is shown to be robust only up to a moderate ground stiffness level of 174 kN/m, while it fails in lower stiffness levels. On the contrary, the proposed controller can produce stable gait in stiffness levels as low as 30 kN/m, which results in a vertical ground penetration of the leg that is deeper than 10% of its rest length. The proposed framework could advance the field of bipedal walking, by generating stable walking trajectories for a wide range of compliant terrains useful for the control of bipeds and humanoids, as well as by improving controllers for prosthetic devices with tunable stiffness.

# Description:
In this code, the 3D Dual-SLIP model is simulated for a set of optimal initial conditions found using the nonlinear optimizer by Liu et al. (2015), but for compliant terrain. The system is initiated using an optimal set of initial conditions found for the ideal rigid terrain, using the nonlinear optimizer by Liu et al. (2015). Initially, the ground stiffness was set to the rigid value of 50 MN/m. Then, after 10 steps the system experienced a one-step unilateral lower stiffness perturbation, after which the ground stiffness was set back to rigid. Specifically, at the 10th Touchdown event the ground stiffness under the leg about to land was lowered to a specific value and was kept constant throughout the whole step, while the ground stiffness of the leg in support remained to rigid. Next, after the completion of that step, the ground stiffness was set to rigid for both legs. In order to handle that perturbation, a biomechanics-inspired controller is proposed. Specifically, at the 10th Touchdown, the stiffness of the leg in support is preadjusted to a higher value, while the stiffness of the leg experiencing the perturbation is also increased throughout the step. At the 11th Midstance step, the stiffness of the leg about to land on rigid terrain is again set to default.

# Instructions
The main code is located in "three_D_Dual_SLIP_compliant_sim_MED_2022.m". The rest of the codes are custom made functions that are called during the execution of the code.

For more information please refer to the following paper:

**[Karakasis, Chrysostomos, Ioannis Poulakakis, and Panagiotis Artemiadis. "Robust Dynamic Walking for a 3D Dual-SLIP Model under One-Step Unilateral Stiffness Perturbations: Towards Bipedal Locomotion over Compliant Terrain." arXiv preprint arXiv:2203.07471 (2022).](https://arxiv.org/abs/2203.07471)**

![perturbation_screenshots_90kNm_final](https://user-images.githubusercontent.com/95447396/159397608-f42770cb-daae-45e7-968f-e54667598d74.png)
