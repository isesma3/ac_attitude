# Attitude Dynamics and Kinematics: 3-2-1 Euler Angles with Quaternions

This MATLAB script simulates the **attitude dynamics and kinematics** of a rigid body using both **Euler angles (3-2-1 sequence)** and **quaternions**. It solves Euler's equations of motion for a rigid body under torque and compares orientation results using two different representations.

---

## Overview

The script models the rotational behavior of a rigid body using:

- **Euler’s rotational equations** to calculate angular velocity evolution.
- **321 Euler angles** and **quaternions** to represent and track the orientation.
- Numerical integration via MATLAB's `ode45` for both dynamics and kinematics.
- Comparison of quaternion-based orientation with Euler angle reconstruction.

This is useful for studying attitude representation, the benefits of quaternions over Euler angles, and understanding rotational dynamics under applied torque.

---

## Features

- Solves nonlinear dynamics for a rotating rigid body
- Implements torque input as a sinusoidal function
- Converts quaternion output to DCM and Euler angles
- Plots:
  - Angular velocity components
  - Quaternion components
  - Orientation via 3-2-1 Euler angles

---

## How to Run

1. Open MATLAB and navigate to the directory containing `project_B.m`.
2. Run:

```matlab
project_B
