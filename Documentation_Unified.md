# Unified Benders cut
Fischetti (2010) proposed a unified Benders cut, in which a single expression can act as both an optimality and a feasibility cut. To construct this unified cut, we examine the following problem.

**1. Primal Subproblem**

$$
\begin{align*}
\min \quad & 0 \\
\text{s.t.} \quad & t^* \geq c^{\top}y \\
& Ax+By \geq b \\
& Cx+Dy \leq d \\
& Ex+Fy = f \\
& x = x^* \\
& y \geq 0
\end{align*}
$$

Infeasibility of the above problem implies that we can generate the unified cut which can effectively eliminate the given separation point $(x^*, t^*)$. 

**2. Primal Subproblem (canonical form) with slack variables**

$$
\begin{align*}
\min \quad & z \\
\text{s.t.} \quad & w_0z - c^{\top}y \geq -t^* & \quad (\pi_0) \\
& w_1z + Ax+By \geq b & \quad (\pi_1) \\
& w_2z -Cx-Dy \geq -d & \quad (\pi_2) \\
& w_3z + Ex+Fy \geq f & \quad (\pi_3)\\
& w_4z -Ex-Fy \geq -f & \quad (\pi_4) \\
& w_5z + x \geq x^* & \quad (\pi_5) \\
& w_6z -x \geq -x^* & \quad (\pi_6) \\
& y \geq 0
\end{align*}
$$

The dual of this subproblem is rooted in the alternative polyhedron, where the indices of an extreme point with strictly positive entries correspond to the primal constraints that contribute to infeasibility. The slack variable $w_i$ is used to normalize the objective value of the dual subproblem, since the dual region would otherwise be conic. Assigning a large $w_i$ implies a stronger penalty if the corresponding constraint contributes to infeasibility. Consequently, when the primal problem is infeasible, we can construct a cut by retrieving the dual values.

**3. Dual Subproblem**

$$
\begin{align*}
\max \quad & -t^*\pi_0 + b^{\top}\pi_1 - d^{\top}\pi_2 + f^{\top}(\pi_3 - \pi_4) + (x^*)^{\top}(\pi_5 - \pi_6) \\
\text{s.t.} \quad & -c^{\top}\pi_0 + B^{\top}\pi_1 - D^{\top}\pi_2 + F^{\top}(\pi_3- \pi_4) \leq 0 \\
& A^{\top}\pi_1 - C^{\top}\pi_2 + E^{\top}(\pi_3- \pi_4) + \pi_5 - \pi_6 = 0\\
& \sum_{i=1}^{6}w_i\bold{1}^{\top}\pi_i + w_0\pi_0 = 1 \\
& \pi_i \geq 0, \quad \forall i
\end{align*}
$$

By leveraging strong duality, we can generate following unified cut.

$$
0 \geq z^* + \pi_0^*t^* - (\pi_5^* - \pi_6^*)^{\top}x^* - \pi_0^*t + (\pi_5^* - \pi_6^*)^{\top}x
$$

