# Pareto-optimal cut
Pareto-optimal cut is introduced by Magnanti and Wong (1981), where it is motivated by the degenaracy of the subproblem. As multiple dual solutions correspond to an optimal objective value, multiple Benders cuts which may have dominance relationship can exist. A cut corresponds to $\pi_1 \in \Pi$ dominates that of $\pi_2 \in \Pi$ if it satisfies the following inequality. 

$$
\begin{aligned} 
(b-x)^{\top}\pi_1 \geq (b-x)^{\top}\pi_2, \quad \forall x \in X
\end{aligned}
$$

It implies that any $\pi \in \Pi$ that generates Pareto-optimal cut cannot be $\pi_2$.

# Magnanti-Wong probem
A dual solution corresponding to a pareto-optimal cut can be obtained by solving the Magnanti-Wong problem. To this end, we should solve primal subproblem first to obtain optimal objective value $\xi^*$. The Magnanti-Wong problem finds a pareto-optimal cut that ensures to have $\xi^*$ at $x^*$, and to be as tight as possible at a core point $x_0$. 

**Dual Subproblem**

$$
\begin{align*}
\max \quad & (b - Ax^*)^{\top}\pi \\
\text{s.t.} \quad & B^{\top}\pi \leq c \\
& \pi \geq 0
\end{align*}
$$

**Magnanti-Wong Problem**

$$
\begin{align*}
\max \quad & (b - Ax_0)^{\top}\pi \\
\text{s.t.} \quad & B^{\top}\pi \leq c \\
& (b - Ax^*)^{\top}\pi = \xi^* \\
& \pi \geq 0
\end{align*}
$$

Since __BendersDecomposition.jl__ implements primal subproblem with the variable-fixing constraint, we need to derive the Magnanti-Wong primal problem that corresponds to our primal subproblem. The following describes the derivation process of the Magnanti-Wong primal problem, where the parentheses mean the corresponding dual variable. 

**1. Primal Subproblem with the variable-fixing constraint**

$$
\begin{align*}
\min \quad & d^{\top}y \\
\text{s.t.} \quad & Ax + By \geq b & \quad (\pi_0) \\
& Cx + Dy = k & \quad (\pi_1) \\
& x = x^*  & \quad (\pi_2)\\
& y \geq 0 
\end{align*}
$$

**2. Dual Subproblem with the variable-fixing constraint**

$$
\begin{align*}
\max \quad & b^{\top}\pi_0 + k^{\top}\pi_1 + (x^*)^{\top}\pi_2 \\
\text{s.t.} \quad & A^{\top}\pi_0 + C^{\top}\pi_1 + \pi_2 = 0 & \quad (x) \\
& B^{\top}\pi_0 + D^{\top}\pi_1 \leq d & \quad (y) \\
& \pi_0 \geq 0 
\end{align*}
$$

**3. Magnanti-Wong Problem**

$$
\begin{align*}
\max \quad & b^{\top}\pi_0 + k^{\top}\pi_1 + (x_0)^{\top}\pi_2 \\
\text{s.t.} \quad & A^{\top}\pi_0 + C^{\top}\pi_1 + \pi_2 = 0 & \quad (x) \\
& B^{\top}\pi_0 + D^{\top}\pi_1 \leq d & \quad (y) \\
& b^{\top}\pi_0 + k^{\top}\pi_1 + (x^*)^{\top}\pi_2 = \xi^* & \quad (z) \\
& \pi_0 \geq 0 
\end{align*}
$$

**4. Magnanti-Wong Primal Problem**

$$
\begin{align*}
\min \quad & d^{\top}y + \xi^*z \\
\text{s.t.} \quad & Ax + By + bz \geq b & \quad (\pi_0) \\
& Cx + Dy + kz = k & \quad (\pi_1) \\
& x + x^*z = x_0 & \quad (\pi_2) \\
& y \geq 0 
\end{align*}
$$

As discussed in previous section, the Benders cut from dual subproblem with the variable-fixing constraint is as follows.

$$
\begin{aligned} 
t \geq d^{\top}y^* - (\pi_2^*)^{\top}x^* + (\pi_2^*)^{\top}x
\end{aligned}
$$

Therefore, a valid pareto-optimal cut can be obtained by using optimal dual value $\pi_2^*$ retrieved from Magnanti-Wong Primal Problem.

# Reference: Magnanti-Wong problem (relaxed version)
We begin this section from Magnanti-Wong problem described in the previous section.

**3-1. Magnanti-Wong Problem (relaxed)**

$$
\begin{align*}
\max \quad & b^{\top}\pi_0 + k^{\top}\pi_1 + (x_0)^{\top}\pi_2 \\
\text{s.t.} \quad & A^{\top}\pi_0 + C^{\top}\pi_1 + \pi_2 = 0 & \quad (x) \\
& B^{\top}\pi_0 + D^{\top}\pi_1 \leq d & \quad (y) \\
& b^{\top}\pi_0 + k^{\top}\pi_1 + (x^*)^{\top}\pi_2 \geq \xi^* -\epsilon & \quad (z) \\
& \pi_0 \geq 0 
\end{align*}
$$

**4-1. Magnanti-Wong Primal Problem (relaxed & not canonical)**

$$
\begin{align*}
\min \quad & d^{\top}y + (\xi^* - \epsilon)z \\
\text{s.t.} \quad & Ax + By + bz \geq b & \quad (\pi_0) \\
& Cx + Dy + kz = k & \quad (\pi_1) \\
& x + x^*z = x_0 & \quad (\pi_2) \\
& y \geq 0, z \leq 0
\end{align*}
$$

The above can be expressed as following canonical form.

**4-2. Magnanti-Wong Primal Problem (relaxed & not canonical)**

$$
\begin{align*}
\min \quad & d^{\top}y - (\xi^* - \epsilon)z \\
\text{s.t.} \quad & Ax + By - bz \geq b & \quad (\pi_0) \\
& Cx + Dy - kz = k & \quad (\pi_1) \\
& x - x^*z = x_0 & \quad (\pi_2) \\
& y,z \geq 0
\end{align*}
$$