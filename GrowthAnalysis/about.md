## Parameters & Variables


### Length $l(t)$:

Must be a strictly increasing function. One example is the simplest linear growth scheme with $l(t) = l_0 + \dot{l} t$.


### Radius $r(s,t)$:

The radius is decided by one function $\tau$ as:

$$ r(s,t) = \tau(l(t) - s) $$

where $\tau$ is the **taper function** defined in $[0, L]$, which links the distance measured from the tip to the radius. $\tau(0)$ is the radius at the tip, the point of primary growth. Basically $\tau(x)$ will be a graph of the final radius distribution at time $T$ and length $L$ with flipped direction. A condition that must be met for $\tau$ is $\tau' \geq 0$, for all $x \in [0, L]$ (because the radius should never decrease). One example may be a linear scheme as $\tau(x) = r_{tip} + \alpha x$ with $alpha$ determining the taper slope.


### Local angle $\theta(s,t)$ and actual curvature $\kappa(s,t)$

$\theta(s,t)$ is the angle between the horizontal ($x$-axis positive direction) and the tangent (material positive direction) at point $s$ at time $t$.
Its spatial derivative is the actual curvature, $\partial_s \theta(s,t) = \kappa(s,t)$.


### Intrinsic curvature $\kappa_g(s,t)$

Intrinsic curvature $\kappa_g(s,t)$ is the curvature in the loadless state. The applied moment on the beam will be related to the difference between the two curvatures, $\kappa - \kappa_g$.

### Other parameters

- Density $\rho$
- Gravitational acceleration $g$
- Young's modulus $E$ ... maybe this should be variable in the future


## General scheme


### Quasi-static equilibrium

For every time $t$ the equilibrium holds. Obtained from the assumptions of Euler-Bernoulli beam, with the second moment of inertia as $I = \pi Er^4/4$

$$ \dfrac{\pi E}{4} \partial_s \left[ r^4 (\kappa - \kappa_g) \right] = \rho g \cos{\theta} \int_s^{l(t)} \pi r^2 \; d\tilde{s} $$
$$ \implies \partial_s \left[ r^4 (\kappa - \kappa_g) \right] = \dfrac{4 \rho g}{E} \cos{\theta} \int_s^{l(t)} r^2 \; d\tilde{s} $$

The boundary conditions are:

- $\theta(0,t) = \theta_0$ (clamped end to the main trunk)

- $\kappa(l(t),t) = \kappa_g(l(t),t)$ (no moment at tip)


### Evolution of intrinsic curvature

From the assumptions of the newly added layers, we have the evolution equation of the intrinsic curvature $\kappa_g$:

$$ \partial_t \kappa_g = 4 \dfrac{\partial_t r}{r} (\kappa - \kappa_g) + 48 \pi \dfrac{\partial_t r}{r^2} \epsilon(\theta,\kappa) $$

with the initial condition $\kappa_g(s,0) = 0$ (no intrinsic curvature in the beginning). $\epsilon$ is a given function which denotes the stress strain created through the creation of reaction wood: it gives some kind of control to the branch (such as gravitropism). To model the recovery of branch posture to some preferred angle $\theta_P$, one can let $\epsilon = \epsilon_{max} \sin{(\theta_P - \theta)}$.

Using these two differential equations (one spacial, one temporal), $\theta(s,t)$ can be obtained to know the shape of the branch.



## Non-dimensionalization

We modify the variables and equations to make it into a non-dimensional scheme. We use the non-dimensional arclength

$$ \eta = \dfrac{s}{l(t)}, \; \eta \in [0,1]$$

to represent the spacial position in the branch; the base is $\eta = 0$, the tip is $\eta = 1$ for any time. Thus, everything can be represented by the pair ($\eta$, $t$) instead of ($s$, $t$). The functions will be non-dimensionalized and change their representation as follows:

- $\Theta(\eta,t) = \theta(s,t)$

- $K(\eta,t) = l(t) \kappa(s,t) = \partial_\eta \Theta(\eta,t)$

- $K_g(\eta, t) = l(t) \kappa_g(s,t)$

- $R(\eta, t) = \dfrac{r(s,t)}{r(0,t)} = \dfrac{\tau(l(t)-s)}{\tau(l(t))}$


### Quasi-static equilibrium (non-dim)

$$ \partial_s \left[ r^4 (\kappa - \kappa_g) \right] = \dfrac{4 \rho g}{E} \cos{\theta} \int_s^{l(t)} r^2 \; d\tilde{s} $$

$$ \implies \left[ R^4 (K - K_g) \right]' = \dfrac{4 \rho g [l(t)]^3}{E [r(0,t)]^2} \cos{\Theta} \int_\eta^{1} R^2 \; d\tilde{\eta} $$

where $\partial_\eta (\cdot) = (\cdot)'$.

We can let $k(t) = [l(t)]^3/[r(0,t)]^2 = [l(t)]^3/[\tau(l(t))]^2$ and finally write

$$ \implies \left[ R^4 (K - K_g) \right]' = 4 \rho g k(t) \cos{\Theta} \int_\eta^{1} R^2 \; d\tilde{\eta} $$


### Evolution of intrinsic curvature (non-dim)

$$ \partial_t \kappa_g = 4 \dfrac{\partial_t r}{r} (\kappa - \kappa_g) + 48 \pi \dfrac{\partial_t r}{r^2} \epsilon(\theta,\kappa) $$

$$ \implies \dot{K_g} = \dfrac{\dot{l}}{l} \left[ K_g + \eta K_g' - \dfrac{R'}{R} \left(4(K-K_g) + \dfrac{48 \pi \epsilon}{R} \dfrac{l}{\tau(l)}\right) \right] $$

### Conditions (non-dim)

- $\Theta(0,t) = \theta_0$

- $K(1,t) = K_g(1,t)$

- $K(\eta,0) = 0$

### Numerical scheme

We solve the equilibrum equation with a ODE solver of scipy (scipy.bvp_solve) treating $t$ as a constant. Then we update $K_g$ with the evolution equation with small time step, and keep repeating these two processes.

After $\Theta(\eta,t)$ is obtained, one can draw the shape of the branch using the relations

$$ x(\eta,t) = l(t) \int_0^{\eta} \cos{\Theta} \, d\tilde{\eta} $$
$$ y(\eta,t) = l(t) \int_0^{\eta} \sin{\Theta} \, d\tilde{\eta} $$

along with the information $r(s,t)$.
