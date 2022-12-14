# Molecular Dynamics

The molecular dynamics are a simulation technique largely used in physics, astrophysics, chemistry and materials science. It consists in integrate the equation of motion of an N classical particle system using the Verlet algorithm. This repository will be organized as follow: First we will present what is the Verlet algorithm, then we will apply it for two problems: A single particle under a central force and N-body gravitational problem.

## The Verlet Algorithm

In another repository we discussed more about differential equations an presented some techniques to solve them, as you can see in this [README](https://github.com/jescott07/solving-differential-equations/blob/main/README.md). If we were in a classical world, i.e. we are studying an not so small system that we would need to use quantum mechanics and also not so massive and fast enough that we needed to use general relativity. Then the equation of motion is given by the second Newton's Law, which is given by the equation:

<p align="center">
$F = m a$.
</p>

Where $F$ is the force being applied under the particle with mass $m$ and $a=\frac{dv}{dt} = \frac{d^2x}{dt^2}$ is the acceleration of that particle. Thus, the equation of motion is a second order differential equation:

<p align="center">
$F = m \ddot x(t)$.
</p>

Where $\ddot x(t)$ represents the second order derivative of x in relation to time (t).

The Verlet algorithm consists of first discretize the time $t$, such as:


<p align="center">
$t_i = t_0 + i\Delta t \qquad i = 1,...,N$.
</p>

Where $t_0$ is the initial time, $\Delta t$ the step size in time and N is the number of iterations. Thus  the numerical solution for each $t_i$ will be an approximation of the real solution $x_i \approx x(t_i)$.

So, our solution space will be discretized, as well as the acceleration, then:

<p align="center">
$a_i = \frac{\Delta^2 x}{\Delta t^2} = \frac{\frac{\Delta x_{i+1}}{\Delta t} - \frac{\Delta x_i}{\Delta t}}{\Delta t}$.
</p>

Where:

<p align="center">
$\frac{\Delta x_{i+1}}{\Delta t} = \frac{x_{i+1} - x_i}{\Delta t}$,
</p>

and

<p align="center">
$\frac{\Delta x_{i}}{\Delta t} = \frac{x_{i} - x_{i-1}}{\Delta t}$.
</p>


Then:

<p align="center">
$a_i = \frac{\frac{x_{i+1} - x_i}{\Delta t} - \frac{x_{i} - x_{i-1}}{\Delta t}}{\Delta t} = \frac{x_{i+1} - 2x_i + x_{i-1}}{\Delta t^2}$.
</p>

So, knowing $x_0, x_1$ and $a_i \forall i$ we can advance the solution as:

<p align="center">
$x_{i+1} = a_i \Delta t^2 + 2 x_i - x_{i-1}$
</p>

Remembering that:

<p align="center">
$a_i = \frac{F_i}{m}$
</p>

If the force is proportional to the velocity, we also need to compute its solution, remembering that:

<p align="center">
$v = \frac{\Delta x}{\Delta t}$
</p>

Then the solution $v_i \approx v(t_i)$ can be computed as:

<p align="center">
$v_i = \frac{x_{i+1} - x_{i-1}}{2\Delta t}$
</p>

So, if we know the mass and the force acting on the particle, then we will also know the acceleration of the particle and then, from the initial conditions, we can determine the position and the velocity of the particle in every discretized time. With this, we can compute the mechanical energy of the system, which is given by the equation:

<p align='center'>
$E_{mec} = K + U$
<p/>

Where, K is the kinetic energy and U is the potential energy. Which are given by the equations:

<p align='center'>
$K = \frac{1}{2} m v^2$,
<p/>

and

<p align='center'>
$U = -\int_{x_i}^{x_f} F_x dx$,
<p/>

which is the potential for a conservative force, in the direction of $x$, between $x_i$ and $x_f$.

We can implement this method through the following pseudo-algorithm:

**Algorithm** *Verlet*

**Input** $f(x, \dot x, t)$: Second order differential equation to be integrated (function), $x(0)$, $\dot x(0) = v(0), t_0$: initial conditions (float), $dt$: Integration step (float), $N$: Number of integration steps (positive integer).

**Output** $x(N)$, $t(N)$: Vectors with size N containing the solutions of the system.

**Compute x(i+1), v(i+1) and t(i+1) for each x(i) and t(i) from $f(x, \dot x, t)$ and $dt$ as discussed above.**

1. Define i = 1.
2. Define $x(0) = x_i$, $v(0) = v_i$ and $t(0) = t_i$ (the initial condition inputs)
3. Define $x(1) = x(0) + v(0) dt + \frac{1}{2}f(x(0),v(0),t_0) dt^2$
4. Define dt2 = $dt^2$
4. **do**

      x(i+1) = 2x(i) - x(i-1) + f(x(i), v(i), t(i))*dt2
      
      v(i+1) = (x(i+1) - x(i-1)) / 2dt
      
      t(i+1) = t(i) + h
      
      i = i + 1
      
    **while** i $\neq$ N + 1
    
 4. Return x(N), v(N) and t(N).
 
## A particle under a central force

A force directed to the center, which decays with the distance ( $r$ ) to the center, i.e. $F \propto \frac{1}{r^\lambda}$ is:

<p align='center'>
  $\vec F = - k' \frac{\hat r}{r^\lambda} = - k' \frac{\vec r}{r^{\lambda + 1}}$
</p>

Where, $\vec r$ is a three dimensional vector, $k'$ is the proportionality factor and $\lambda$ is the decay factor.

Thus, from the Newton's second law the acceleration will be:

<p align='center'>
$\vec a = \vec{\ddot{r}} = - \frac{k'}{m}\frac{\vec r}{r^{\lambda + 1}} \equiv -k \frac{\vec r}{r^{\lambda + 1}}$
<p/>

From this force, the energy is:

<p align='center'>
$K_i = \frac{1}{2}v_i^2$
</p>

and for $\lambda > 1$:

<p align='center'>
$U_i = \frac{-k}{\lambda - 1} r_i^{\lambda - 1}$
</p>

With that $\vec{a}$ and using the Verlet algorithm presented above, we computed the position $x_i$, the velocity $v_i$ and the energy of the system for every discretized time $t_i$, setting: $\Delta t = 0,01$, $k = 1000$, $\vec{r_0} = -100 \hat i + 100 \hat j - 100 \hat k$, $\vec{v_0} = 1\hat i + 1 \hat j + 1 \hat k$ and varying $\lambda$.

First, lets investigate the Newtonian gravitation where $\lambda = 2$. The space solution is shown in the figure below:

![](ex1/space_solution_lambda2.png)

And the solution in the $x \times y$ plane is:

![](ex1/x_times_y_plane_lambda2.png)

From these results we see that for the Newtonian gravitation, we have a closed elliptical orbit as we spected. And for the energy we have:

![](ex1/energy_lambda2.png)

Thus, the kinect and the potential energy oscillates, however, the mechanical energy is conserved since there is no dissipative force in the system.

For $\lambda = 1.5$ we have the following space solution:

![](ex1/space_solution_lambda15.png)

and for $\lambda = 2.5$:

![](ex1/space_solution_lambda25.png)

Notice that for $\lambda = 2.5$ the particle acts like a free particle, i.e. it is ejected from the orbit. This because as $F \propto \lambda^{-1}$ the force under the particle is weaker if we compare to others values of $\lambda$. So, to see the orbit, we need a greater value of $k$, thus lets define $k=10000$ and run the program again. So, we have the solution:

![](ex1/space_solution_lambda25_2.png)

To see in details how the program was build see the [README](ex1/README.md).

## n-body Gravitational System

Now we are gone study a very interesting problem, mainly in astrophysics. A n-body gravitational problem is a system, like a star cluster, with N point masses under gravitational forces. From the section above, we saw that for Newtonian gravitation, the potential is $U(r) \propto r^{-1}$, where $r$ is the distance between the masses. Notice that for $r \rightarrow 0$, $U(r) \rightarrow \infty$ which will lead to numerical problems. Then we need to use a modified potential.

In a physical system, the limit $r \rightarrow 0$ could lead to big changes in the whole system, such as a colision of stars. In this example, we will not consider this type of implications, instead we wanna find a potential where $U(r) \rightarrow 0$ for $r \rightarrow 0$ and for $r \gg 1$, $U(r) \propto r^{-1}$, i.e. is equal to the Newtonian gravitational potential. So, lets consider the following modified potential:

<p align='center'>
$U(r) = - m_1 m_2 \frac{1 - e^{-r^2}}{r}$
</p>

That fulfills all the prerequisites we wanted for a modified potential. But, to use the Verlet method we need the gravitational force, not the potential, so lets determine it. Remembering that:

<p align='center'>
$F(r) = - \nabla U(r)$
</p>

<p align='center'>
$\Rightarrow F(r) = - m_1 m_2 \nabla (\frac{1 - e^{-r^2}}{r}) = \frac{m_1 m_2}{r^2}[1 - e^{-r^2}(1+2r^2)]\hat r$
</p>

For a pair of particles $i$ and $j$:

<p align='center'>
$\vec r_{ij} = (\vec r_j - \vec r_i)$
</p>

Which is the difference between the two particles. Their module is them:

<p align='center'>
$r_{ij} = |\vec r_{ij}| = |(\vec r_j - \vec r_i)|$
</p>

Then the force that particle $i$ does on particle $j$ is:

<p align='center'>
$\vec F_{ij} = \frac{\vec r_j - \vec r_i}{r_{ij}^3} m_i m_j [1 - e^{-r^2}(1 + 2 r_{ij}^2)]$
</p>

Remembering that:

<p align='center'>
$\hat r_{ij} = \frac{(\vec r_j - \vec r_i)}{r_{ij}}$
</p>

And the total energy is:

<p align='center'>
$E = \frac{1}{2} \Sigma_{i=1}^{N} m v_i^2 + \Sigma_{i < j} U(r_{ij})$
</p>

To simplify the problem we will set the same mass for all particles and define $m = 1$.

For the initial conditions, we will set a random position inside a box with volume $L^3$ where $L \gg 1$. Thus, the position of the first particle will be:

<p align='center'>
$x_1 = L(R - 0.5)$
</p>

Where $R$ is a random number between 0 and 1. We will do the same for the N particles. The initial velocities will be set as null for all particles, so the center of mass velocity will also be zero as well as the total angular momentum.

To analyse how the system evolves, we defined $L = 10$, the number of particles $p = 10$ and we did an animation with the position of the particles in different times:

<p align='center'>

![](ex2/ex2_animation_l10_p10.gif)

</p>

To analyse how the energy of the system evolves, we defined a bigger system with more particles, so we defined $L = 100$, $p = 100$ and plotted the potential energy, the kinect energy and the total energy in each time space as show in the figure below.

![](ex2/ex2_energy_l100_p100.png)

From this image we see that at the beginning the kinetic and potential energy of the system varies a lot and then they stabilize around a value after a certain time. This occurs because at first the system particles are randomly positioned inside the box of volume $L^3$ and from that initial condition the system will evolve to more stable orbits.

To see in details how the program was built, see the [README](ex2/README.md).








