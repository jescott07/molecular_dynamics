# Molecular Dynamics

The molecular dynamics is a simulation technique largely used in physics, astrophysics, chemistry and materials science. It consists in integrate the equation of motion of N classical particles system using the Verlet algorithm. This repository will be organized as follow: First we will present what is the Verlet algorithm, then we will apply it for two problems: A particle under a central force and N-body gravitational problem.

## The Verlet Algorithm

In another repository we discussed more about differential equations an presented some techniques to solve them, as you can see in this [README](https://github.com/jescott07/solving-differential-equations/blob/main/README.md). If we where in a classical world, i.e. we are studing an not so small systens that we would need to use quantum mechanics and also not so massive and fast enoth that we needed to use general relativity. Thus the equation of motion is given by the second Newton's Law, wich is given by the equation:

<p align="center">
$F = m a$
</p>

Where $F$ is the force being apply under the particle with mass $m$ and $a=\frac{dv}{dt} = \frac{d^2x}{dt^2}$ is the aceleration of that particle. Thus, the equation of motion ia a second order differential equation:

<p align="center">
$F = m \ddot x(t)$
</p>

Where $\ddot x(t)$ represents the second order derivative of x in relation to time (t).

The Verlet algorithm consist at first discretize the time $t$, such:


<p align="center">
$t_i = t_0 + i\Delta t \qquad i = 1,...,N$
</p>

Where $t_0$ is the inicial time, $\Delta t$ the step size in time and N is the number of iterations. Thus  the numerical solution for each $t_i$ is a approximation of the real solution $x_i \approx x(t_i)$, then our solution space will be discretized as well as the acceleration, so:

<p align="center">
$a_i = \frac{\Delta^2 x}{\Delta t^2} = \frac{\frac{\Delta x_{i+1}}{\Delta t} - \frac{\Delta x_i}{\Delta t}}{\Delta t}$
</p>

Where:

<p align="center">
$\frac{\Delta x_{i+1}}{\Delta t} = \frac{x_{i+1} - x_i}{\Delta t}$
</p>

and

<p align="center">
$\frac{\Delta x_{i}}{\Delta t} = \frac{x_{i} - x_{i-1}}{\Delta t}$
</p>


Then:

<p align="center">
$a_i = \frac{\frac{x_{i+1} - x_i}{\Delta t} - \frac{x_{i} - x_{i-1}}{\Delta t}}{\Delta t} = \frac{x_{i+1} - 2x_i + x_{i-1}}{\Delta t^2}$
</p>

So, knowing $x_0, x_1$ and $a_i \forall i$ we can advance the solution as:

<p align="center">
$x_{i+1} = a_i \Delta t^2 + 2 x_i - x_{i-1}$
</p>

Remembering that:

<p align="center">
$a_i = \frac{F_i}{m}$
</p>

So, if we know the mass and the force acting on the particle then we will also know the acceleration of the particle and then we can determine the position of the particle in every discretized time if we know the incial condition.

We can implement this method through the following pseudo-algorithm:

**Algorithm** *Verlet*

**Input** $f(x, \dot x, t)$: Second order differential equation to be integrated (function), $x(0)$, $\dot x(0) = v(0), t_0$: initial conditions (float), $dt$: Integration step (float), $N$: Number of integration steps (positive integer).

**Output** $x(N)$, $t(N)$: Vectors with size N containing the solutions of the system.

**Compute x(i+1), t(i+1) for each x(i) and t(i) from $f(x, \dot x, t)$ and $dt$ as discussed above.**

1. Define i = 1.
2. Define $x(0) = x_i$, $t(0) = t_i$ (the initial condition inputs)
3. Define $x(1) = x(0) + v(0) dt + \frac{1}{2}f(x(0),v(0),t_0)$ dt^2
4. **do**

      x(i+1) = 2x(i) - x(i-1) + f(x(i), v(i), t(i))*dt**2
      
      v(i+1) = (x(i+1) - x(i-1)) / 2dt
      
      t(i+1) = t(i) + h
      
      i = i + 1
      
    **while** i $\neq$ N + 1
    
 4. Return x(N) and t(N).

## A particle under a central force

A force directed to the center wich decay with the distance (r) to the center, i.e. $\propto \frac{1}{r^\lambda}$ is:

<p align='center'>
  $\vec F = - k \frac{\hat r}{r^\lambda} = - k \frac{\vec r}{r^{\lambda + 1}}$
</p>

Where, $\vec r = x \hat i + y \hat j + z \hat k$, $k = m$ and $\lambda$ is the decay factor. For Newtonian Gravitation $\lambda = 2$.
