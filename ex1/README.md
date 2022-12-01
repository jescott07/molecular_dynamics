# A Particle Under a Central Force

To execute this example we build a module in Fortran and run it in Python using f2py. To know more about this see the f2py section on this [README](https://github.com/jescott07/solving-differential-equations/blob/main/conv/README.md).

As [discussed before](../README.md) we wanna compute two 3D vectors, the position and the velocity, besides the energy of the system to every discretized time. With the initial parameters t0 (initial time), dt (step size in time), k and lb ( $\lambda$ ) which are parameters of the equation of motion. So, first, let's define these variables:

```fortran
integer, parameter :: n=100000

real(8),dimension(3,n) :: x,v

real(8),dimension(n)  :: t,ec,ep,et

real(8),dimension(3) :: xi,vi

real(8) :: t0,dt,k,lb
```

Where, n is the number of steps and xi, and vi are an 3D vector with the initial position and velocity respectively.

Then we need to build a subroutine that will run the Verlet algorithm:

```fortran
subroutine vl()

    real(8) :: dt2

    integer :: i

    dt2 = dt*dt


    x(:,1) = xi(:)
    v(:,1) = vi(:)
    x(:,2) = xi(:) + vi(:)*dt + 0.5d0*f(xi(:))*dt2
    t(1) = ti
    t(2) = ti + dt

    do i=3,n
        x(:,i) = 2*x(:,i-1) - x(:,i-2) + f(x(:,i-1))*dt2
        v(:,i-1) = (x(:,i) - x(:,i-2)) / (2.0d0*dt)
        t(i) = t(i-1) + dt
    enddo

    v(:,n) = (x(:,n) - x(:,n-1))/dt

    do i=1,n
        ec(i) = ecf(v(:,i))
        ep(i) = epf(x(:,i))
        et(i)  = ec(i)+ep(i)         
    end do

end subroutine vl
```

where ecf(v) and epf(x) are the functions that compute the kinect and the potential energy respectevely, and f(x) is the function to be solved by the Verlet algorithm, i.e. $f(x) = \frac{\vec F}{m}$. This functions are built as follows:

```fortran
function ecf(xf)
    real(8) :: ecf
    real(8),dimension(3) :: xf
    ecf = 0.5d0*norm2(xf)**2 
    return
end function ecf

function epf(xf)
    real(8) :: epf
    real(8),dimension(3) :: xf    
    epf = -k / ((lb-1.0d0)*(norm2(xf)**(lb-1.0d0)))
    return
end function epf

function f(xf)

    real(8),dimension(3) :: f,xf

    f(:) = -k*xf(:) / (norm2(xf)**(lb+1.0d0))

    return

end function f
```

Now we need to build a program in Python run this model, named $\texttt{verlet}$ and plot the results.

To do this, first we need to wrap the Fortran program to build an extension module that will be imported on Python. To do this, let's write in the terminal:

```
$ f2py -c ex1.f90 -m ex1
```

Where, $\texttt{ex1.f90}$ is the name of the program and $\texttt{ex1}$ is the name of the extension module that we wanna create, which contain the $\texttt{verlet}$ module that we build (do not confuse extension module with module).

Then, from the extension module $\texttt{ex1}$ lets import the $\texttt{verlet}$ module, as well as other useful libraries.

```python3
from ex1 import verlet

import matplotlib.pyplot as plt

import numpy as np

from mpl_toolkits.mplot3d import Axes3D
```

Now, lets define the initial conditions and the parameters:

```python3
verlet.dt = 0.01
verlet.t0 = 0.
verlet.k = 1000
verlet.lb = 1.5
verlet.xi = [-100,100,-100]
verlet.vi = [1,1,1]
```

run the $\texttt{vl}$ subroutine and plot the results:

```python3
verlet.vl()

x = verlet.x[0,:]
y = verlet.x[1,:]
z = verlet.x[2,:]
t = verlet.t[:]
ec = verlet.ec[:]
ep = verlet.ep[:]
et = verlet.et[:]

fig = plt.figure()
ax = fig.gca(projection='3d')
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')
ax.plot(x, y, z)
ax.set_title('Space Solution for ' r'$\lambda = 1.5$')
ax.view_init(elev=10., azim=-25)
plt.savefig('space_solution_lambda15.png',dpi=300)
plt.close()

plt.figure()
plt.title('Solution in the 'r'$x \times y$' ' plane for ' r'$\lambda = 1.5$')
plt.xlabel('X')
plt.ylabel('Y')
plt.plot(x,y)
plt.savefig('x_times_y_plane_lambda15.png',dpi=300)
plt.close()

plt.figure()
plt.title('Energy ' r'$\times$' ' Time for ' r'$\lambda = 1.5$' )
plt.xlabel('Time')
plt.ylabel('Energy')
plt.ylim(min(ep) - 5, max(ec) + 15)
plt.plot(t,ec,label='Kinetic Energy')
plt.plot(t,ep,label='Potential Energy')
plt.plot(t,et,label='Total Energy')
plt.legend(fontsize='small')
plt.savefig('energy_lambda15.png',dpi=300)
plt.close()
```

We can do the same for others initial conditions and parameters.
