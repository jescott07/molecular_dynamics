# n-body Gravitational System

As explained in the first example of this repository (see the [README](../ex1/README.md)), we will create a [module](src/ex2.f90) in Fortran and run it in Python. First, lets define the variables:

```fortran
integer, parameter :: n=1000000

integer, parameter :: L=10

integer, parameter :: p=10

integer :: seed

real(8), dimension(3,n,p) :: x,v

real(8),dimension(n) :: t,ec,ep,et

real(8) :: dt,ti

COMMON /rng/seed
```

Almost all the variables are as like the first example, but here we declare L which is the width of the box and the number of particles p. The common parameter seed is used to generate the random number as we will see.

Now, lets define a subroutine $\texttt{main()}$ that execute the Verlet algorithm and compute the energy.

```fortran
subroutine main()

    integer :: i,j,k

    real(8) :: dt2,r

    real(8),dimension(3) :: dx

    !Defining the initial conditions for each particle
    do i=1,p 
        x(1,1,i) = (L*(drandom()-0.5d0))
        x(2,1,i) = (L*(drandom()-0.5d0))
        x(3,1,i) = (L*(drandom()-0.5d0))
        v(1,1,i) = 0.0d0
        v(2,1,i) = 0.0d0
        v(3,1,i) = 0.0d0
    enddo

    ! For the second position:

    do i = 1,p
        x(:,2,i) = x(:,1,i) + v(:,1,i)*dt + 0.5d0*fs(1,i)*dt2
    enddo

    t(1) = ti
    t(2) = ti + dt

    ! Now let's run Verlet's algorithm for each particle at each point to compute x and v.
    dt2 = dt*dt

    do i=3,n
        do j=1,p
            x(:,i,j) = 2*x(:,i-1,j) - x(:,i-2,j) + fs(i-1,j)*dt2
            v(:,i-1,j) = (x(:,i,j) - x(:,i-2,j)) / (2.0d0*dt)
            t(i) = t(i-1) + dt
        enddo
    enddo

    ! The velocity in i=n is:

    do j=1,p
        v(:,n,j) = (x(:,n,j) - x(:,n-1,j)) / h
    enddo

    ! Computing the kinetic energy

    do i=1,n
        ec(i) = 0.0d0
        do j=1,p
            ec(i) = ec(i) + ecf(v(:,i,j))
        enddo
    enddo

    ! Computing the potential energy
    do i=1,n
        ep(i) = 0.0d0
        do j=1,p
            do k=1,p
                if (k .lt. j) then
                    dx(:) = x(:,i,j) - x(:,i,k)
                    r = norm2(dx)
                    ep(i) = ep(i) + epf(r)
                endif
            enddo
        enddo
    enddo

    ! The total energy in each point is:

    do i=1,n
        et(i) = ec(i) + ep(i)
    enddo

end subroutine main
```

What is most different about this example is the functions $\texttt{fs()}$, $\texttt{ecf}$ and $\texttt{epf}$ which compute the acceleration, the kinetic energy and the potential energy respectively, as show below:

```fortran
function fs(n_index,p_index)

    integer, intent(in) :: n_index,p_index

    integer :: i

    real(8),dimension(3) :: fp,dx,fs 

    real(8) :: r

    fs(:) = 0.0d0 

    do i=1,p
        if (i .ne. p_index) then
            dx(:) = x(:,n_index,i) - x(:,n_index,p_index)
            r = norm2(dx)
            fp(:) = (dx(:) / r**3) * (1.0d0 - (EXP(-r*r)*(1+2.0d0*r*r))) 
            fs(:) = fs(:) + fp(:)           
        endif
    enddo
    return
end function fs

function ecf(xf)
    real(8) :: ecf
    real(8),dimension(3) :: xf
    ecf = 0.5d0*norm2(xf)**2 
    return
end function ecf

function epf(xf)
    real(8) :: epf, xf  
    epf = -(1.0d0 - exp(-xf**2))/xf
    return
end function epf
```

Finely the function that computes the random number was taken from Konstantino's book: Computational Physics (2014).

```fortran
FUNCTION drandom()
    REAL(8)            :: drandom, dr
    INTEGER, PARAMETER :: a = 16807
    INTEGER, PARAMETER :: m = 2147483647
    INTEGER, PARAMETER :: q = 127773
    INTEGER, PARAMETER :: r = 2836
    REAL(8), PARAMETER :: f = 1.0d0/m
    INTEGER            :: p, seed
    COMMON /rng/seed

    101 CONTINUE
    p = seed/q
    seed = a*(seed-q*p)-r*p
    IF(seed .lt. 0) seed = seed + m
    dr = f*seed
    IF(dr .le. 0.0d0 .or. dr .ge. 1.0d0) GOTO 101
    drandom = dr
    RETURN
END FUNCTION drandom
```

To analyse the results, lets create two files in python, [animation.py](animation.py) and [energy.py](energy.py) which will analyse the space solution and the energy of the system respectively. For the animation part, lets import the module crated in Fortran and some useful libraries.

```python3
from grav_2 import grav as gv

from matplotlib import pyplot as plt

from celluloid import Camera

import numpy as np
```

Now lets define our parameters and run the $\texttt{main()}$ subroutine:

```python3 
gv.seed = 1369420
gv.ti = 0
gv.dt = 0.001

n = int(gv.n)

gv.main()

t = gv.t[:]
```

With $\texttt{Camera}$ from $\texttt{celluloid}$ and $\texttt{pyplot}$ from $\texttt{matplotlib}$ lets take a snap of the solution in different times and save it as a mp4 video:

```python3
fig = plt.figure()
ax = fig.gca(projection='3d')
camera = Camera(fig)
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')

X = np.zeros(int(gv.p))
Y= np.zeros(int(gv.p))
Z = np.zeros(int(gv.p))

k = 0
while k <= n:
    for i in range(len(X)): # Gravar a posição em j daz p partículas.
        x = gv.x[0,:,i]
        y = gv.x[1,:,i]
        z = gv.x[2,:,i]
        X[i] = x[k-1]
        Y[i] = y[k-1]
        Z[i] = z[k-1]
    k += 10000
    ax.scatter(X,Y,Z,marker='o',c='black')
    camera.snap()
animation = camera.animate()
animation.save('ex2_animation_l10_p10.mp4')
```
