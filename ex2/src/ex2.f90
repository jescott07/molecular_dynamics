module grav

    integer, parameter :: n=1000000

    integer, parameter :: L=10

    integer, parameter :: p=10

    integer :: seed

    real(8), dimension(3,n,p) :: x,v

    real(8),dimension(n) :: t,ec,ep,et

    real(8) :: dt,ti

    COMMON /rng/seed

    contains

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

    ! Gerador de número aleatório tomado do livro computacional physics (Konstantinos, 2014).
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


end module
