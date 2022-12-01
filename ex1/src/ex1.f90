module verlet 

    integer, parameter :: n=100000

    real(8),dimension(3,n) :: x,v

    real(8),dimension(n)  :: t,ec,ep,et

    real(8),dimension(3) :: xi,vi

    real(8) :: t0,dt,k,lb

    contains

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

end module
