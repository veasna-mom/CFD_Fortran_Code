program FTCS
    implicit none
    integer :: N, i, nsteps, max_iter
    real :: L, T0, alpha, dx, dt, tol, k, cp, rho, Fo
    real :: T_base, T_end, time
    real, allocatable :: T(:), T_new(:), x(:)
    integer :: output_interval

    ! ----------------- Parameters -----------------
    L = 1.0
    T0 = 25.0
    dx = 0.1
    dt = 1.0
    N = 11  ! Number of points (0.0 to 1.0)
    k = 400.0
    cp = 385.0
    rho = 8000.0
    alpha = k/(rho*cp)
    Fo = alpha*dt/(dx*dx)
    T_base = 400.0
    T_end = 400.0
    tol = 1.0e-4
    max_iter = 900
    output_interval = 100

    ! ----------------- Allocate arrays -----------------
    allocate(T(N), T_new(N), x(N))

    ! Initialize temperature and positions
    T = T0
    T_new = T0
    do i = 1, N
        x(i) = (i-1)*dx
    end do

    ! Apply initial boundary conditions
    T(1) = T_base
    T(N) = T_end
    T_new(1) = T_base
    T_new(N) = T_end

    ! ----------------- Open Tecplot file -----------------
    open(unit=10, file='temperature.dat', status='replace')
    write(10,*) 'TITLE = "1D Heat Conduction - FTCS"'
    write(10,*) 'VARIABLES = "x" "T"'

    ! ----------------- Write initial condition -----------------
    time = 0.0
    write(10,'(A,F6.2,A,I0,A,F6.2)') 'ZONE T="time=', time, '", I=', N, ', DATAPACKING=POINT, SOLUTIONTIME=', time
    do i = 1, N
        write(10,'(F8.4,1X,F8.4)') x(i), T(i)
    end do

    ! ----------------- Time Loop -----------------
    do nsteps = 1, max_iter
        time = nsteps * dt

        ! Update interior nodes (explicit FTCS)
        do i = 2, N-1
            T_new(i) = T(i) + Fo*(T(i+1) - 2.0*T(i) + T(i-1))
        end do

        ! Apply boundary conditions
        T_new(1) = T_base
        T_new(N) = T_end

        ! Write to Tecplot file every output_interval steps
        if (mod(nsteps, output_interval) == 0 .or. nsteps == max_iter) then
            write(10,'(A,F6.2,A,I0,A,F6.2)') 'ZONE T="time=', time, '", I=', N, ', DATAPACKING=POINT, SOLUTIONTIME=', time
            do i = 1, N
                write(10,'(F8.4,1X,F8.4)') x(i), T_new(i)
            end do
        end if

        ! Check convergence BEFORE updating T
        if (maxval(abs(T_new - T)) < tol) then
            print *, "Solution converged after", nsteps, "time steps"
            if (mod(nsteps, output_interval) /= 0) then
                write(10,'(A,F6.2,A,I0,A,F6.2)') 'ZONE T="time=', time, '", I=', N, ', DATAPACKING=POINT, SOLUTIONTIME=', time
                do i = 1, N
                    write(10,'(F8.4,1X,F8.4)') x(i), T_new(i)
                end do
            end if
            exit
        end if

        ! Update temperature for next step
        T = T_new
    end do

    ! ----------------- Close file and deallocate arrays -----------------
    close(10)
    deallocate(T, T_new, x)

    print *, "Simulation complete. Tecplot file 'temperature.dat' created."

end program FTCS
