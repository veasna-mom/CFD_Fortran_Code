program FTCS_2D
    implicit none
    integer :: Nx, Ny, i, j, nsteps, max_iter
    real :: Lx, Ly, T0, alpha, dx, dy, dt, tol, k, cp, rho, Fox, Foy
    real :: T_base, T_end, time
    real, allocatable :: T(:,:), T_new(:,:), x(:), y(:)
    integer :: output_interval

    ! ----------------- Parameters -----------------
    Lx = 1.0
    Ly = 1.0
    T0 = 25.0
    dx = 0.05
    dy = 0.05
    dt = 0.001
    Nx = int(Lx/dx) + 1
    Ny = int(Ly/dy) + 1
    k = 400.0
    cp = 385.0
    rho = 8000.0
    alpha = k/(rho*cp)
    Fox = alpha*dt/(dx*dx)
    Foy = alpha*dt/(dy*dy)
    T_base = 400.0
    T_end  = 400
    tol = 1.0e-4
    max_iter = 50000
    output_interval = 50

    ! Stability check
    if (Fox + Foy > 0.5) then
        print *, "Warning: Stability condition violated. Reduce dt."
        stop
    end if

    ! ----------------- Allocate arrays -----------------
    allocate(T(Nx,Ny), T_new(Nx,Ny), x(Nx), y(Ny))

    ! Initialize temperature and positions
    T = T0
    T_new = T0
    do i = 1, Nx
        x(i) = (i-1)*dx
    end do
    do j = 1, Ny
        y(j) = (j-1)*dy
    end do

    ! Apply initial boundary conditions (Dirichlet)
    T(1,:)   = T_base
    T(Nx,:)  = T_end
    T(:,1)   = T_base
    T(:,Ny)  = T_end
    T_new = T

    ! ----------------- Open Tecplot file -----------------
    open(unit=10, file='temperature2D.dat', status='replace')
    write(10,*) 'TITLE = "2D Heat Conduction - FTCS"'
    write(10,*) 'VARIABLES = "x" "y" "T"'

    ! ----------------- Write initial condition -----------------
    time = 0.0
    write(10,'(A,F6.3,A,I0,A,I0,A,F6.3)') &
         'ZONE T="time=', time, '", I=', Nx, ', J=', Ny, &
         ', DATAPACKING=POINT, SOLUTIONTIME=', time
    do j = 1, Ny
        do i = 1, Nx
            write(10,'(F8.4,1X,F8.4,1X,F8.4)') x(i), y(j), T(i,j)
        end do
    end do

    ! ----------------- Time Loop -----------------
    do nsteps = 1, max_iter
        time = nsteps * dt

        ! Update interior nodes
        do j = 2, Ny-1
            do i = 2, Nx-1
                T_new(i,j) = T(i,j) + Fox*(T(i+1,j)-2.0*T(i,j)+T(i-1,j)) &
                                       + Foy*(T(i,j+1)-2.0*T(i,j)+T(i,j-1))
            end do
        end do

        ! Apply boundary conditions
        T_new(1,:)   = T_base
        T_new(Nx,:)  = T_end
        T_new(:,1)   = T_base
        T_new(:,Ny)  = T_end

        ! Write to Tecplot file every output_interval steps
        if (mod(nsteps, output_interval) == 0 .or. nsteps == max_iter) then
            write(10,'(A,F6.3,A,I0,A,I0,A,F6.3)') &
                 'ZONE T="time=', time, '", I=', Nx, ', J=', Ny, &
                 ', DATAPACKING=POINT, SOLUTIONTIME=', time
            do j = 1, Ny
                do i = 1, Nx
                    write(10,'(F8.4,1X,F8.4,1X,F8.4)') x(i), y(j), T_new(i,j)
                end do
            end do
        end if

        ! Check convergence
        if (maxval(abs(T_new - T)) < tol) then
            print *, "Solution converged after", nsteps, "time steps"
            exit
        end if

        ! Update for next step
        T = T_new
    end do

    ! ----------------- Close file -----------------
    close(10)
    deallocate(T, T_new, x, y)

    print *, "Simulation complete. Tecplot file 'temperature2D.dat' created."

end program FTCS_2D
