program dykstra_solution_l2
    use iso_fortran_env
    implicit none

    !Scalar arguments
    integer(int32)    :: n, i, j, k, iter, max_iter, n_real, allocerr, ny, nz, nextnnyi, m_small
    integer(int64)    :: l, countny
    real(real64)      :: tol, element1, element2, stp, initial_time, final_time, s, avgny, res
    
    !Array arguments
    integer(int32)          :: p(3)
    integer(int64), pointer :: ind_y(:), ind_z(:)
    real(real64)            :: x_tmp(3), x_res(3), y_old(3), cntrs(3)
    real(real64), pointer   :: x(:), y(:), z(:), a(:), w(:)
    
    !Open input file
    character(len = 80)     :: inputfnm
    write(*, *) 'Enter the input file: '
    read(*,'(A)') inputfnm  
    inputfnm = ''

    open(10, file = trim(adjustl(inputfnm)), status = "old", action = "read")
    
    !Read the value of n_real from the file
    read(10, *) n_real

    !m = (n_real * (n_real - 1)) / 2 * (n_real - 2)
    n = n_real * (n_real - 1) / 2
    m_small = 357913941
    
    !Allocate matrices A and W
    allocate(a(n), w(n), x(n), y(m_small), ind_y(m_small), z(m_small), ind_z(m_small), stat = allocerr)

    !Check if allocation succeeded
    if ( allocerr .ne. 0 ) then
        write(*, *) 'Allocation error in main program.'
         stop
    end if

    k = 0
    do j = 1, n_real
        do i = 1, n_real
            read(10, *) element1, element2
            if (i .lt. j) then
                k = k + 1
                a(k) = element1
                w(k) = element2
            end if
        end do
    end do
    
    if (k .ne. n) then
        write(*, *) 'There is something wrong here (2).'
        stop
    end if

    close(10)
    
    call CPU_TIME(initial_time)

    !Initial value
    where (a(1 : n) .eq. 0.0d0)
        x(1:n) = 0.0d0
    elsewhere 
        x(1 : n) = a(1 : n) * w(1 : n) 
    end where

    cntrs = (/ -1.0d0, 1.0d0, 1.0d0 /)
    
    !Maximum number of iterations
    max_iter = 1000000
    
    !Set iter
    iter = 0
    
    !Tolerance criteria
    tol = 1.0d-8

    !Initial number of y
    ny = 0
    countny = 0
    
100 continue

    countny = countny + ny
    nz = 0
    l = 0
    iter = iter + 1
    stp = 0.0d0

    if (ny .eq. 0) then
       nextnnyi = - 1
    else
       nextnnyi = 1
    end if
    
    do i = 1, n_real - 1
        do j = i + 1, n_real
            do k = 1, n_real
                if (k .ne. i .and. k .ne. j) then
                     l = l + 1

                     if (nextnnyi .eq. - 1) then
                        y_old(1 : 3) = 0.0d0
                     else
                        if (l .lt. ind_y(nextnnyi)) then
                           y_old(1 : 3) = 0.0d0

                        else !if (l .eq. ind_y(nextnnyi)) then
                           y_old(1 : 3) = y(nextnnyi) * cntrs
                           nextnnyi = nextnnyi + 1

                           if (nextnnyi .gt. ny) then
                              nextnnyi = - 1
                           end if
                        end if
                     end if

                     if (k .lt. i .and. i .lt. j) then
                        !xij <= xki + xkj
                        call index_vec(i, j, p(1))
                        call index_vec(k, i, p(2))
                        call index_vec(k, j, p(3))
                
                     else if (i .lt. k .and. k .lt. j) then
                        !xij <= xik + xkj
                        call index_vec(i, j, p(1))
                        call index_vec(i, k, p(2))
                        call index_vec(k, j, p(3))
                        
                     else if (i .lt. j .and. j .lt. k) then
                        !xij <= xik + xjk
                        call index_vec(i, j, p(1))
                        call index_vec(i, k, p(2))
                        call index_vec(j, k, p(3))
                     end if
 
                    x_tmp(1 : 3) = x(p(1 : 3)) - y_old(1 : 3)

                    call project_triangular_inequality(x_tmp, res)
                    
                    if (res .ne. 0.0d0) then
                       nz = nz + 1

                       if (nz .gt. m_small) then
                          write(*, *) 'Too many non-null y. Increase m_small and re-run.'
                          stop
                       end if

                       ind_z(nz) = l
                       z(nz) = res
                    end if

                    x_res(1 : 3) = res * cntrs
                    
                    stp = stp + sum((y_old(1 : 3) - x_res(1 : 3)) ** 2)
                    
                    x(p(1 : 3)) = x_tmp(1 : 3) + x_res(1 : 3)
                end if
             end do
        end do
    end do

    write(*, *) "At", iter, "the stop criteria is", stp

    if (stp .le. tol) then
        x(1 : n) = x(1 : n) 
       
        s = norm2(w * (x - a))

        call CPU_TIME(final_time)

        avgny = countny / iter
        call write_results(inputfnm, n_real, iter, final_time - initial_time, avgny, s)
        return
    end if

    if (iter .gt. max_iter) then
        print *, 'Reached maximum iterations without convergence.'
        go to 500
    endif

    ny = nz
    y(1 : ny) = z(1 : nz)
    ind_y(1 : ny) = ind_z(1 : nz)

    go to 100
    
500 continue

    deallocate(a, w, y, x, stat = allocerr)

    if (allocerr .ne. 0) then
        write(*, *) 'Deallocation error in main program.'
         stop
    end if

contains

   
    subroutine project_triangular_inequality(v, s)
        use iso_fortran_env
        implicit none 

        !Scalar arguments
        real(real64) :: adotv

        !Array arguments
        real(real64), intent(in)  :: v(3) 
        real(real64), intent(out) :: s
       
        !Inner product
        adotv = v(1) - v(2) - v(3)

        !Check if satisfies the inequality with b = 0
        s = 0.0d0

        if (adotv .gt. 0.0d0) then
           s = adotv / 3 !Norm of a
        endif
    end subroutine project_triangular_inequality

    subroutine index_vec(i, j, k)
        use iso_fortran_env
        implicit none 

        !Scalar arguments
        integer(int32), intent(in)  :: i, j
        integer(int32), intent(out) :: k
     
        k = ((j - 1) * j) / 2 - j + i + 1
     
    end subroutine index_vec  
    
    subroutine write_results(file, n_real, n_iter, final_time, avgny, objfnc)
        use iso_fortran_env
        implicit none

        !Scalar arguments
        integer(int32)    :: unit_number, size1, n_real, n_iter
        real(real64)      :: final_time, objfnc, avgny

        !Character arguments
        character(len = 80) :: file
        character(len = 80) :: outputfnm

        !Prompt for the output file
        write(*,*) 'Enter the output file: '
        read(*,'(A)') outputfnm
        !outputfnm = ''

        unit_number = 10

        open(unit = unit_number, file = outputfnm, action = "write", position = "append", status = "unknown")

        ! Check if the file is empty and add headers
        inquire(file = outputfnm, size = size1)
        if (size1 == 0) then
            write(unit_number, '(A)') "File       | Size of A | Iterations |" // &
                                    "    CPU Time     | Average size of y | Objective Function"
        end if

        write(unit_number, '(A10,3X,I12,3X,I12,3X,F15.4,3X,F20.6,3X,F20.6)') &
            trim(file), n_real, n_iter, final_time, avgny, objfnc

        close(unit_number)
    end subroutine write_results
    
end program dykstra_solution_l2
