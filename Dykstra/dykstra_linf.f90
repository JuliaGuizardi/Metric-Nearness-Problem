program dykstra_solution_l1
    use iso_fortran_env
    implicit none

    ! Scalar arguments
    integer(int32)    :: n, i, j, k, iter, max_iter, n_real, allocerr, m_small, ny, nz, nextnnyi
    integer(int64)    ::  l, countny
    real(real64)      :: tol,  gamma, adotx, theta, stp, element1, element2, initial_time, final_time, y_old, & 
                     eps, linear_of, quadratic_of, dual_of, max_inf, cmpl, optm, by

    ! Array arguments
    integer(int32)          :: p(3)
    integer(int64), pointer :: ind_y(:), ind_z(:)
    real(real64)            :: cntrs(3)
    real(real64), pointer   :: x(:), y(:), z(:), winv(:), weight(:), a(:), aypc(:)

    !Open input file
    character(len = 80)     :: inputfnm
    write(*, *) 'Enter the input file: '
    read(*,'(A)') inputfnm  
    !inputfnm = 'test3_8.txt'  

    open(10, file = trim(adjustl(inputfnm)), status = "old", action = "read")
    
    !Read the value of n_real from the file
    read(10, *) n_real

    n = (n_real * (n_real - 1)) / 2 + 1
    m_small = (n_real * (n_real - 1)) / 2 * (n_real - 2) + 2 * (n - 1)

    !Allocate 
    allocate(winv(n), x(n), aypc(n), a(n - 1), weight(n - 1), &
            y(m_small), ind_y(m_small), z(m_small), ind_z(m_small), stat = allocerr)

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
                weight(k) = element2
            end if
        end do
    end do

    if (k .ne. n - 1) then
        write(*, *) 'There is something wrong here (2).'
        stop
    end if

    close(10)

    call CPU_TIME(initial_time)

    ! Set gamma
    gamma = 2.0d0

    winv(1 : n - 1) =  gamma * 1.0d0 
    winv(n) = gamma * 1.0d0 

    ! Initial value of x
    x(1 : n - 1) = 0.0d0 
    x(n) = - winv(n) 

    !Maximum number of iterations
    max_iter = 100000
    
    !Set iter
    iter = 0
    
    !Tolerance criteria
    tol = 1.0d-8
    eps = 1.0d-6

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
                        y_old = 0.0d0
                     else
                        if (l .lt. ind_y(nextnnyi)) then
                           y_old = 0.0d0

                        else !if (l .eq. ind_y(nextnnyi)) then
                           y_old = y(nextnnyi)
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

                    cntrs = (/ 1.0d0, - 1.0d0, - 1.0d0 /)

                    ! Correction step
                    x(p(1 : 3)) = x(p(1 : 3)) + y_old * winv(p(1 : 3)) * cntrs(1 : 3)

                    ! Projection step
                    adotx = dot_product(cntrs(1 : 3), x(p(1 : 3)))

                    theta = 0.0d0
                    if (adotx .gt. 0.0d0) then
                        nz = nz + 1

                        if (nz .gt. m_small) then
                            write(*, *) 'Too many non-null y. Increase m_small and re-run.'
                            stop
                        end if

                        theta = adotx / sum(winv(p(1 : 3))) 
            
                        ! Update dual variables
                        ind_z(nz) = l
                        z(nz) = theta 

                        x(p(1 : 3)) = x(p(1 : 3)) - theta * winv(p(1 : 3)) * cntrs(1 : 3)

                    end if
                    
                    ! Stop criteria 
                    stp = max(stp, abs(y_old - theta))

                end if
            end do
        end do
    end do

    do i = 1, n - 1
        l = l + 1

        if (nextnnyi .eq. - 1) then
            y_old = 0.0d0
        else
            if (l .lt. ind_y(nextnnyi)) then
               y_old = 0.0d0

            else !if (l .eq. ind_y(nextnnyi)) then
               y_old = y(nextnnyi) 
               nextnnyi = nextnnyi + 1

               if (nextnnyi .gt. ny) then
                  nextnnyi = - 1
               end if
            end if
        end if

       p(1 : 2) = (/i, n/) 

       cntrs(1 : 2) = (/weight(i), -1.0d0/)

       ! Correction step
       x(p(1 : 2)) = x(p(1 : 2)) + y_old * winv(p(1 : 2)) * cntrs(1 : 2)

       ! Projection step
       adotx = dot_product(cntrs(1 : 2), x(p(1 : 2)))
        
       theta = 0.0d0
       if (adotx .gt. weight(i) * a(i)) then
            nz = nz + 1

            if (nz .gt. m_small) then
                 write(*, *) 'Too many non-null y. Increase m_small and re-run.'
                stop
            end if

            theta = (adotx - weight(i) * a(i)) / (winv(i) * weight(i) ** 2 + winv(n))

            ! Update dual variables
            ind_z(nz) = l
            z(nz) = theta

            x(p(1 : 2)) = x(p(1 : 2)) - theta * winv(p(1 : 2)) * cntrs(1 : 2)

        end if

       ! Stop criteria 
       stp = max(stp, maxval(abs(y_old - theta) * cntrs(1 : 2)))

       l = l + 1

       if (nextnnyi .eq. - 1) then
            y_old = 0.0d0
        else
            if (l .lt. ind_y(nextnnyi)) then
            y_old = 0.0d0

            else !if (l .eq. ind_y(nextnnyi)) then
            y_old = y(nextnnyi) 
            nextnnyi = nextnnyi + 1

            if (nextnnyi .gt. ny) then
                nextnnyi = - 1
            end if
            end if
        end if

       cntrs(1 : 2) = (/- weight(i), -1.0d0/)

       ! Correction step
       x(p(1 : 2)) = x(p(1 : 2)) + y_old * winv(p(1 : 2)) * cntrs(1 : 2)

       ! Projection step
       adotx = dot_product(cntrs(1 : 2), x(p(1 : 2)))

       theta = 0.0d0
       if (adotx .gt. - weight(i) * a(i)) then          
            nz = nz + 1

            if (nz .gt. m_small) then
                write(*, *) 'Too many non-null y. Increase m_small and re-run.'
                stop
            end if

            theta = (adotx + weight(i) * a(i)) / (winv(i) * weight(i) ** 2 + winv(n))

            ! Update dual variables
            ind_z(nz) = l
            z(nz) = theta
            x(p(1 : 2)) = x(p(1 : 2)) - theta * winv(p(1 : 2)) * cntrs(1 : 2)

        end if

       ! Stop criteria 
       stp = max(stp, maxval(abs(y_old - theta) * cntrs(1 : 2)))

    end do

    !!!! End of the update !!!!

    if (ny .eq. 0) then
        nextnnyi = - 1
    else
        nextnnyi = 1
    end if

    l = 0

    cmpl = 0.0d0
    max_inf = 0.0d0
    aypc(1 : n - 1) = 0.0d0
    aypc(n) = 1.0d0
    by = 0.0d0

    do i = 1, n_real - 1
        do j = i + 1, n_real
            do k = 1, n_real
                if (k .ne. i .and. k .ne. j) then
                    l = l + 1

                    if (nextnnyi .eq. - 1) then
                        y_old = 0.0d0
                     else
                        if (l .lt. ind_y(nextnnyi)) then
                           y_old = 0.0d0

                        else !if (l .eq. ind_y(nextnnyi)) then
                           y_old = z(nextnnyi)
                           nextnnyi = nextnnyi + 1

                           if (nextnnyi .gt. nz) then
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

                    cntrs = (/ 1.0d0, - 1.0d0, - 1.0d0 /)

                    adotx = dot_product(cntrs(1 : 3), x(p(1 : 3)))

                    max_inf = max(max_inf, adotx)

                    cmpl = max(cmpl, min(-(adotx), y_old))


                    if (y_old .ne. 0.0d0) then 

                        aypc(p(1 : 3)) = aypc(p(1 : 3)) + cntrs(1 : 3) * y_old

                    end if

                end if
            end do
        end do
    end do

    do i = 1, n - 1
        l = l + 1

        if (nextnnyi .eq. - 1) then
            y_old = 0.0d0
            else
                if (l .lt. ind_y(nextnnyi)) then
                    y_old = 0.0d0

                else !if (l .eq. ind_y(nextnnyi)) then
                    y_old = z(nextnnyi)
                    nextnnyi = nextnnyi + 1

                    if (nextnnyi .gt. nz) then
                        nextnnyi = - 1 
                    end if
                end if
        end if

       p(1 : 2) = (/i, n/) 
        
       cntrs(1 : 2) = (/weight(i), -1.0d0/)

       adotx = dot_product(cntrs(1 : 2), x(p(1 : 2)))

       max_inf = max(max_inf, adotx - weight(i) * a(i))

       cmpl = max(cmpl, min(-(adotx - weight(i) * a(i)), y_old))

       if (y_old .ne. 0.0d0) then 
            
            aypc(p(1 : 2)) = aypc(p(1 : 2)) + cntrs(1 : 2) * y_old

            by = by + weight(i) * a(i) * y_old
            
       end if

       l = l + 1

       if (nextnnyi .eq. - 1) then
            y_old = 0.0d0
            else
                if (l .lt. ind_y(nextnnyi)) then
                    y_old = 0.0d0

                else !if (l .eq. ind_y(nextnnyi)) then
                    y_old = z(nextnnyi)
                    nextnnyi = nextnnyi + 1

                    if (nextnnyi .gt. nz) then
                        nextnnyi = - 1 
                    end if
                end if
        end if

       cntrs(1 : 2) = (/- weight(i), -1.0d0/)

       adotx = dot_product(cntrs(1 : 2), x(p(1 : 2)))

       max_inf = max(max_inf, adotx + weight(i) * a(i))

       cmpl = max(cmpl, min(-(adotx + weight(i) * a(i)), y_old))

       if (y_old .ne. 0.0d0) then 

            aypc(p(1 : 2)) = aypc(p(1 : 2)) + cntrs(1 : 2) * y_old

            by = by - weight(i) * a(i) * y_old

       end if

    end do

    ! c^T * x
    linear_of = x(n)
    ! c^T * x + (1 / 2 * gamma) * x^T * W * x
    quadratic_of = (0.5d0 / gamma) * sum(x(1 : n) ** 2) + linear_of
    ! - b^T * y - (1 / 2 * gamma) * (A^T * y + c) * W^-1 *  (A^T * y + c)
    dual_of = - by - (0.5d0 * gamma) * sum(aypc(1 : n) ** 2)

    ! A^T * y + c + gamma * W * x
    optm = maxval(abs(aypc(1 : n) + (1.0d0 / gamma) * x(1 : n)))

    write(*, *) "At", iter, "inf=", max_inf, 'compl=', cmpl, 'optimality=', optm, 'gap=', (dual_of - quadratic_of) / dual_of

   
    if (stp .le. tol .and. max_inf .le. tol) then
    
        call CPU_TIME(final_time)

        call write_parameters(inputfnm, iter, final_time - initial_time, stp, max_inf, cmpl, optm,  &
                                (dual_of - quadratic_of) / dual_of,linear_of, gamma)
        go to 500
    end if

    if (iter .gt. max_iter) then
        print *, 'Reached maximum iterations without convergence.'
        
        call CPU_TIME(final_time)

        call write_parameters(inputfnm, iter, final_time - initial_time, stp, max_inf, cmpl, optm,  &
                                (dual_of - quadratic_of) / dual_of, linear_of, gamma)

        go to 500
    endif        

    ny = nz
    y(1 : ny) = z(1 : nz)
    ind_y(1 : ny) = ind_z(1 : nz)

    go to 100

500 continue  

    deallocate(winv, x, aypc, a, weight, y, ind_y, z, ind_z, stat = allocerr)
    
    if (allocerr .ne. 0) then
        write(*, *) 'Deallocation error in main program.'
        stop
    end if

contains

    subroutine index_vec(i, j, k)
        use iso_fortran_env
        implicit none 

        !Scalar arguments
        integer(int32), intent(in)  :: i, j
        integer(int32), intent(out) :: k
    
        k = ((j - 1) * j) / 2 - j + i + 1
    
    end subroutine index_vec  

    subroutine write_parameters(file, n_iter, time, stp, max_inf, cmpl, optm, gap, linear_of, gamma)
        use iso_fortran_env
        implicit none

        !Scalar arguments
        integer(int32)    :: unit_number, size1, n_iter
        real(real64)      :: max_inf, cmpl, optm, gap, stp,  gamma, time, linear_of

        !Character arguments
        character(len = 80) :: file
        character(len = 80) :: outputfnm

        !Prompt for the output file
        write(*,*) 'Enter the output file: '
        read(*,'(A)') outputfnm
        !outputfnm = 'dykstra_linf.txt'

        unit_number = 10

        open(unit = unit_number, file = outputfnm, action = "write", position = "append", status = "unknown")

        ! Check if the file is empty and add headers
        inquire(file = outputfnm, size = size1)
        if (size1 == 0) then
            write(unit_number, '(A)') "File       | Iterations | CPU Time | Stop Criteria |" // &
                                    " Max Infeasibility | Complementarity |" // &
                                    " Optimality | Dual Gap | Linear function | Gamma"
        end if

        write(unit_number, '(A10,    6X,    I8,     3X,     1P,F8.2,    8X,     1P,E8.1,     12X, &
           & 1P,E8.1,     10X,      1P,E8.1,    5X,    1P,E8.1,    3X,       1P,E8.1,    6X,   &
           & 0P,F12.6,    3X,   0P,F5.1)') &
            trim(file), n_iter, time, stp, max_inf, cmpl, optm, gap, linear_of, gamma
           

        close(unit_number)
    end subroutine write_parameters
    
end program dykstra_solution_l1
  
