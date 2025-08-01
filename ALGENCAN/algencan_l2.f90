! ******************************************************************
! ******************************************************************

program algencan_solution_l2

    implicit none
 
    !Common blocks
    integer               :: m
    integer, pointer      :: triples(:, :)
    real(kind = 8), pointer :: a(:), w(:)
    common /commonblock/ a, w, triples, m
    
    ! Local scalars
    logical      :: checkder
    integer      :: allocerr, hnnzmax, inform, jcnnzmax, n, nvparam, j, i, k, cont, n_real 
    real(kind = 8) :: cnorm, efacc, efstain, eoacc, eostain, epsfeas, epsopt, &
                    f, nlpsupn, snorm, element1, element2, s, final_time, initial_time
 
    ! Local subroutines
    character(len = 15)     :: strtmp
    character(len = 80)     :: specfnm, outputfnm, vparam(10)
    logical               :: coded(11)
    logical,      pointer :: equatn(:), linear(:)
    real(kind = 8), pointer :: l(:), lambda(:), u(:), x(:)
    
    ! External subroutines
    external :: myevalf, myevalg, myevalh, myevalc, myevaljac, myevalhc, &
                myevalfc, myevalgjac, myevalgjacp, myevalhl, myevalhlp
 
    ! Read input data
    character(len = 80) :: inputfnm
    write(*, *) 'Enter the input file: '
    read(*, '(A)') inputfnm  

    open(10, file = trim(adjustl(inputfnm)), status = "old", action = "read")
                
    ! Read the value of n_real from the file
    read(10, *) n_real
 
    m = (n_real * (n_real - 1)) / 2 * (n_real - 2)
    n = n_real * (n_real - 1) / 2
 
    ! Allocate matrices A and W
    allocate(a(n), w(n), triples(3, m), stat = allocerr)

    ! Check if allocation succeeded
    if (allocerr .ne. 0) then
        write(*, *) 'Allocation error in main program (1).'
        stop
    end if
    
    cont = 0
    do i = 1, n_real - 1
       do j = i + 1, n_real
          do k = 1, n_real
             if ( k .ne. i .and. k .ne. j ) then
                cont = cont + 1
                if ( cont .gt. m ) then
                   write(*, *) 'There is something wrong here (1).'
                   stop
                end if
 
                triples(1, cont) = i
                triples(2, cont) = j
                triples(3, cont) = k
             end if
          end do
       end do
    end do
 
    if ( cont .ne. m ) then
       write(*, *) 'There is something wrong here (2).'
       stop
    end if
 
    ! Read matrices A and W from the file
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
        write(*, *) 'There is something wrong here (3).'
        stop
    end if

 
    close(10)
 
    ! Set lower bounds, upper bounds, and initial guess
    allocate(x(n), l(n), u(n), stat = allocerr)

    if ( allocerr .ne. 0 ) then
       write(*,*) 'Allocation error in main program (2).'
       stop
    end if
 
    l(1 : n) =   0.0d0
    u(1 : n) =   1.0d+20
 
    ! Iniciate x = a
    x(1 : n) = a(1 : n)
 
    allocate(equatn(m), linear(m), lambda(m), stat = allocerr)

    if ( allocerr .ne. 0 ) then
       write(*, *) 'Allocation error in main program (3).'
       stop
    end if
 
    equatn(1 : m) = .false.
    lambda(1 : m) = 0.0d0

    linear(1 : m) = .true.
 
    ! Coded subroutines
    coded(1 : 6)  = .true.  ! fsub, gsub, hsub, csub, jacsub, hcsub
    coded(7 : 11) = .false. ! fcsub,gjacsub,gjacpsub,hlsub,hlpsub
 
    ! Upper bounds on the number of sparse-matrices non-null elements
    jcnnzmax = 3 * m 
    hnnzmax  = n
 
    ! Checking derivatives?
    checkder = .false.
 
    ! Parameters setting
    epsfeas   = 1.0d-08
    epsopt    = 1.0d-08
 
    efstain   = sqrt(epsfeas)
    eostain   = epsopt ** 1.5d0
 
    efacc     = sqrt(epsfeas)
    eoacc     = sqrt(epsopt)
 
    outputfnm = ''
    specfnm   = ''
 
    nvparam = 2
    vparam(1) = 'TRUNCATED-NEWTON-LINE-SEARCH-INNER-SOLVER'
    vparam(2) = 'SKIP-ACCELERATION-PROCESS'
 
    call CPU_TIME(initial_time)

    ! Optimize
    call algencan(myevalf, myevalg, myevalh, myevalc, myevaljac, myevalhc, &
                 myevalfc, myevalgjac, myevalgjacp, myevalhl, myevalhlp, jcnnzmax, &
                 hnnzmax, epsfeas, epsopt, efstain, eostain, efacc, eoacc, outputfnm, &
                 specfnm, nvparam, vparam, n, x, l, u, m, lambda, equatn, linear, coded, &
                 checkder, f, cnorm, snorm, nlpsupn, inform)

    s = norm2(w *(x - a))   

    call CPU_TIME(final_time)

    call write_results(inputfnm, n_real, m, final_time - initial_time, s)
 
    deallocate(x, l, u, lambda, equatn, linear, a, w, stat = allocerr)

    if ( allocerr .ne. 0 ) then
        write(*, *) 'Deallocation error in main program.'
         stop
    end if
 
    stop
 
end program algencan_solution_l2
 
subroutine myevalf(n, x, f, flag)
   
    implicit none
 
    ! Scalar arguments
    integer,      intent(in)  :: n
    integer,      intent(out) :: flag
    real(kind = 8), intent(out) :: f
 
    ! Array arguments
    real(kind = 8), intent(in) :: x(n)

    integer               :: n_real, n_triples
    integer, pointer      :: triples(:, :)
    real(kind = 8), pointer :: a(:), w(:)
    common /commonblock/ a, w, triples, n_real, n_triples
 
    flag = 0
    
    f = 0.5d0 * sum(w(1 : n) ** 2 * (x(1 : n) - a(1 : n)) ** 2)
   
end subroutine myevalf
 
 
subroutine myevalg(n, x, g, flag)
 
    implicit none
 
    ! Scalar arguments
    integer, intent(in)  :: n
    integer, intent(out) :: flag
 
    ! Array arguments
    real(kind = 8), intent(in)  :: x(n)
    real(kind = 8), intent(out) :: g(n)

    ! Common blocks
    integer               :: m
    integer, pointer      :: triples(:, :)
    real(kind = 8), pointer :: a(:), w(:)
    common /commonblock/ a, w, triples, m

    flag = 0
    
    g(1 : n) = (x(1 : n) - a(1 : n)) * w(1 : n) ** 2
 
end subroutine myevalg
 
 
subroutine myevalh(n, x, hrow, hcol, hval, hnnz, lim, lmem, flag)
 
    implicit none
 
    ! Scalar arguments
    logical,      intent(out) :: lmem
    integer,      intent(in)  :: lim, n
    integer,      intent(out) :: flag, hnnz
 
    ! Array arguments
    integer, intent(out)        :: hcol(lim), hrow(lim)
    real(kind = 8), intent(in)  :: x(n)
    real(kind = 8), intent(out) :: hval(lim)

    ! Common blocks
    integer               :: m
    integer, pointer      :: triples(:, :)
    real(kind = 8), pointer :: a(:), w(:)
    common /commonblock/ a, w, triples, m

    integer :: i
 
    flag = 0
    lmem = .false.
 
    hnnz = n

    if ( hnnz .gt. lim ) then
        lmem = .true.
        return
    end if

    hcol(1 : n) = (/ (i, i = 1, n) /)
    hrow(1 : n) = hcol(1 : n)
    hval(1 : n) = w(1 : n) ** 2 

end subroutine myevalh
 
 
subroutine myevalc(n, x, ind, c, flag)
 
    implicit none
 
    ! Scalar arguments
    integer                   :: i, j, k, p_1, p_2, p_3
    integer, intent(in)       :: ind, n
    integer, intent(out)      :: flag
    real(kind = 8), intent(out) :: c
 
    ! Array arguments
    real(kind = 8), intent(in)  :: x(n)
 
    ! Common blocks
    integer               :: m
    integer, pointer      :: triples(:, :)
    real(kind = 8), pointer :: a(:), w(:)
    common /commonblock/ a, w, triples, m
    
    flag = 0
 
    if ( 1 .le. ind .and. ind .le. m ) then
       
       i = triples(1, ind)
       j = triples(2, ind)
       k = triples(3, ind)
 
       if ( k .lt. i .and. i .lt. j ) then
          ! xij <= xki + xkj
          call index_vec(i, j, p_1)
          call index_vec(k, i, p_2)
          call index_vec(k, j, p_3)
 
       else if ( i .lt. k .and. k .lt. j ) then
          ! xij <= xik + xkj
          call index_vec(i, j, p_1)
          call index_vec(i, k, p_2)
          call index_vec(k, j, p_3)
          
       else if ( i .lt. j .and. j .lt. k ) then
          ! xij <= xik + xjk
          call index_vec(i, j, p_1)
          call index_vec(i, k, p_2)
          call index_vec(j, k, p_3)
       else
          write(*, *) 'There is something wrong in myevalc.'
          stop
       end if
       
       c = x(p_1) - x(p_2) - x(p_3)
 
    end if 
end subroutine myevalc
 
 
subroutine myevaljac(n, x, ind, jcvar, jcval, jcnnz, lim, lmem, flag)
 
    implicit none

    ! Scalar arguments
    logical, intent(out) :: lmem
    integer, intent(in)  :: ind, lim, n
    integer, intent(out) :: flag, jcnnz
 
    ! Array  arguments
    integer, intent(out)      :: jcvar(lim)
    real(kind = 8), intent(in)  :: x(n)
    real(kind = 8), intent(out) :: jcval(lim)

    integer :: i, j, k, p_1, p_2, p_3
 
    ! Common blocks
    integer               :: m
    integer, pointer      :: triples(:, :)
    real(kind = 8), pointer :: a(:), w(:)
    common /commonblock/ a, w, triples, m
    
    flag = 0
    lmem = .false.
 
    if ( 1 .le. ind .and. ind .le. m ) then
        jcnnz = 3

        if ( jcnnz .gt. lim ) then
            lmem = .true.
            return
        end if

        i = triples(1, ind)
        j = triples(2, ind)
        k = triples(3, ind)

        if ( k .lt. i .and. i .lt. j ) then
            ! xij <= xki + xkj
            call index_vec(i, j, p_1)
            call index_vec(k, i, p_2)
            call index_vec(k, j, p_3)

        else if ( i .lt. k .and. k .lt. j ) then
            ! xij <= xik + xkj
            call index_vec(i, j, p_1)
            call index_vec(i, k, p_2)
            call index_vec(k, j, p_3)
            
        else if ( i .lt. j .and. j .lt. k ) then
            ! xij <= xik + xjk
            call index_vec(i, j, p_1)
            call index_vec(i, k, p_2)
            call index_vec(j, k, p_3)
        else
            write(*,*) 'There is something wrong here in myevaljac.'
            stop
        end if

        jcvar(1 : 3) = (/ p_1, p_2, p_3 /)
        jcval(1 : 3) = (/ 1.0d0, -1.0d0, -1.0d0 /)

    end if 
 
end subroutine myevaljac
 
 
subroutine myevalhc(n, x, ind, hcrow, hccol, hcval, hcnnz, lim, lmem, flag)
 
    implicit none

    ! Scalar arguments
    logical, intent(out) :: lmem
    integer, intent(in)  :: ind, lim, n
    integer, intent(out) :: flag, hcnnz
 
    ! Array arguments
    integer, intent(out)      :: hccol(lim), hcrow(lim)
    real(kind = 8), intent(in)  :: x(n)
    real(kind = 8), intent(out) :: hcval(lim)
 
 
    flag = 0
    lmem = .false.  
 
    hcnnz = 0.0d0
        
end subroutine myevalhc
 
 
subroutine myevalfc(n, x, f, m, c, flag)
 
    implicit none
 
    ! Scalar arguments
    integer, intent(in)       :: m, n
    integer, intent(out)      :: flag
    real(kind = 8), intent(out) :: f
 
    ! Array arguments
    real(kind = 8), intent(in)  :: x(n)
    real(kind = 8), intent(out) :: c(m)
 
    flag = - 1
 
end subroutine myevalfc
 
 
subroutine myevalgjac(n, x, g, m, jcfun, jcvar, jcval, jcnnz, lim, lmem, flag)
 
    implicit none
 
    ! Scalar arguments
    logical, intent(out) :: lmem
    integer, intent(in)  :: lim, m, n
    integer, intent(out) :: flag, jcnnz
 
    ! Array arguments
    integer,      intent(out) :: jcfun(lim), jcvar(lim)
    real(kind = 8), intent(in)  :: x(n)
    real(kind = 8), intent(out) :: g(n), jcval(lim)
 
    flag = - 1
 
end subroutine myevalgjac
 
 
subroutine myevalgjacp(n, x, g, m, p, q, work, gotj, flag)
 
    implicit none
 
    ! Scalar arguments
    logical,   intent(inout) :: gotj
    integer,   intent(in)    :: m, n
    integer,   intent(out)   :: flag
    character, intent(in)    :: work
 
    ! Array arguments
    real(kind = 8), intent(in)    :: x(n)
    real(kind = 8), intent(inout) :: p(m), q(n)
    real(kind = 8), intent(out)   :: g(n)
 
    flag = - 1
 
end subroutine myevalgjacp
 
 
subroutine myevalhl(n, x, m, lambda, sf, sc, hlrow, hlcol, hlval, hlnnz, lim, lmem, flag)
 
    implicit none
 
    ! Scalar arguments
    logical, intent(out)     :: lmem
    integer, intent(in)      :: lim, m, n
    integer, intent(out)     :: flag, hlnnz
    real(kind = 8), intent(in) :: sf
 
    ! Array arguments
    integer, intent(out)      :: hlcol(lim), hlrow(lim)
    real(kind = 8), intent(in)  :: lambda(m), sc(m), x(n)
    real(kind = 8), intent(out) :: hlval(lim)
 
    flag = - 1
 
end subroutine myevalhl
 
 
subroutine myevalhlp(n, x, m, lambda, sf, sc, p, hp, goth, flag)
 
    implicit none
 
    ! Scalar arguments
    logical, intent(inout)   :: goth
    integer, intent(in)      :: m, n
    integer, intent(out)     :: flag
    real(kind = 8), intent(in) :: sf
 
    ! Array arguments
    real(kind = 8), intent(in)  :: lambda(m), p(n), sc(m), x(n)
    real(kind = 8), intent(out) :: hp(n)
 
    flag = - 1
 
end subroutine myevalhlp
 
subroutine index_vec(i, j, k)

    implicit none 

    integer, intent(in)  :: i, j
    integer, intent(out) :: k
 
    k = ((j - 1) * j) / 2 - j + i + 1
    return 
 
end subroutine index_vec
 
subroutine write_results(file, n_real, m, final_time, objfnc)

    implicit none

    character(len = 80) :: file
    integer             :: unit_number, size1, m, n_real
    real(kind = 8)      :: final_time, objfnc
    character(len = 80) :: outputfnm
 
    write(*, *) 'Enter the output file: '
    read(*, '(A)') outputfnm
 
    unit_number = 10
 
    open(unit = unit_number, file = outputfnm, action = "write", position = "append", status = "unknown")
 
    inquire(file = outputfnm, size = size1)
    if (size1 == 0) then
        write(unit_number, '(A)') "File       | Size of A | Number of constraints |" // &
                              "    CPU Time    | Objective Function"
    end if
 
    write(unit_number, '(A10, 3X, I10, 3X, I24, 3X, F12.4, 3X, F18.6)') file, n_real, m, final_time, objfnc
 
    close(unit_number)
end subroutine write_results
   
