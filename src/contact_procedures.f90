module contact_procedures

  use params
  use contact_object

  implicit none

contains

  function method_required_link_params(a_p) result(valid_params)

    implicit none

    type(analysis_params), intent(in) :: a_p

    integer(int32) :: N_required

    logical :: valid_params

    valid_params = .true.

    N_required = -1

    select case (a_p%method)

    case(1)

       N_required = 1

    case(2)

       N_required = 1
       
    case(3)

       N_required = 1
       
    case(4)

       N_required = 1

    case(6)

       N_required = 2

    case(7)

       N_required = 3

    case(8)

       N_required = 2

    end select

    if (a_p%N_link_params .ne. N_required) then
       write(unit_log,"(A,I2,A,I2,A)")' ERROR: method #',a_p%method,&
            ' requires ',N_required,&
            ' parameters'
       valid_params = .false.
    end if

  end function method_required_link_params

  subroutine contact_contribution(c_m,a_p,x)

    implicit none

    type(contact_mat), intent(inout) :: c_m
    type(analysis_params), intent(in) :: a_p
    real(rp), intent(in) :: x(3,a_p%N)

    select case (a_p%method)

    case(1)

       call method_1(c_m,a_p,x)

    case(2)
       
    case(3)

       call method_3(c_m,a_p,x)
       
    case(4)

       call method_4(c_m,a_p,x)

    case(6)

       call method_6(c_m,a_p,x)

    case(7)

       call method_7(c_m,a_p,x)

    case(8)

       call method_8(c_m,a_p,x)

    end select

  end subroutine contact_contribution

  !distance matrix between CoMs of CG beads
  subroutine method_1(c_m,a_p,x)

    implicit none

    type(contact_mat), intent(inout) :: c_m
    type(analysis_params), intent(in) :: a_p
    real(rp), intent(in) :: x(3,a_p%N)

    integer(int32) :: i_pair, b1, b2

    real(rp) :: dx(1:3)

    if (a_p%run_flag .eq. 1) then
    
       do i_pair = 1,c_m%N_pair

          b1 = c_m%pairs(1,i_pair)
          b2 = c_m%pairs(2,i_pair)

          if (b1 .ne. b2) then

             dx = sum(x(:,c_m%map(1,b1):c_m%map(2,b1)))/c_m%N_mono(b1)
             dx = dx - sum(x(:,c_m%map(1,b2):c_m%map(2,b2)))/c_m%N_mono(b2)

             c_m%rval(i_pair) = c_m%ival(i_pair) + norm2(dx)

          end if

       end do

    elseif (a_p%run_flag .eq. 2) then

       !$OMP PARALLEL DO &
       !$OMP DEFAULT(PRIVATE) &
       !$OMP SHARED(c_m,a_p,x)

       do i_pair = 1,c_m%N_pair

          b1 = c_m%pairs(1,i_pair)
          b2 = c_m%pairs(2,i_pair)

          if (b1 .ne. b2) then

             dx = sum(x(:,c_m%map(1,b1):c_m%map(2,b1)))/c_m%N_mono(b1)
             dx = dx - sum(x(:,c_m%map(1,b2):c_m%map(2,b2)))/c_m%N_mono(b2)

             c_m%rval(i_pair) = c_m%ival(i_pair) + norm2(dx)

          end if

       end do

    end if

  end subroutine method_1

  !contact count of CoMs of CG beads less than r_cut
  subroutine method_3(c_m,a_p,x)

    implicit none

    type(contact_mat), intent(inout) :: c_m
    type(analysis_params), intent(in) :: a_p
    real(rp), intent(in) :: x(3,a_p%N)

    integer(int32) :: i_pair, b1, b2

    real(rp) :: dx(1:3)

    if (a_p%run_flag .eq. 1) then
    
       do i_pair = 1,c_m%N_pair

          b1 = c_m%pairs(1,i_pair)
          b2 = c_m%pairs(2,i_pair)

          dx = sum(x(:,c_m%map(1,b1):c_m%map(2,b1)))/c_m%N_mono(b1)
          dx = dx - sum(x(:,c_m%map(1,b2):c_m%map(2,b2)))/c_m%N_mono(b2)

          if (norm2(dx) .lt. a_p%link_params(1)) c_m%ival(i_pair) = c_m%ival(i_pair) + 1

       end do

    elseif (a_p%run_flag .eq. 2) then

       !$OMP PARALLEL DO &
       !$OMP DEFAULT(PRIVATE) &
       !$OMP SHARED(c_m,a_p,x)

       do i_pair = 1,c_m%N_pair

          b1 = c_m%pairs(1,i_pair)
          b2 = c_m%pairs(2,i_pair)

          dx = sum(x(:,c_m%map(1,b1):c_m%map(2,b1)))/c_m%N_mono(b1)
          dx = dx - sum(x(:,c_m%map(1,b2):c_m%map(2,b2)))/c_m%N_mono(b2)

          if (norm2(dx) .lt. a_p%link_params(1)) c_m%ival(i_pair) = c_m%ival(i_pair) + 1

       end do

    end if

  end subroutine method_3


  !contact count of CoMs of CG beads less than r_cut
  subroutine method_4(c_m,a_p,x)

    implicit none

    type(contact_mat), intent(inout) :: c_m
    type(analysis_params), intent(in) :: a_p
    real(rp), intent(in) :: x(3,a_p%N)

    integer(int32) :: i_pair, b1, b2, j1, j2

    real(rp) :: dx(1:3)

    !if (norm2(dx) .lt. a_p%link_params(1)) c_m%ival(i_pair) = c_m%ival(i_pair) + 1

    if (a_p%run_flag .eq. 1) then
    
       do i_pair = 1,c_m%N_pair

          b1 = c_m%pairs(1,i_pair)
          b2 = c_m%pairs(2,i_pair)

          if (b1 .ne. b2) then

             do j1 = c_m%map(1,b1),c_m%map(2,b1)

                do j2 = c_m%map(1,b2),c_m%map(2,b2)

                   dx = x(:,j1) - x(:,j2)

                   if (norm2(dx) .lt. a_p%link_params(1)) &
                        c_m%ival(i_pair) = c_m%ival(i_pair) + 1

                end do

             end do

          else

             do j1 = c_m%map(1,b1),c_m%map(2,b1)

                do j2 = c_m%map(1,b2),c_m%map(2,b2)

                   if (j1 .lt. j2) then
                   
                      dx = x(:,j1) - x(:,j2)

                      if (norm2(dx) .lt. a_p%link_params(1)) &
                           c_m%ival(i_pair) = c_m%ival(i_pair) + 1

                   end if

                end do

             end do

          end if

       end do

    elseif (a_p%run_flag .eq. 2) then

       !$OMP PARALLEL DO &
       !$OMP DEFAULT(PRIVATE) &
       !$OMP SHARED(c_m,a_p,x)

       do i_pair = 1,c_m%N_pair

          b1 = c_m%pairs(1,i_pair)
          b2 = c_m%pairs(2,i_pair)

          if (b1 .ne. b2) then

             do j1 = c_m%map(1,b1),c_m%map(2,b1)

                do j2 = c_m%map(1,b2),c_m%map(2,b2)

                   dx = x(:,j1) - x(:,j2)

                   if (norm2(dx) .lt. a_p%link_params(1)) &
                        c_m%ival(i_pair) = c_m%ival(i_pair) + 1

                end do

             end do

          else

             do j1 = c_m%map(1,b1),c_m%map(2,b1)

                do j2 = c_m%map(1,b2),c_m%map(2,b2)

                   if (j1 .lt. j2) then
                   
                      dx = x(:,j1) - x(:,j2)

                      if (norm2(dx) .lt. a_p%link_params(1)) &
                           c_m%ival(i_pair) = c_m%ival(i_pair) + 1

                   end if

                end do

             end do

          end if

       end do

    end if

  end subroutine method_4

  !contact count of CoMs of CG beads less than r_cut with correction for backbone
  subroutine method_5(c_m,a_p,x)

    implicit none

    type(contact_mat), intent(inout) :: c_m
    type(analysis_params), intent(in) :: a_p
    real(rp), intent(in) :: x(3,a_p%N)

    integer(int32) :: i_pair, b1, b2, j1, j2, excl

    real(rp) :: dx(1:3)

    excl = int(a_p%link_params(1)/2)

    if (a_p%run_flag .eq. 1) then
    
       do i_pair = 1,c_m%N_pair

          b1 = c_m%pairs(1,i_pair)
          b2 = c_m%pairs(2,i_pair)

          do j1 = c_m%map(1,b1),c_m%map(2,b1)

             do j2 = c_m%map(1,b2),c_m%map(2,b2)

                if (min_dist(a_p%N,j1,j2) .gt. excl) then

                   dx = x(:,j1) - x(:,j2)

                   if (norm2(dx) .lt. a_p%link_params(1)) c_m%ival(i_pair) = c_m%ival(i_pair) + 1

                end if

             end do

          end do

       end do

    elseif (a_p%run_flag .eq. 2) then

       !$OMP PARALLEL DO &
       !$OMP DEFAULT(PRIVATE) &
       !$OMP SHARED(c_m,a_p,x,excl)

       do i_pair = 1,c_m%N_pair

          b1 = c_m%pairs(1,i_pair)
          b2 = c_m%pairs(2,i_pair)

          do j1 = c_m%map(1,b1),c_m%map(2,b1)

             do j2 = c_m%map(1,b2),c_m%map(2,b2)

                if (min_dist(a_p%N,j1,j2) .gt. excl) then

                   dx = x(:,j1) - x(:,j2)

                   if (norm2(dx) .lt. a_p%link_params(1)) c_m%ival(i_pair) = c_m%ival(i_pair) + 1

                end if

             end do

          end do

       end do

    end if

  end subroutine method_5

  !sum of probabilities P=(sigma/d)^n of constituent beads
  subroutine method_6(c_m,a_p,x)

    implicit none

    type(contact_mat), intent(inout) :: c_m
    type(analysis_params), intent(in) :: a_p
    real(rp), intent(in) :: x(3,a_p%N)

    integer(int32) :: i_pair, b1, b2, j1, j2

    real(rp) :: dx(1:3)

    if (a_p%run_flag .eq. 1) then
    
       do i_pair = 1,c_m%N_pair

          b1 = c_m%pairs(1,i_pair)
          b2 = c_m%pairs(2,i_pair)

          if (b1 .ne. b2) then

             do j1 = c_m%map(1,b1),c_m%map(2,b1)

                do j2 = c_m%map(1,b2),c_m%map(2,b2)

                   dx = x(:,j1) - x(:,j2)

                   c_m%rval(i_pair) = c_m%rval(i_pair) + &
                        (a_p%link_params(1)/norm2(dx))**a_p%link_params(2)

                end do

             end do

          else

             do j1 = c_m%map(1,b1),c_m%map(2,b1)

                do j2 = c_m%map(1,b2),c_m%map(2,b2)

                   if (j1 .lt. j2) then
                   
                      dx = x(:,j1) - x(:,j2)

                      c_m%rval(i_pair) = c_m%rval(i_pair) + &
                           (a_p%link_params(1)/norm2(dx))**a_p%link_params(2)

                   end if

                end do

             end do

          end if

       end do

    elseif (a_p%run_flag .eq. 2) then

       !$OMP PARALLEL DO &
       !$OMP DEFAULT(PRIVATE) &
       !$OMP SHARED(c_m,a_p,x)

       do i_pair = 1,c_m%N_pair

          b1 = c_m%pairs(1,i_pair)
          b2 = c_m%pairs(2,i_pair)

          if (b1 .ne. b2) then

             do j1 = c_m%map(1,b1),c_m%map(2,b1)

                do j2 = c_m%map(1,b2),c_m%map(2,b2)

                   dx = x(:,j1) - x(:,j2)

                   c_m%rval(i_pair) = c_m%rval(i_pair) + &
                        (a_p%link_params(1)/norm2(dx))**a_p%link_params(2)

                end do

             end do

          else

             do j1 = c_m%map(1,b1),c_m%map(2,b1)

                do j2 = c_m%map(1,b2),c_m%map(2,b2)

                   if (j1 .lt. j2) then
                   
                      dx = x(:,j1) - x(:,j2)

                      c_m%rval(i_pair) = c_m%rval(i_pair) + &
                           (a_p%link_params(1)/norm2(dx))**a_p%link_params(2)

                   end if

                end do

             end do

          end if

       end do

    end if

  end subroutine method_6

  !hill-type function
  subroutine method_7(c_m,a_p,x)

    implicit none

    type(contact_mat), intent(inout) :: c_m
    type(analysis_params), intent(in) :: a_p
    real(rp), intent(in) :: x(3,a_p%N)

    integer(int32) :: i_pair, b1, b2, j1, j2

    real(rp) :: dx(1:3)

    if (a_p%run_flag .eq. 1) then
    
       do i_pair = 1,c_m%N_pair

          b1 = c_m%pairs(1,i_pair)
          b2 = c_m%pairs(2,i_pair)

          if (b1 .ne. b2) then

             do j1 = c_m%map(1,b1),c_m%map(2,b1)

                do j2 = c_m%map(1,b2),c_m%map(2,b2)

                   dx = x(:,j1) - x(:,j2)

                   c_m%rval(i_pair) = c_m%rval(i_pair) + &
                           hill_fxn(norm2(dx),&
                           a_p%link_params(1),&
                           a_p%link_params(2),&
                           a_p%link_params(3))

                end do

             end do

          else

             do j1 = c_m%map(1,b1),c_m%map(2,b1)

                do j2 = c_m%map(1,b2),c_m%map(2,b2)

                   if (j1 .lt. j2) then
                   
                      dx = x(:,j1) - x(:,j2)

                      c_m%rval(i_pair) = c_m%rval(i_pair) + &
                           hill_fxn(norm2(dx),&
                           a_p%link_params(1),&
                           a_p%link_params(2),&
                           a_p%link_params(3))

                   end if

                end do

             end do

          end if

       end do

    elseif (a_p%run_flag .eq. 2) then

       !$OMP PARALLEL DO &
       !$OMP DEFAULT(PRIVATE) &
       !$OMP SHARED(c_m,a_p,x)

       do i_pair = 1,c_m%N_pair

          b1 = c_m%pairs(1,i_pair)
          b2 = c_m%pairs(2,i_pair)

          if (b1 .ne. b2) then

             do j1 = c_m%map(1,b1),c_m%map(2,b1)

                do j2 = c_m%map(1,b2),c_m%map(2,b2)

                   dx = x(:,j1) - x(:,j2)

                   c_m%rval(i_pair) = c_m%rval(i_pair) + &
                           hill_fxn(norm2(dx),&
                           a_p%link_params(1),&
                           a_p%link_params(2),&
                           a_p%link_params(3))

                end do

             end do

          else

             do j1 = c_m%map(1,b1),c_m%map(2,b1)

                do j2 = c_m%map(1,b2),c_m%map(2,b2)

                   if (j1 .lt. j2) then
                   
                      dx = x(:,j1) - x(:,j2)

                      c_m%rval(i_pair) = c_m%rval(i_pair) + &
                           hill_fxn(norm2(dx),&
                           a_p%link_params(1),&
                           a_p%link_params(2),&
                           a_p%link_params(3))

                   end if

                end do

             end do

          end if

       end do

    end if

  end subroutine method_7

  !tanh-type function
  subroutine method_8(c_m,a_p,x)

    implicit none

    type(contact_mat), intent(inout) :: c_m
    type(analysis_params), intent(in) :: a_p
    real(rp), intent(in) :: x(3,a_p%N)

    integer(int32) :: i_pair, b1, b2, j1, j2

    real(rp) :: dx(1:3)

    if (a_p%run_flag .eq. 1) then
    
       do i_pair = 1,c_m%N_pair

          b1 = c_m%pairs(1,i_pair)
          b2 = c_m%pairs(2,i_pair)

          if (b1 .ne. b2) then

             do j1 = c_m%map(1,b1),c_m%map(2,b1)

                do j2 = c_m%map(1,b2),c_m%map(2,b2)

                   dx = x(:,j1) - x(:,j2)

                   c_m%rval(i_pair) = c_m%rval(i_pair) + &
                           tanh_fxn(norm2(dx),&
                           a_p%link_params(1),&
                           a_p%link_params(2))

                end do

             end do

          else

             do j1 = c_m%map(1,b1),c_m%map(2,b1)

                do j2 = c_m%map(1,b2),c_m%map(2,b2)

                   if (j1 .lt. j2) then
                   
                      dx = x(:,j1) - x(:,j2)

                      c_m%rval(i_pair) = c_m%rval(i_pair) + &
                           tanh_fxn(norm2(dx),&
                           a_p%link_params(1),&
                           a_p%link_params(2))

                   end if

                end do

             end do

          end if

       end do

    elseif (a_p%run_flag .eq. 2) then

       !$OMP PARALLEL DO &
       !$OMP DEFAULT(PRIVATE) &
       !$OMP SHARED(c_m,a_p,x)

       do i_pair = 1,c_m%N_pair

          b1 = c_m%pairs(1,i_pair)
          b2 = c_m%pairs(2,i_pair)

          if (b1 .ne. b2) then

             do j1 = c_m%map(1,b1),c_m%map(2,b1)

                do j2 = c_m%map(1,b2),c_m%map(2,b2)

                   dx = x(:,j1) - x(:,j2)

                   c_m%rval(i_pair) = c_m%rval(i_pair) + &
                           tanh_fxn(norm2(dx),&
                           a_p%link_params(1),&
                           a_p%link_params(2))

                end do

             end do

          else

             do j1 = c_m%map(1,b1),c_m%map(2,b1)

                do j2 = c_m%map(1,b2),c_m%map(2,b2)

                   if (j1 .lt. j2) then
                   
                      dx = x(:,j1) - x(:,j2)

                      c_m%rval(i_pair) = c_m%rval(i_pair) + &
                           tanh_fxn(norm2(dx),&
                           a_p%link_params(1),&
                           a_p%link_params(2))
                      
                   end if

                end do

             end do

          end if

       end do

    end if

  end subroutine method_8

  !hill function for method 7
  pure function hill_fxn(x,x0,sigma,n) result(f)

    implicit none

    real(rp), intent(in) :: x, x0, sigma, n

    real(rp) :: f

    f = x - x0
    
    f = (sigma - x0)/f
    
    f = f**n

    f = f/(1.0+f)

  end function hill_fxn

  !tanh function for method 8
  pure function tanh_fxn(x,x0,sigma) result(f)

    implicit none

    real(rp), intent(in) :: x, x0, sigma

    real(rp) :: f

    f = (x - x0)/sigma

    f = 1.0+tanh(-f)

  end function tanh_fxn

  !function to calculate the minimum distance and directionality between two nodes on the ring
  pure function min_dist(N,j,k) result(d)

    implicit none
    
    integer(int32), intent(in) :: N, j, k
    integer(int32) :: half, d

    half=N/2

    if (k .ge. j) then
       if (k - j .ge. half) then
          d = j + n - k
       else
          d = k - j
       end if
    else
       if (j - k .ge. half) then
          d = k + n - j
       else
          d= j - k
       end if
    end if
    
  end function min_dist
  
end module
