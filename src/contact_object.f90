module contact_object

  use params

  implicit none

  !derived type to hold the analysis parameters
  type analysis_params
     integer(int32) :: run_flag, min_rep, max_rep, N_reps
     integer(int32) :: method, N, N_base, CG, map_flag, data_size, N_link_params
     real(rp) :: r
     real(rp), allocatable :: link_params(:)
     character(len=240) :: data_label, data_dir, map_file
  end type analysis_params

  !derived type to hold the symmetric pairwise matrix as it is populated
  type contact_mat
     integer(int32) :: N_pair, N_pair_base, N_CG, N_base_CG
     integer(int32), allocatable :: pairs(:,:), pairs_base(:,:), map(:,:), N_mono(:), base(:)
     integer(int32), allocatable :: ival(:), ival_base(:)
     real(rp), allocatable :: rval(:), rval_base(:)
  end type contact_mat

contains

  subroutine new_analysis_params(a_p,&
       run_flag,&
       min_rep,max_rep,&
       method,&
       N_link_params,&
       r,link_params,&
       CG,N,map_flag,&
       data_label,data_dir,map_file,&
       data_size)

    implicit none

    type(analysis_params), intent(out) :: a_p

    integer(int32), intent(in) :: run_flag, min_rep, max_rep
    integer(int32), intent(in) :: method, N, CG, map_flag, data_size, N_link_params

    real(rp), intent(in) :: r, link_params(1:N_link_params)

    character(len=240), intent(in) :: data_label, data_dir, map_file

    a_p%run_flag = run_flag
    a_p%min_rep = min_rep
    a_p%max_rep = max_rep
    a_p%N_reps = max_rep - min_rep + 1
    a_p%method = method
    a_p%N = N
    a_p%CG = CG
    a_p%map_flag = map_flag
    a_p%data_size = data_size
    
    a_p%r = r

    a_p%N_link_params = N_link_params

    if (allocated(a_p%link_params)) deallocate(a_p%link_params)
    allocate(a_p%link_params(1:N_link_params))
    a_p%link_params = link_params

    a_p%data_label = data_label
    a_p%data_dir = data_dir
    a_p%map_file = map_file

  end subroutine new_analysis_params

  subroutine new_contact_mat(c_m,a_p)

    implicit none

    type(contact_mat), intent(inout) :: c_m
    type(analysis_params), intent(in) :: a_p

    integer(int32) :: CG_rem, i_CG
    integer(int32) :: i, j, c

    if (allocated(c_m%pairs)) deallocate(c_m%pairs)
    if (allocated(c_m%map)) deallocate(c_m%map)
    if (allocated(c_m%base)) deallocate(c_m%base)
    if (allocated(c_m%N_mono)) deallocate(c_m%N_mono)
    if (allocated(c_m%ival)) deallocate(c_m%ival)
    if (allocated(c_m%rval)) deallocate(c_m%rval)
    if (allocated(c_m%ival_base)) deallocate(c_m%ival_base)
    if (allocated(c_m%rval_base)) deallocate(c_m%rval_base)

    !automatic coarse-graining if no information is provided
    if (a_p%CG .gt. 0) then
       CG_rem = modulo(a_p%N,a_p%CG)
       c_m%N_CG = a_p%N/a_p%CG
       if (CG_rem .gt. 0) c_m%N_CG = c_m%N_CG + 1
       
       !create mapping
       allocate(c_m%N_mono(1:c_m%N_CG))
       allocate(c_m%map(1:2,1:c_m%N_CG))

       if (CG_rem .gt. 0) then
          do i_CG=1,(c_m%N_CG-1)
             c_m%map(:,i_CG) = (/(i_CG-1)*a_p%CG+1,(i_CG)*a_p%CG/)
          end do
          c_m%map(:,c_m%N_CG) = (/(c_m%N_CG-1)*a_p%CG+1,(c_m%N_CG-1)*a_p%CG+CG_rem/)
       else
          do i_CG=1,c_m%N_CG
             c_m%map(:,i_CG) = (/(i_CG-1)*a_p%CG+1,(i_CG)*a_p%CG/)
          end do
       end if
       
    end if

    !calculate number of monomers per mapped CG locus
    
    c_m%N_mono = c_m%map(2,:) - c_m%map(1,:) + 1

    !create pairwise indices
    
    c_m%N_pair = c_m%N_CG*(c_m%N_CG-1)/2 + c_m%N_CG
    allocate(c_m%pairs(1:2,1:c_m%N_pair))

    c = 1
    do i=1,c_m%N_CG
       do j=i,c_m%N_CG
          c_m%pairs(:,c) = (/i,j/)
          c = c + 1
       end do
    end do

    !allocate value arrays (integer and real)

    allocate(c_m%ival(1:c_m%N_pair))
    allocate(c_m%rval(1:c_m%N_pair))

    c_m%ival = 0
    c_m%rval = 0.0d0

  end subroutine new_contact_mat

  !subroutine to map CG results back to base
  pure subroutine map_to_base(c_m)

    implicit none

    type(contact_mat), intent(inout) :: c_m

    integer(int32) :: i, j, i_pair

    real(rp) :: A_base(1:c_m%N_base_CG,1:c_m%N_base_CG)
    real(rp) :: A(1:c_m%N_CG,1:c_m%N_CG)

    call full_mat(A,c_m)

    A_base = A(1:c_m%N_base_CG,1:c_m%N_base_CG)

    do i=c_m%N_base_CG+1,c_m%N_CG

       A_base(c_m%base(i),c_m%base(i)) = A_base(c_m%base(i),c_m%base(i)) + A(i,i)

       do j=1,i-1

          A_base(c_m%base(i),c_m%base(j)) = A_base(c_m%base(i),c_m%base(j)) + A(i,j)
          A_base(c_m%base(j),c_m%base(i)) = A_base(c_m%base(j),c_m%base(i)) + A(i,j)
          
       end do
       
    end do

    call sym_mat_base(c_m,A_base)

  end subroutine map_to_base

  !subroutine to convert symmetric form to full matrix
  pure subroutine full_mat(A,c_m)

    implicit none

    type(contact_mat), intent(in) :: c_m
    real(rp), intent(out) :: A(1:c_m%N_CG,1:c_m%N_CG)

    integer(int32) :: i_pair, b1, b2

    do i_pair=1,c_m%N_pair

       b1 = c_m%pairs(1,i_pair)
       b2 = c_m%pairs(2,i_pair)

       if (b1 .eq. b2) then
          A(b1,b2) = c_m%rval(i_pair)
       else
          A(b1,b2) = c_m%rval(i_pair)
          A(b2,b1) = c_m%rval(i_pair)
       end if

    end do

  end subroutine full_mat

  !subroutine to convert full matrix to symmetric form
  pure subroutine sym_mat(c_m,A)

    implicit none

    type(contact_mat), intent(inout) :: c_m
    real(rp), intent(in) :: A(1:c_m%N_CG,1:c_m%N_CG)

    integer(int32) :: i_pair, b1, b2

    do i_pair=1,c_m%N_pair

       b1 = c_m%pairs(1,i_pair)
       b2 = c_m%pairs(2,i_pair)
       
       c_m%rval(i_pair) =  A(b1,b2)

    end do
    
  end subroutine sym_mat

  !subroutine to convert full matrix to symmetric form
  pure subroutine sym_mat_base(c_m,A)

    implicit none

    type(contact_mat), intent(inout) :: c_m
    real(rp), intent(in) :: A(1:c_m%N_base_CG,1:c_m%N_base_CG)

    integer(int32) :: i_pair, b1, b2

    do i_pair=1,c_m%N_pair_base

       b1 = c_m%pairs_base(1,i_pair)
       b2 = c_m%pairs_base(2,i_pair)

       c_m%rval_base(i_pair) =  A(b1,b2)

    end do

  end subroutine sym_mat_base

  !routine to perform a binary search for the indices in a sorted array
  !that bracket the target
  pure subroutine binary_search(lo,hi,A,N,tar)

    implicit none

    integer, intent(in) :: N, tar
    integer, intent(in) :: A(1:N)
    integer, intent(inout) :: lo, hi

    integer :: mid
    logical :: target_found, find

    if (N .gt. 1) then

       !select midpoint
       lo=1
       hi=N

       target_found=.false.

       do while((target_found .eqv. .false.) .and. (lo .lt. hi))

          mid=lo+(hi-lo)/2

          if (A(mid) .eq. tar) then
             target_found=.true.

             !target is above the midpoint
          elseif (A(mid) .lt. tar) then
             lo=mid+1
             !R=N

             !target is below the midpoint
          elseif (A(mid) .gt. tar) then
             !lo=1
             if (mid .ne. 1) then
                hi=mid-1
             else
                hi=1
             end if

          end if
       end do

       if ((lo .eq. hi) .and. (A(lo) .eq. tar)) then
          mid=lo
          target_found=.true.
       end if

       !target has been found, now search for full interval with target value
       if (target_found .eqv. .true.) then

          lo=mid
          hi=mid

          !search for the lower bound of the interval containing the target value
          find=.true.
          do while(find .eqv. .true.)
             if (lo .gt. 1) then
                if (A(lo-1) .eq. tar) then
                   lo=lo-1
                else
                   find=.false.
                end if
             else
                find=.false.
             end if
          end do

          !search for the upper bound of the interval containing the target value
          find=.true.
          do while(find .eqv. .true.)
             if (hi .lt. N) then
                if (A(hi+1) .eq. tar) then
                   hi=hi+1
                else
                   find=.false.
                end if
             else
                find=.false.
             end if
          end do

       else !if the target was not found, return -1's
          lo=-1
          hi=-1
       end if

    else

       !check if only element is the target
       if (A(1) .eq. tar) then
          lo=1
          hi=1
       else !if the target was not found, return -1's
          lo=-1
          hi=-1
       end if

    end if !end conditional checking if target was found

  end subroutine binary_search

end module contact_object
