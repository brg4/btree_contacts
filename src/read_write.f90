module read_write

  !module used to handle reading the input file and
  !writing the restart and output files
  
  use params
  use contact_object

  implicit none

contains

  !routine to load the input file
  subroutine load_input_file(a_p,input_file)

    implicit none

    type(analysis_params), intent(inout) :: a_p
    character(len=240), intent(in) :: input_file

    integer(int32) :: unit_input, stat
    integer(int32) :: run_flag, min_rep, max_rep, N_reps, map_flag
    integer(int32) :: method, N, CG, data_size, N_link_params
    integer(int32) :: i_link_param

    real(rp) :: r
    real(rp), allocatable :: link_params(:)
    
    character(len=240) :: temp_line, arg_specifier, arg
    character(len=240) :: data_label, data_dir, map_file

    unit_input=offset_units

    run_flag = 0
    min_rep = 0
    max_rep = 0
    
    open(UNIT=unit_input,FILE=trim(input_file),ACTION='read',STATUS='old')

        do

       read(unit_input,"(A)",IOSTAT=stat)temp_line
       !write(unit_log,*)temp_line

       if (IS_IOSTAT_END(stat)) then
          exit
       else

          !take argument from terms with "="
          if (index(temp_line,"=") .ne. 0) then
             arg_specifier=trim(adjustl(temp_line(:index(temp_line,"=")-1)))
             arg=trim(adjustl(temp_line(index(temp_line,"=")+1:)))
          else
             arg_specifier=trim(adjustl(temp_line))
          end if

          select case(arg_specifier)

          case ('run')

             read(arg,"(I5)")run_flag

          case ('min_rep')

             read(arg,"(I5)")min_rep

          case ('max_rep')

             read(arg,"(I5)")max_rep

          case ('method')

             read(arg,"(I5)")method
             
          case ('r')

             read(arg,*)r

          case ('N_link_params')

             read(arg,"(I5)")N_link_params

             if (allocated(link_params)) deallocate(link_params)
             allocate(link_params(1:N_link_params))

             !read the link parameter list
             do i_link_param=1,N_link_params
                
                read(unit_input,"(A)",IOSTAT=stat)temp_line

                if (index(temp_line,"=") .ne. 0) then
                   arg=trim(adjustl(temp_line(index(temp_line,"=")+1:)))
                   read(arg,*)link_params(i_link_param)
                end if

             end do

          case ('CG')

             read(arg,"(I5)")CG

          case ('N')

             read(arg,"(I6)")N

          case ('map_flag')

             read(arg,"(I5)")map_flag

          case ('data_label')

             read(arg,"(A)")data_label

          case ('data_dir')

             read(arg,"(A)")data_dir
             
          case ('map_file')

             read(arg,"(A)")map_file

          case ('data_size(bytes)')

             read(arg,"(I5)")data_size
             
          end select

       end if
       
    end do

    !close the input file
    close(unit_input)
    
    call new_analysis_params(a_p,&
       run_flag,&
       min_rep,max_rep,&
       method,&
       N_link_params,&
       r,link_params,&
       CG,N,map_flag,&
       data_label,data_dir,map_file,&
       data_size)
    
  end subroutine load_input_file

  subroutine new_contact_mat_from_map_file(a_p,c_m)

    implicit none

    type(analysis_params), intent(inout) :: a_p
    type(contact_mat), intent(inout) :: c_m

    integer(int32) :: unit_input
    integer(int32) :: i_CG, i_temp, i, j, c, N_ter

    character(len=240) :: temp_line

    if (allocated(c_m%pairs)) deallocate(c_m%pairs)
    if (allocated(c_m%map)) deallocate(c_m%map)
    if (allocated(c_m%base)) deallocate(c_m%base)
    if (allocated(c_m%N_mono)) deallocate(c_m%N_mono)
    if (allocated(c_m%ival)) deallocate(c_m%ival)
    if (allocated(c_m%rval)) deallocate(c_m%rval)
    if (allocated(c_m%ival_base)) deallocate(c_m%ival_base)
    if (allocated(c_m%rval_base)) deallocate(c_m%rval_base)

    unit_input = offset_units

    open(UNIT=unit_input,FILE=trim(a_p%map_file),ACTION='read',STATUS='old')

    !read the base number of points
    read(unit_input,"(A)")temp_line
    temp_line=trim(adjustl(temp_line(index(temp_line,"=")+1:)))
    read(temp_line,"(I7)")a_p%N_base
    
    !read the number of points
    read(unit_input,"(A)")temp_line
    temp_line=trim(adjustl(temp_line(index(temp_line,"=")+1:)))
    read(temp_line,"(I7)")a_p%N

    !read the coarse-graining factor
    read(unit_input,"(A)")temp_line
    temp_line=trim(adjustl(temp_line(index(temp_line,"=")+1:)))
    read(temp_line,"(I6)")a_p%CG

    !read the number of base loci
    read(unit_input,"(A)")temp_line
    temp_line=trim(adjustl(temp_line(index(temp_line,"=")+1:)))
    read(temp_line,"(I7)")c_m%N_base_CG

    !read the total number of loci
    read(unit_input,"(A)")temp_line
    temp_line=trim(adjustl(temp_line(index(temp_line,"=")+1:)))
    read(temp_line,"(I7)")c_m%N_CG

    !skip two lines
    read(unit_input,*)
    read(unit_input,*)

    !read the number of terminal chromosomes
    read(unit_input,"(A)")temp_line
    temp_line=trim(adjustl(temp_line(index(temp_line,"=")+1:)))
    read(temp_line,"(I7)")N_ter

    !skip lines equal to the number of terminal chromosomes + 2
    do c=1,(N_ter+2)
       read(unit_input,*)
    end do
    
    !create mapping
    allocate(c_m%N_mono(1:c_m%N_CG))
    allocate(c_m%base(1:c_m%N_CG))
    allocate(c_m%map(1:2,1:c_m%N_CG))

    !read the details of the CG loci
    do i_CG=1,c_m%N_CG
       read(unit_input,*)i_temp,c_m%base(i_CG),c_m%map(1,i_CG),c_m%map(2,i_CG)
       !write(unit_log,*)i_CG,c_m%base(i_CG),c_m%map(1,i_CG),c_m%map(2,i_CG)
    end do

    !close the map file
    close(unit_input)

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

    !create pairwise indices for base
    
    c_m%N_pair_base = c_m%N_base_CG*(c_m%N_base_CG-1)/2 + c_m%N_base_CG
    allocate(c_m%pairs_base(1:2,1:c_m%N_pair_base))

    c = 1
    do i=1,c_m%N_base_CG
       do j=i,c_m%N_base_CG
          c_m%pairs_base(:,c) = (/i,j/)
          c = c + 1
       end do
    end do

    !allocate value arrays (integer and real)

    allocate(c_m%ival(1:c_m%N_pair))
    allocate(c_m%rval(1:c_m%N_pair))
    allocate(c_m%ival_base(1:c_m%N_pair_base))
    allocate(c_m%rval_base(1:c_m%N_pair_base))

    
    c_m%ival = 0
    c_m%rval = 0.0d0

    c_m%ival_base = 0
    c_m%rval_base = 0.0d0

  end subroutine new_contact_mat_from_map_file

  subroutine read_coords_binary(x,unit,file,m,n,data_size)

    implicit none

    real(rp), intent(out) :: x(m,n)
    integer(int16), intent(in) :: unit
    integer(int32), intent(in) :: m, n, data_size
    character(len=*), intent(in) :: file

    integer(int32) :: i

    !open the file
    open(UNIT=unit,&
         FILE=trim(adjustl(file)),&
         STATUS='old',&
         FORM='unformatted',&
         RECL=m*data_size,&
         ACCESS='direct')

    do i=1,n
       read(unit,rec=i)x(:,i)
    end do

    !close the file
    close(unit)
    
  end subroutine read_coords_binary

  subroutine write_cm_pairs_binary(c_m,unit,file)

    implicit none

    type(contact_mat), intent(in) :: c_m
    integer(int32), intent(in) :: unit
    character(len=*), intent(in) :: file

    !open the file
    open(UNIT=unit,&
         FILE=trim(adjustl(file)),&
         STATUS='replace',&
         FORM='unformatted',&
         RECL=2*c_m%N_pair*int32,&
         ACCESS='direct')

    write(unit,rec=1)c_m%pairs

    !close the file
    close(unit)

  end subroutine write_cm_pairs_binary

  subroutine write_cm_rval_binary(c_m,unit,file)

    type(contact_mat), intent(in) :: c_m
    integer(int32), intent(in) :: unit
    character(len=*), intent(in) :: file

    !open the file
    open(UNIT=unit,&
         FILE=trim(adjustl(file)),&
         STATUS='replace',&
         FORM='unformatted',&
         RECL=c_m%N_pair*rp,&
         ACCESS='direct')

    write(unit,rec=1)c_m%rval

    !close the file
    close(unit)
    
  end subroutine write_cm_rval_binary

    subroutine write_cm_pairs_base_binary(c_m,unit,file)

    implicit none

    type(contact_mat), intent(in) :: c_m
    integer(int32), intent(in) :: unit
    character(len=*), intent(in) :: file

    !open the file
    open(UNIT=unit,&
         FILE=trim(adjustl(file)),&
         STATUS='replace',&
         FORM='unformatted',&
         RECL=2*c_m%N_pair_base*int32,&
         ACCESS='direct')

    write(unit,rec=1)c_m%pairs_base

    !close the file
    close(unit)

  end subroutine write_cm_pairs_base_binary

  subroutine write_cm_rval_base_binary(c_m,unit,file)

    type(contact_mat), intent(in) :: c_m
    integer(int32), intent(in) :: unit
    character(len=*), intent(in) :: file

    !open the file
    open(UNIT=unit,&
         FILE=trim(adjustl(file)),&
         STATUS='replace',&
         FORM='unformatted',&
         RECL=c_m%N_pair_base*rp,&
         ACCESS='direct')

    write(unit,rec=1)c_m%rval_base

    !close the file
    close(unit)
    
  end subroutine write_cm_rval_base_binary

end module read_write
