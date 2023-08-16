program contact_calc

  use params
  use contact_object
  use contact_procedures
  use read_write
  use omp_lib

  implicit none

  logical :: log_specified, input_specified
  logical :: output_dir_specified, output_label_specified
  logical :: chain_flag

  integer(int32) :: num_threads
  integer(int16) :: i_arg, num_cl_args
  integer(int16) :: unit_rep
  integer(int32) :: i_rep
  integer(int32) :: i

  real(rp) :: res
  real(rp), allocatable :: x(:,:)
  
  character(len=240) :: input_file, log_file, output_dir
  character(len=240) :: output_label, output_label_rep, rep_id, data_label_rep
  character(len=240) :: arg_specifier, arg

  type(analysis_params) :: a_p
  type(contact_mat) :: c_m

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!! MAIN PROGRAM !!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  input_specified = .false.
  log_specified = .false.
  chain_flag = .false.
  output_dir_specified = .false.
  output_label_specified = .false.  

  num_cl_args=command_argument_count()

  if (num_cl_args .gt. 0) then
     
     do i_arg=1,num_cl_args
        
        call getarg(i_arg,arg)

        !take argument from terms with "="                                                                                      
        if (index(arg,"=") .ne. 0) then                                                                                         
           arg_specifier=trim(adjustl(arg(:index(arg,"=")-1)))                                                                  
           arg=trim(adjustl(arg(index(arg,"=")+1:)))                                                                            
        else                                                                                                                    
           arg_specifier=trim(adjustl(arg))                                                                                     
        end if

        select case (arg_specifier)
           !test if input is specified
        case ("--input_file", "--i_f")

           input_file=trim(arg)
           write(6,*)"input_file = ",trim(input_file)
           input_specified=.true.

           !test if output_dir is specified
        case ("--output_dir", "--o_d")

           output_dir=trim(arg)
           write(6,*)"output_dir = ",trim(output_dir)
           output_dir_specified=.true.

           !test if output_label is specified
        case ("--output_label", "--o_l")

           output_label=trim(arg)
           write(6,*)"output_label = ",trim(output_label)
           output_label_specified=.true.

           !test if log file was specified
        case ("--log", "--l")

           log_file=trim(arg)
           write(6,*)"log_file = ",trim(log_file)
           log_specified=.true.

           !test if chain_flag was present specified
        case ("--chain_flag", "--c_f")

           chain_flag=.true.
           write(6,*)"chain_flag = ",chain_flag

           !test if the number of threads were specified
        case ("--num_threads", "--n_t")

           read(arg,"(I3)")num_threads
           num_threads=min(num_threads,max_num_threads)
           
        end select
        
     end do
     
  end if

  !if the log unit is not stdout, open a log file at that unit number
  if (unit_log .ne. 6) then
     if (log_specified .eqv. .false.) then
        log_file="./run.log"
     end if
     open(UNIT=unit_log,FILE=trim(log_file),ACTION='write',STATUS='replace')
  end if

  write(unit_log,*)'[contact_calc] BEGIN'

  if (output_dir_specified .eqv. .false.) output_dir = "../output/"

  if (output_label_specified) then
     output_label = "_"//trim(output_label)
  else
     output_label = ""
  end if
  

  write(unit_log,*)' [contact_calc] read analysis parameters from input file'

  !read input file
  call load_input_file(a_p,input_file)

  !test for valid method parameters
  if (.not. method_required_link_params(a_p)) stop

  !create a new contact matrix depending on the method
  if (a_p%map_flag .eq. 1) then
     write(unit_log,*)' [contact_calc] read mapping from file'
     call new_contact_mat_from_map_file(a_p,c_m)
  else
     write(unit_log,*)' [contact_calc] generate default mapping'
     call new_contact_mat(c_m,a_p)
  end if

  write(unit_log,*)' [contact_calc] SYSTEM PARAMETERS'
  write(unit_log,*)'   run_flag = ',a_p%run_flag
  write(unit_log,*)'   min_rep = ',a_p%min_rep
  write(unit_log,*)'   max_rep = ',a_p%max_rep
  write(unit_log,*)'   N_reps = ',a_p%N_reps
  write(unit_log,*)'   method = ',a_p%method
  write(unit_log,*)'   r = ',a_p%r
  write(unit_log,*)'   N_link_params = ',a_p%N_link_params
  write(unit_log,*)'   CG = ',a_p%CG
  write(unit_log,*)'   N = ',a_p%N
  write(unit_log,*)'   N_CG = ',c_m%N_CG
  write(unit_log,*)'   data_label = ',trim(a_p%data_label)
  write(unit_log,*)'   data_dir = ',trim(a_p%data_dir)
  write(unit_log,*)'   data_size = ',a_p%data_size
  
  if ((a_p%N_reps .gt. 0) .and. (a_p%run_flag .ne. 0)) then

     if (allocated(x)) deallocate(x)
     allocate(x(1:3,1:a_p%N))
     
     do i_rep = 1,a_p%N_reps

        write(rep_id,"(A,I5.5,A)")'(',(a_p%min_rep+i_rep-1),')'
        rep_id=trim(rep_id)

        write(unit_log,*)' [contact_calc] '//trim(rep_id)//' BEGIN'

        unit_rep = offset_units

        !read the datafile

        data_label_rep = trim(a_p%data_label)//"_rep"
        write(data_label_rep,"(A,I5.5)")trim(data_label_rep),(a_p%min_rep+i_rep-1)

        if (chain_flag .eqv. .true.) then
           call read_coords_binary(x,unit_rep,&
                trim(a_p%data_dir)//'x_chain_'//trim(data_label_rep)//'.bin',&
                3,a_p%N,a_p%data_size)
        else
           call read_coords_binary(x,unit_rep,&
                trim(a_p%data_dir)//'x_'//trim(data_label_rep)//'.bin',&
                3,a_p%N,a_p%data_size)
        end if

        !calculate contributions to contact matrix

        call contact_contribution(c_m,a_p,x)

        !write the current contact matrix

        output_label_rep = trim(output_label)//"_rep"
        write(output_label_rep,"(A,I5.5)")trim(output_label_rep),(a_p%min_rep+i_rep-1)

        write(unit_log,*)' [contact_calc] '//trim(rep_id)//' END'

     end do

     if ((a_p%method .eq. 3) .or. (a_p%method .eq. 4)) then
        c_m%rval = 1.0d0*c_m%ival
        if (a_p%map_flag .eq. 1) c_m%rval_base = 1.0d0*c_m%ival_base
     else
        c_m%rval = c_m%rval/a_p%N_reps
        if (a_p%map_flag .eq. 1) c_m%rval_base = c_m%rval_base/a_p%N_reps
     end if

     write(unit_log,*)c_m%rval(1)

     call write_cm_pairs_binary(c_m,offset_units,&
          trim(output_dir)//'pairs'//trim(output_label)//'.bin')

     call write_cm_rval_binary(c_m,offset_units,&
          trim(output_dir)//'rval'//trim(output_label)//'.bin')

     if (a_p%map_flag .eq. 1) then

        call map_to_base(c_m)
        
        call write_cm_pairs_base_binary(c_m,offset_units,&
             trim(output_dir)//'pairs_mapped'//trim(output_label)//'.bin')

        call write_cm_rval_base_binary(c_m,offset_units,&
             trim(output_dir)//'rval_mapped'//trim(output_label)//'.bin')
        
     end if

  end if
  
  write(unit_log,*)'[contact_calc] END'

  !if not writing to stdout, then close the log file
  if (unit_log .ne. 6) then
     close(unit_log)
  end if

end program contact_calc
