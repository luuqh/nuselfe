!===============================================================================
!===============================================================================
! Routine to read in param.in
!===============================================================================
!===============================================================================
subroutine get_param(varname,vartype,ivarvalue,varvalue1,varvalue2)
  ! Get a parameter from param.in
  ! Inputs:
  !        varname: parameter name (string no longer than 90)
  !        vartype: parameter value type (0: 2-char string; 1: integer; 2: float)
  ! Outputs:
  !        ivarvalue: integer output;
  !        varvalue1: float output;
  !        varvalue2: 2-char string output.
  ! Format rules for param.in:
  ! (1) Lines beginning with "!" are comments; blank lines are ignored;
  ! (2) one line for each parameter in the format: keywords= value;
  !     keywords are case sensitive; spaces allowed between keywords and "=" and value;
  !     comments starting with "!"  allowed after value;
  ! (3) value is an integer, double, or 2-char string; for double, any of the format is acceptable:
  !     40 40. 4.e1
  !     Use of decimal point in integers is OK but discouraged.
  use elfe_glbl, only : rkind,errmsg
  use elfe_msgp, only : parallel_abort,myrank
  implicit real(rkind)(a-h,o-z), integer(i-n)

  character(*),intent(in) :: varname
  integer,intent(in) :: vartype
  integer,intent(out) :: ivarvalue
  real(rkind),intent(out) :: varvalue1
  character(len=2),intent(out) :: varvalue2

  character(len=90) :: line_str,str_tmp,str_tmp2

  str_tmp2=adjustl(varname)
  lstr_tmp2=len_trim(str_tmp2)
!  print*, varname !,str_tmp2(1:lstr_tmp2)

  ! Scan param.in
  open(15,file='param.in',status='old')
  rewind(15)
  line=0
  do
    line=line+1
    read(15,'(a)',end=99)line_str
    line_str=adjustl(line_str) !place blanks at end
    len_str=len_trim(line_str)
    if(len_str==0.or.line_str(1:1)=='!') cycle

    loc=index(line_str,'=')
    loc2=index(line_str,'!')
    if(loc2/=0.and.loc2-1<loc+1) call parallel_abort('READ_PARAM: ! before =')

    str_tmp=''
    str_tmp(1:loc-1)=line_str(1:loc-1) !keyword
    str_tmp=trim(str_tmp)
    lstr_tmp=len_trim(str_tmp)
    
    if(str_tmp(1:lstr_tmp)==str_tmp2(1:lstr_tmp2)) then
       if(loc2/=0) then
         str_tmp2=line_str(loc+1:loc2-1)
       else
         str_tmp2=line_str(loc+1:len_str)
       endif
       str_tmp2=adjustl(str_tmp2)
       str_tmp2=trim(str_tmp2)
       if(vartype==0) then !string
         varvalue2=str_tmp2(1:2)
#ifdef DEBUG
         if(myrank==0) write(86,*)varname,' = ',varvalue2
#endif
       else if(vartype==1) then !integer
         read(str_tmp2,*)ivarvalue
#ifdef DEBUG
         if(myrank==0) write(86,*)varname,' = ',ivarvalue
#endif
       else if(vartype==2) then !float
         read(str_tmp2,*)varvalue1
#ifdef DEBUG
         if(myrank==0) write(86,*)varname,' = ',real(varvalue1)
#endif
       else
         write(errmsg,*)'read_param: unknown type:',vartype
         call parallel_abort(errmsg)
       endif
       exit
    endif
  enddo !scan param.in
  
!   print*, 'Found it on line: ',line
   close(15)
   return

99  close(15)
   write(errmsg,*)'Failed to find parameter:',varname
   call parallel_abort(errmsg)

end subroutine get_param

