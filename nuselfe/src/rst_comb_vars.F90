MODULE rst_comb_vars_mod

  implicit real(4)(a-h,o-z),integer(i-n)
  parameter(nbyte=1)
  parameter(mnp=80000)
  parameter(mne=160000)
  parameter(mnv=200)  
  character*30 file63
  character*12 it_char
  character*48 start_time,version,variable_nm,variable_dim
  character*48 data_format
  character(72)     :: fgb,fgb2,fdb  ! Processor specific global output file name
  integer           :: lfgb,lfdb       ! Length of processor specific global output file name
  character(len= 4) :: a_4
  allocatable ne(:),np(:),nsproc(:),ihot_len(:)
  allocatable ztot(:),sigma(:),outb(:,:,:),eta2(:)
  allocatable i34(:),nm(:,:),nm2(:,:),js(:,:),xctr(:),yctr(:),dpe(:)
  allocatable x(:),y(:),dp(:),kbp(:),iplg(:,:),ielg(:,:),kbp01(:,:)
  allocatable xcj(:),ycj(:),dps(:)
  dimension   kfp(mnp),hmod(mnp),cs(mnv),z(0:mnv,mnp),htot(mnp),delz(mnv,mnp)
  
END MODULE 