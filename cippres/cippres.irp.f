

program cippres
  use general
  implicit none
  BEGIN_DOC
! CIPPRES stands for Configuration Interaction Plugin for Photoionized and Resonant Electronic States.
! ORMAS-CI calculations are performed to compute bound and approximate continuum states. Hamiltonian and Dipole coupling matrix elements can then be computed between these states. Finally, Stieltjes imaging technique is applied to recover the correct continuum coupling matrix elements.
  END_DOC

  character(len=lenmax) :: finput
  integer :: ilen, jlen
  logical :: file_e

  integer :: i, j


  if(ifcsf==1) then

    do i = 1, n_ciruns_cippres
!      print*,'ncsfs',i, n_csf_cippres(i) 
      print*,'CI run = ', i
      print*,'CI eigval =', eigvalues_cippres(1:n_csf_cippres(i),i)
    enddo

  else

! READ XML file containing the ORMAS info
   print*, 'Reads ', finput_cippres
   print*, 'And generate the list of CSFs '

   finput=finput_cippres
   jlen=index(finput,' ')
   inquire( file="./"//finput(1:jlen-1), exist=file_e )
   if ( file_e .eqv. .false. ) then
    write(*,*) finput(1:jlen-1), " does not exist"
    stop
   endif
   call generate_csfs(finput)
   call ezfio_set_cippres_ifcsf(1)

  endif

end program cippres
 
subroutine generate_csfs(finput)
  use general
  implicit none
  integer                        :: i,j
  character(len=lenmax), intent(in) :: finput

! the python script will write the info and list of CSFs, for as many CIruns as given in finput, into header$irun.txt and list$irun.txt ($irun=1,2,3,....)
  call system('python $QP_ROOT/plugins/local/cippres/scripts/generate_csfs.py '//trim(finput))

! save the info and list of CSFs into EZFIO
  open(unit=10,file='parser.txt') 
   read(10,*) n_ciruns_cippres
   call ezfio_set_cippres_n_ciruns_cippres(n_ciruns_cippres)
   read(10,*) n_csf_max
   call ezfio_set_cippres_n_csf_max(n_csf_max)
   read(10,*) n_det_max_csf
   call ezfio_set_cippres_n_det_max_csf(n_det_max_csf)
  close(10)

end subroutine generate_csfs

