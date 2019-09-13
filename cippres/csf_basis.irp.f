use bitmasks ! you need to include the bitmasks_module.f90 features
use general

 BEGIN_PROVIDER [integer(bit_kind), csf_basis, (N_int,2, n_det_max_csf, n_csf_max ,n_ciruns_cippres )]
&BEGIN_PROVIDER [double precision, coef_det_csf_basis, (n_det_max_csf, n_csf_max ,n_ciruns_cippres )]
&BEGIN_PROVIDER [integer, n_csf_cippres, (n_ciruns_cippres )]
&BEGIN_PROVIDER [integer, n_sta_cippres, (n_ciruns_cippres )]
&BEGIN_PROVIDER [double precision, prttol_cippres, (n_ciruns_cippres )]
&BEGIN_PROVIDER [integer, n_det_csf_cippres, (n_csf_max,n_ciruns_cippres )]
 implicit none
  integer :: irun
  integer :: i,j,k
  integer :: nalpha, nbeta
  character(len=lenmax)             :: fileigvec

 integer, allocatable :: occ_a(:,:)
 integer, allocatable :: occ_b(:,:)
 character(1) :: ca, cb

 allocate(occ_a(mo_num,n_det_max_csf), occ_b(mo_num,n_det_max_csf)) 

! read the info for each CI run
  do irun = 1, n_ciruns_cippres

   open(unit=21,file='header'//achar(48+irun)//'.txt')
   open(unit=22,file='list'//achar(48+irun)//'.txt')

    read(21,*)n_sta_cippres(irun),prttol_cippres(irun)
    read(21,*)nalpha,nbeta
    if(nalpha/=elec_alpha_num) then
      print*,"Input inconsistent, nalpha in cirun",irun
    endif
    if(nbeta/=elec_beta_num) then
      print*,"Input inconsistent, nbeta in cirun",irun
    endif
    read(21,*)n_csf_cippres(irun)

   do i = 1, n_csf_cippres(irun)
    read(22,*)n_det_csf_cippres(i,irun) 
    do j = 1, n_det_csf_cippres(i,irun)
     read(22,*)coef_det_csf_basis(j,i,irun),ca,(occ_a(k,j),k=1,nalpha),cb,(occ_b(k,j),k=1,nbeta)
     call create_det(nalpha,nbeta,occ_a(1,j),occ_b(1,j),csf_basis(1,1,j,i,irun))
    enddo
   enddo 

   close(21)
   close(22)
  enddo

!! see Manu  call ezfio_set_cippres_csf_basis(csf_basis)

  deallocate(occ_a,occ_b)
END_PROVIDER 


BEGIN_PROVIDER [double precision, H_matrix_cippres, (n_csf_max,n_csf_max,n_ciruns_cippres)]
 implicit none
 integer :: irun
 integer :: i,j,k,l
 double precision :: hij
 H_matrix_cippres = 0.d0
 do irun = 1, n_ciruns_cippres
  do i = 1, n_csf_cippres(irun) ! first loop on the csf of the space ispace 
   do j = 1, n_csf_cippres(irun)
    do k = 1, n_det_csf_cippres(i,irun) ! then on the determinants belonging to the ith CSF of space ispace
     do l = 1, n_det_csf_cippres(j,irun)
      call i_H_j( csf_basis(1,1,k,i,irun) , csf_basis(1,1,l,j,irun),N_int,hij)
      H_matrix_cippres(j,i,irun) += hij * coef_det_csf_basis(k,i,irun) * coef_det_csf_basis(l,j,irun) 
     enddo
    enddo
   enddo
  enddo
 enddo
END_PROVIDER 


 BEGIN_PROVIDER [double precision, eigvectors_cippres, (n_csf_max,n_csf_max,n_ciruns_cippres)]
&BEGIN_PROVIDER [double precision, eigvalues_cippres, (n_csf_max,n_ciruns_cippres)]
 implicit none
 integer :: irun, j
 double precision, allocatable :: eigval(:),eigvec(:,:),hmat(:,:)

 do irun = 1, n_ciruns_cippres

  allocate(eigval(n_csf_cippres(irun)),eigvec(n_csf_cippres(irun),n_csf_cippres(irun)),hmat(n_csf_cippres(irun),n_csf_cippres(irun)))

!  print*,'nico',irun,H_matrix_cippres(1:n_csf_cippres(irun),1:n_csf_cippres(irun),irun)
  hmat(:,:) = H_matrix_cippres(1:n_csf_cippres(irun),1:n_csf_cippres(irun),irun)

  call lapack_diagd(eigval,eigvec,hmat,n_csf_cippres(irun),n_csf_cippres(irun)) 
  eigvalues_cippres(:,irun) = eigval(:)

!  do j = 1, n_csf_cippres(irun)
!   print*,'eigval',irun,j,eigval(j)
!  enddo

  deallocate(eigval,eigvec,hmat)

!  do j = 1, n_csf_cippres(irun)
!   eigvectors_cippres(j,irun) = eigvec(j,irun_que_tas_choisi)
!  enddo
 enddo
! call ezfio_set_cippres_eigvectors_cippres(eigvectors_cippres)
END_PROVIDER 
