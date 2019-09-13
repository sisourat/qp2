program tests
  use general
  implicit none
  BEGIN_DOC
! reads the wf and computes the dipolar moment of that wf
  END_DOC
  read_wf = .True.
  touch read_wf 
 !call routine_E
 !call play_with_wf
  print*,'lenmax = ',lenmax
end

subroutine routine_E
 implicit none
 BEGIN_DOC
! computes the Energy of the WF in the EZFIO
 END_DOC
 use bitmasks ! you need to include the bitmasks_module.f90 features
 integer :: i,j,k,istate
 integer(bit_kind), allocatable :: det_i(:,:),det_j(:,:)
 allocate(det_i(N_int,2))
 allocate(det_j(N_int,2))

 double precision, allocatable :: accu(:)
 allocate(accu(N_states))
 double precision :: hij
 accu = 0.d0
 do istate = 1, N_states
  do i = 1, N_det
   ! first copy 
   do k = 1, N_int
   ! alpha = 1 / beta = 2
    det_i(k,1) = psi_det(k,1,i)
    det_i(k,2) = psi_det(k,2,i)
   enddo
   do j = 1, N_det
    ! second copy 
    do k = 1, N_int
     det_j(k,1) = psi_det(k,1,j)
     det_j(k,2) = psi_det(k,2,j)
    enddo
    ! you get the H matrix element 
    call i_H_j(det_i,det_j,N_int,hij)
!   call i_H_j(psi_det(1,1,i),psi_det(1,1,j),N_int,hij) ! equivalent but without copy
    ! you get the coef 
    accu(istate) += hij * psi_coef(j,istate) * psi_coef(i,istate)
   enddo
  enddo
 enddo
 print*,'accu               = ',accu(:) 
 print*,'ref_bitmask_energy = ',ref_bitmask_energy

end

subroutine  routine_D
 implicit none
 BEGIN_DOC
! computes the Dipole moment of the WF in the EZFIO
 END_DOC
 use bitmasks ! you need to include the bitmasks_module.f90 features
 integer :: i,j,k,istate,degree
 integer(bit_kind), allocatable :: det_i(:,:),det_j(:,:)
 double precision, allocatable :: accu(:)
 double precision :: dij
 integer, allocatable           :: occ(:,:)
 integer                        :: n_occ_ab(2)
 integer :: iorb,jorb
 integer :: exc(0:2,2,2)
 double precision :: phase
 integer :: h1,p1,h2,p2,s1,s2

 allocate(occ(N_int*bit_kind_size,2))

 allocate(det_i(N_int,2))
 allocate(det_j(N_int,2))
 allocate(accu(N_states))

 accu = 0.d0
 do istate = 1, N_states
  do i = 1, N_det
   ! first copy 
   do k = 1, N_int
   ! alpha = 1 / beta = 2
    det_i(k,1) = psi_det(k,1,i)
    det_i(k,2) = psi_det(k,2,i)
   enddo
   do j = 1, N_det
    ! second copy 
    do k = 1, N_int
     det_j(k,1) = psi_det(k,1,j)
     det_j(k,2) = psi_det(k,2,j)
    enddo
    call get_excitation_degree(det_i,det_j,degree,N_int)
    ! Filter the determinants which are singly excited wr to one another
    if(degree .gt.1)cycle

    if(i==j)then ! diagonal part
    call bitstring_to_list_ab(det_i, occ, n_occ_ab, N_int)
   !print*,'alpha electrons orbital occupancy'
   !do i = 1, n_occ_ab(1) ! browsing the alpha electrons
   ! print*,occ(i,1)
   !enddo
   !print*,'beta  electrons orbital occupancy'
   !do i = 1, n_occ_ab(2) ! browsing the beta  electrons
   ! print*,occ(i,2)
   !enddo
    dij = 0.d0
    do k = 1, n_occ_ab(1) ! browsing the alpha electrons
     iorb = occ(k,1) ! orbital of the ith alpha electron
     dij += mo_dipole_x(iorb,iorb) 
    enddo
    do k = 1, n_occ_ab(2) ! browsing the beta electrons
     iorb = occ(k,2) ! orbital of the ith beta electron
     dij += mo_dipole_x(iorb,iorb) 
    enddo
    else ! extra diagonal part
     call get_single_excitation(det_i,det_j,exc,phase,N_int) 
     call decode_exc(exc,degree,h1,p1,h2,p2,s1,s2)
     iorb = h1
     jorb = p1
     dij = phase * mo_dipole_x(iorb,jorb)
     ! you get the excitation operator btw det_i and det_j
     ! you get the coef 
    endif
    accu(istate) += dij * psi_coef(j,istate) * psi_coef(i,istate)
   enddo
  enddo
 enddo
 print*,'< psi | D | psi >  = ',accu(:)

end

subroutine play_with_wf
 implicit none
 use bitmasks ! you need to include the bitmasks_module.f90 features
 integer :: i,j,k
 integer :: n_det_tmp,n_a,n_b
 integer, allocatable :: occ_a(:,:)
 integer, allocatable :: occ_b(:,:)
 integer(bit_kind), allocatable :: psi_tmp(:,:,:)
 double precision, allocatable :: psi_coef_tmp(:,:)

 print*,''
 print*,''
 print*,'PAssing through the global variables as defined in the EZFIO data'
 print*,''
 print *,  'N_det = ', N_det
 print*,'******************************'
 print *,  'Energies  of the states:'
 do i = 1,N_states
   print *,  i, CI_energy(i)
 enddo
 print*,''
 print*,''
 print*,''

 n_det_tmp = 2
 n_a = elec_alpha_num
 n_b = elec_beta_num

 allocate( occ_a(mo_num,n_det_tmp), occ_b(mo_num,n_det_tmp)) 
 allocate(psi_tmp(N_int,2,n_det_tmp))
 allocate(psi_coef_tmp(n_det_tmp,n_det_tmp))

 integer :: idet
 ! first det is the first orbital 
 idet = 1
 occ_a( 1 , idet) = 1
 occ_b( 1 , idet) = 1
 call create_det(n_a,n_b,occ_a(1,idet),occ_b(1,idet),psi_tmp(1,1,idet))
 ! first det is some excited orbital 4 
 idet = 2
 occ_a( 1 , idet) = 2
 occ_b( 1 , idet) = 2
 call create_det(n_a,n_b,occ_a(1,idet),occ_b(1,idet),psi_tmp(1,1,idet))

 double precision, allocatable :: hmat(:,:),eigval(:),eigvec(:,:)
 double precision :: hij
 allocate(hmat(n_det_tmp,n_det_tmp),eigval(n_det_tmp),eigvec(n_det_tmp,n_det_tmp))
 print*,'Computing the H mat and printing the determinants '
 do i = 1, n_det_tmp
  print*,'i = ',i
  call debug_det(psi_tmp(1,1,i),N_int)
  do j = 1, n_det_tmp
   call i_H_j(psi_tmp(1,1,i),psi_tmp(1,1,j),N_int,hij)
   hmat(j,i) = hij
  enddo
 enddo
 print*,''
 print*,'H matrix '
 print*,''
 do i = 1, n_det_tmp
  write(*,'(100(F10.5,X))')hmat(i,:)
 enddo
 call lapack_diagd(eigval,eigvec,hmat,n_det_tmp,n_det_tmp) 
 print*,''
 print*,''
 do i = 1, n_det_tmp
  print*,'eigval = ',i,eigval(i)
 enddo

 ! Reinitialize the WF to zero 
 psi_det = 0_bit_kind 
 psi_coef = 0.d0 
 ! you copy the new set of determinants into psi_det
 do i = 1, n_det_tmp
  do j = 1, N_int 
   psi_det(j,1,i) = psi_tmp(j,1,i)
   psi_det(j,2,i) = psi_tmp(j,2,i)
  enddo
 enddo
 N_det = n_det_tmp
 psi_coef(1,1) = 1.d0
 touch N_det psi_coef psi_det
 
 print*,''
 print*,''
 print*,'PAssing through the global variables ...'
 print*,''
 print *,  'N_det = ', N_det
 print*,'******************************'
 print *,  'Energies  of the states:'
 do i = 1,N_states
   print *,  i, CI_energy(i)
 enddo

 
end
