
program GB_Bmatrix_1505
  implicit none
!--- for subroutine bmat
  integer :: nf, nbl, nba, nto, nob, nod, nlbe, natom, nint
  integer,allocatable :: IBL(:,:),IBA(:,:),ITO(:,:),IOB(:,:),IOD(:,:),ILBE(:,:)
  double precision :: t(3,3)
  double precision,allocatable :: BM(:,:), X(:),Y(:),Z(:)
!--- for building A matrix
  integer :: lwork,info
  integer, allocatable :: indx(:)
  double precision, allocatable :: amat(:,:),bubt(:,:),u(:,:),work(:), eigval(:) &
       , eigvalmat(:,:)
!--- my own
!  double precision,parameter :: PI=3.14159265358979d0
  integer :: itmp, i, j, ifcart=101, ifnacmec=102, ifibl=103, ifiba=104, ifito=105 &
       , ifiob=106, ifiod=107, ifnrdef=108 &
       , ifbmat=201, ifamat=202, ifbmatr=203, itmp2, nintnr
  integer,allocatable :: iarraytmp(:)
  double precision :: rtmp
  double precision,allocatable :: mass(:), bmnr(:,:), tmat(:,:) &
       , rarraytmp(:)
  character(len=2) :: atomsym
  logical :: booltmp
!--- open files
  open(ifcart, file='cart.txt', action='read')
  open(ifibl, file='ibl.txt', action='read')
  open(ifiba, file='iba.txt', action='read')
  open(ifito, file='ito.txt', action='read')
  open(ifiob, file='iob.txt', action='read')
  open(ifiod, file='iod.txt', action='read')
  open(ifbmat, file='bmat.txt', action='write')
  open(ifamat, file='amat.txt', action='write')
!--- (the nacme part has been moved to nacme_cart2int.f90)
!  open(ifnacmec, file='nacme_cart.txt', action='read')
!--- read in number of items
  read(ifcart, *) natom
  read(ifibl, *) nbl
  read(ifiba, *) nba
  read(ifito, *) nto
  read(ifiob, *) nob
  read(ifiod, *) nod
  nlbe = 0; nint = nbl+nba+nto+nob+nod+nlbe; nf = nint
!--- allocate memory
  allocate(ibl(2,nf),iba(3,nf),ito(4,nf),iob(4,nf),iod(4,nf) &
       ,ilbe(3,nf),bm(nint,natom*3))
  allocate(x(natom),y(natom),z(natom),mass(natom*3),u(natom*3,natom*3))
!--- (the nacme part has been moved to nacme_cart2int.f90)
  ! allocate(nacmec(natom*3),nacmei(nint))
!--- read in cartesian coord
  do i=1, natom
     read(ifcart, *) atomsym, rtmp, x(i), y(i), z(i)
     mass(i*3-2) = rtmp; mass(i*3-1) = rtmp; mass(i*3) = rtmp
  end do
!--- read in definition of internal coord
  do i=1, nbl
     read(ifibl, *) ibl(1,i), ibl(2,i)
  end do
  do i=1, nba
     read(ifiba, *) iba(1,i), iba(2,i), iba(3,i)
  end do
  do i=1, nto
     read(ifito, *) ito(1,i), ito(2,i), ito(3,i), ito(4,i)
  end do
  do i=1, nob
     read(ifiob, *) iob(1,i), iob(2,i), iob(3,i), iob(4,i)
  end do
  do i=1, nod
     read(ifiod, *) iod(1,i), iod(2,i), iod(3,i), iod(4,i)
  end do
!--- read in cartesian nacme
!--- (the nacme part has been moved to nacme_cart2int.f90)
  ! do i=1, natom
  !    read(ifnacmec, *) itmp, nacmec(i*3-2), nacmec(i*3-1), nacmec(i*3)
  ! end do
!--- calculate B matrix. Subroutine bmat from frompolyrate.f90 of mstor. 
!--- currently only ibl(bond len), iba(bond angle), ito(torsion), iob(oop bend) are used
!--- ilbe and t are dummy
  call bmat(NF,NBL,NBA,NTO,NOB,NOD,NLBE,IBL,IBA,ITO,IOB,IOD,ILBE &
       ,X,Y,Z,NATOM,NINT,BM,T)
!--- write B matrix
  open(ifbmatr, file='bmatr.txt', action='write')
  write(ifbmatr,*) nint, natom*3
  do i=1, nint
     do j=1, natom*3
        write(ifbmatr,'(f15.6)',advance='no') bm(i,j)
     end do
     write(ifbmatr,*)
  end do
!--- if intnrdef.txt exists, transform B matrix to nonredundant coordinates
  inquire(file='intnrdef.txt', exist=booltmp)
  if(booltmp) then
     open(ifnrdef, file='intnrdef.txt', action='read')
     read(ifnrdef, *) nintnr
     ! if(nintnr.ne.natom*3-6) then
     !    write(*, *) 'Number of nonredundant int. coord. != natom*3-6. Stop.'
     !    stop
     ! end if
     allocate(bmnr(nintnr,natom*3),tmat(nintnr,nint))
     tmat = 0d0
     do i=1, nintnr
        read(ifnrdef, *) itmp, itmp2
        allocate(iarraytmp(itmp2),rarraytmp(itmp2))
        read(ifnrdef, *) iarraytmp(:)
        read(ifnrdef, *) rarraytmp(:)
        do j=1, itmp2
           tmat(itmp,iarraytmp(j)) = rarraytmp(j)
        end do
        deallocate(iarraytmp,rarraytmp)
     end do
     bmnr = matmul(tmat, bm)
     deallocate(bm)
     nint = nintnr
     allocate(bm(nintnr,natom*3))
     bm = bmnr
     deallocate(bmnr)
  end if
!  write(*,*) 'debug1 before write bmat'
!--- write B matrix
  write(ifbmat,*) nint, natom*3
  do i=1, nint
     do j=1, natom*3
        write(ifbmat,'(f15.6)',advance='no') bm(i,j)
     end do
     write(ifbmat,*)
  end do
!=== form A = u*B^T*(B*u*B^T)^-1. Code partly adapted from vib.f90 of mstor. 
!--- calc B*u*B^T
  u(:,:) = 0.0d0
  do i = 1, natom
     u(3*i-2,3*i-2) = 1.0d0/mass(i)
     u(3*i-1,3*i-1) = u(3*i-2,3*i-2)
     u(3*i  ,3*i  ) = u(3*i-2,3*i-2)
  enddo
  allocate(indx(nint),amat(natom*3,nint),bubt(nint,nint))
  bubt(:,:)=matmul(bm(:,:),matmul(u(:,:),transpose(bm(:,:))))
!--- Scheme 2:
!--- calculate eigenvalues and eigenvectors of bubt in order to invert it
!--- This is slow but needed for the "generalized inverse" 
!--- if the internals are redundant
  !000 intel mkl routine dsyev(jobz, uplo, n, a, lda, w, work, lwork, info)
  ! allocate(eigval(nint),eigvalmat(nint,nint))
  ! lwork = -1
  ! allocate(work(1))
  ! call dsyev('V', 'U', nint, bubt, nint, eigval, work, lwork, info)
  ! lwork=int(work(1))
  ! deallocate (work)
  ! allocate( work(lwork) )
  ! call dsyev('V', 'U', nint, bubt, nint, eigval, work, lwork, info)
  ! if (info.ne.0 ) then
  !    write(6,*) "dsyev exit abnormally"
  !    stop
  ! endif
  ! deallocate (work)
!--- invert nonsingular eigenvalues
!--- Such "generalized inverse" can deal with redundant internals
!--- See Pulay and Fogarasi, JCP 1992, 96, 2856.
  ! eigvalmat = 0d0
  ! do i = 1, nint
  !    if (abs(eigval(i))>1d-7) then
  !       eigvalmat(i,i) = 1d0/eigval(i)
  !    else
  !       eigvalmat(i,i) = eigval(i)
  !    end if
  ! end do
  
!--- calculate generalized inverse of B*u*B^T
  ! bubt = matmul(bubt, matmul(eigvalmat, transpose(bubt)))
!--- Scheme 2 end

!--- Scheme 1: direct inverse of B*u*B^T
!--- Seems not able to handle redundant internals
!--- following from vib.f90 of mstor
  call dgetrf(nint,nint,bubt,nint,indx,info)
  if (info.ne.0 ) then
     write(6,*) "dgetrf exit abnormally"
     stop
  endif
  lwork = -1
  allocate(work(1))
  call dgetri(nint,bubt,nint,indx,work,lwork,info)
  lwork=int(work(1))
  deallocate (work)
  allocate( work(lwork) )
  call dgetri(nint,bubt,nint,indx,work,lwork,info)
  deallocate(work)
  if (info.ne.0 ) then
     write(6,*) "dgetri exit abnormally"
     stop
  endif
!--- Scheme 1 end
!--- calculate and write A matrix
  amat(:,:)=matmul( (matmul(u(:,:),transpose(bm(:,:)))),bubt(:,:))
  write(ifamat,*) natom*3, nint
  do i=1, natom*3
     do j=1, nint
        write(ifamat,'(f30.10)',advance='no') amat(i,j)
!        write(ifamat,*) amat(i,:)
     end do
     write(ifamat,*)
  end do
  ! write(*,*)
  ! do i=1, natom*3
  !    write(*,'(f6.2)',advance='no') amat(i,12)
  !    if(mod(i, 3).eq.0) write(*,*)
  ! end do
  ! write(*,*)
  ! do i=1, natom*3
  !    write(*,'(f6.2)',advance='no') amat(i,39)
  !    if(mod(i, 3).eq.0) write(*,*)
  ! end do
!--- calculate internal nacme
!--- (the nacme part has been moved to nacme_cart2int.f90)
  ! nacmei(:) = matmul(transpose(amat(:,:)), nacmec(:))
  ! do i=1, nint
  !    write(*,*) nacmei(i)
  ! end do

end program GB_Bmatrix_1505
