subroutine vib(coord,masses,totmass,fc,fcm,freq,natoms,n3tm,nf,ibl,iba,ito,ilbe,nbl,nba,nto,nlbe,ntor, &
               ratio,nnegfreq,dbfreq)
!
! Subroutine to calculate normal mode frequencies, internal-coordinate torsion frequencies, 
! and frequencies with torsions projected out using Wilson GF method
!
use constants
implicit none
integer :: natoms,n3tm,nbl,nba,nto,nlbe,ntor,nf,nintn,nontor
integer :: i,j,lwork,info,jmin,jmax,jcol,nnegfreq
integer :: ibl(2,nf),iba(3,nf),ito(4,nf),ilbe(3,nf)
integer, allocatable :: indx(:)
double precision :: p1,p2,totmass,ratio
double precision :: coord(natoms,3),masses(natoms),t(3,3)
double precision :: fc(n3tm,n3tm),fcm(n3tm,n3tm),freq(n3tm)
double precision :: dbfreq(n3tm-ntor)
character(len=1) :: space(n3tm)

double precision, allocatable :: pfreq(:),fcp(:,:),prt(:,:)
double precision, allocatable :: eigen(:),unity(:,:),bm(:,:),gmat(:,:),fmat(:,:),gf(:,:),u(:,:)
double precision, allocatable :: alambd(:),amat(:,:),bubt(:,:),binv(:,:)
double precision, allocatable :: work(:),wr(:),wi(:),vl(:),vr(:)
double precision, allocatable :: usvd(:,:),vtsvd(:,:),ssvd(:),bmtmp(:,:)
double precision, allocatable :: gint(:,:),fint(:,:),g11(:,:),g12(:,:),g21(:,:),g22(:,:)
 
 nintn = nf - ntor 
!nintn = nf
space(:) = ' '
dbfreq(:) = 0d0

allocate ( pfreq(n3tm),fcp(n3tm,n3tm),prt(n3tm,n3tm) )
allocate ( fmat(nintn,nintn),gf(nintn,nintn),u(n3tm,n3tm),alambd(nintn),amat(n3tm,nintn))
allocate ( eigen(n3tm),unity(n3tm,n3tm),bm(nintn,n3tm),gmat(nintn,nintn))
allocate ( bubt(nintn,nintn),binv(nintn,nintn))
allocate ( indx(nintn),wr(nintn),wi(nintn),vl(nintn),vr(nintn) )
allocate (usvd(nintn,nintn),vtsvd(n3tm,n3tm),ssvd(nintn),bmtmp(nintn,n3tm))

unity(:,:) = 0.d0
do i = 1, n3tm
  unity(i,i) = 1.0d0
enddo
! The coordinates have been shifted in main subroutine,so below is commented.
!com(:) = 0.d0
!do i=1,natoms
!  com(:)=com(:)+masses(i)*coord(i,:)
!enddo
!com(:)=com(:)/totmass
! Shift input coordinates so that the center of mass (COM) is at 0,0,0; This is needed for the RT projection to follow
!do i=1,natoms
!  coord(i,:)=coord(i,:)-com(:)
!enddo

!project out translation and rotation
call gen_prt(coord,masses,totmass,natoms,n3tm,prt)
prt(:,:) = unity(:,:) - prt(:,:)
fcp(:,:) = matmul(prt(:,:),matmul(fcm(:,:),prt(:,:)))

lwork = -1
allocate (work (1) )
call dsyev('N','L',n3tm,fcp,n3tm,eigen,work,lwork,info)
lwork=int(work(1))
deallocate (work)
allocate( work(lwork) )
call dsyev('N','L',n3tm,fcp,n3tm,eigen,work,lwork,info)
deallocate (work)
if(info.ne.0) then
  write(6,*) 'Error in diagolization of force constant matrix'
  stop
endif

nnegfreq=0
p1 = 1.d0
do i = 1, n3tm
  freq(i) = sign(1.d0,eigen(i))*sqrt(abs(eigen(i)))*au2cm
! if(i.gt.6.and.freq(i).gt.1.d-1) p1 = p1*freq(i)
  if(i.gt.6.and.freq(i).gt.1.d-1) p1 = p1*freq(i)/au2cm
  if(freq(i) .lt. 0.d0) space(i)="i"
  if(freq(i) .lt. -1.d0) nnegfreq=nnegfreq+1
enddo

write(6,*)
write(6,*) "  Normal-mode Frequencies (cm-1) "
write(6,*) "  (Overall translation and rotation projected out)"
 jcol = 6
 jmin = 1
 do 
  jmax = min(jcol,n3tm)
  write(6,'(1x,6(f12.4,a))')  (abs(freq(j)),space(j),j=jmin,jmax) 
  jmin = jmax + 1
  jcol = jcol + 6
  if(jmax .eq. n3tm) exit
 enddo
!
! Use GF method to do 3N-6-t frequencies
!
! 1. Form B matrix
!
bm(:,:) = 0.d0
 call bmat(nf,nbl,nba,nto-ntor,nlbe,ibl,iba,ito,ilbe,coord(:,1),coord(:,2),coord(:,3), &
         natoms,nintn,bm,t)
!
! 2 Form G matrix and F matrix
!
u(:,:) = 0.0d0
do i = 1, natoms
   u(3*i-2,3*i-2) = 1.0d0/masses(i)
   u(3*i-1,3*i-1) = u(3*i-2,3*i-2)
   u(3*i  ,3*i  ) = u(3*i-2,3*i-2)
enddo
gmat(:,:)=matmul(bm(:,:),matmul(u(:,:),transpose(bm(:,:))))
bubt(:,:)=gmat(:,:)
!
! Option 1: following Pulay and Fogarasi, mass is introduced in A matrix
!
!   CONSTRUCT (BM*U*BMT)-1  AND A matrix
 call dgetrf(nintn,nintn,bubt,nintn,indx,info)
 if (info.ne.0 ) then
    write(6,*) "dgetrf exit abnormally"
    stop
 endif
 lwork = -1
 allocate(work(1))
 call dgetri(nintn,bubt,nintn,indx,work,lwork,info)
 lwork=int(work(1))
 deallocate (work)
 allocate( work(lwork) )
 call dgetri(nintn,bubt,nintn,indx,work,lwork,info)
 deallocate(work)
 if (info.ne.0 ) then
    write(6,*) "dgetri exit abnormally"
    stop
 endif
 amat(:,:)=matmul( (matmul(u(:,:),transpose(bm(:,:)))),bubt(:,:))
 fmat(:,:)=matmul(matmul(transpose(amat(:,:)),fc(:,:)),amat(:,:))
 gf = matmul(gmat,fmat)
 
! diagonaliz GF matrix
 lwork = -1
 allocate(work(1))
  call dgeev('N','N',nintn,gf,nintn,WR,WI,VL,nintn,VR,nintn,WORK,LWORK,INFO)
 lwork=int(work(1))
 deallocate (work)
 allocate( work(lwork) )
 call dgeev('N','N',nintn,gf,nintn,WR,WI,VL,nintn,VR,nintn,WORK,LWORK,INFO)
 deallocate (work)
 if(info.ne.0) then
    write(6,*) "dgeev exited abnormally in subroutine vib"
 endif

 space(:) = " "
 call shell(nintn,wr) !sort the frequencies
 p2 = 1.0d0
 do i = 1, nintn
!  pfreq(i) = sign(1.d0,wr(i))*sqrt(abs(wr(i)))*au2cm
   dbfreq(i)= sign(1.d0,wr(i))*sqrt(abs(wr(i)))
   pfreq(i) = dbfreq(i)*au2cm 
!  if(pfreq(i).gt.1.d-1) p2 = p2*pfreq(i)
   if(pfreq(i).gt.1.d-1) p2 = p2*pfreq(i)/au2cm
   if(pfreq(i).lt.0.d0) space(i)="i"
 enddo

write(6,*) "  Torsion-projected frequencies/non-torsional frequencies (cm-1)"
 jcol = 6
 jmin = 1
 do 
  jmax = min(jcol,nintn)
  write(6,'(1x,6(f12.4,a))')  (abs(pfreq(j)),space(j),j=jmin,jmax) 
  jmin = jmax + 1
  jcol = jcol + 6
  if(jmax .eq. nintn) exit
 enddo
 ratio = p1/p2
end subroutine vib

subroutine gen_prt(coord,masses,totmass,natoms,n3tm,prt)
!
! Project out translations and instantaneous rotations from Hessian
! coord is assumed to be centered at the center of mass
!
implicit none
integer :: natoms,n3tm,i,j,ip,jp,ia,ib,ja,jb,ic,jc,ix,jx,jend
double precision :: masses(natoms),coord(natoms,3),prt(n3tm,n3tm)
double precision :: rot(3,3),tensor(3,3,3),totmass,summ

DATA TENSOR / 5*0.0D0,-1.0D0,0.0D0,1.0D0,3*0.0D0,1.0D0,3*0.0D0,  &
                  -1.0D0,3*0.0D0,-1.0D0,0.0D0,1.0D0,5*0.0D0 /

! Calculate the moment of inertia 
!
rot(:,:) = 0.0d0
prt(:,:) = 0.0d0
do i=1,natoms
  rot(1,1)= rot(1,1) + masses(i)*(coord(i,2)**2+coord(i,3)**2)
  rot(2,2)= rot(2,2) + masses(i)*(coord(i,1)**2+coord(i,3)**2)
  rot(3,3)= rot(3,3) + masses(i)*(coord(i,1)**2+coord(i,2)**2)
  rot(1,2)= rot(1,2) - masses(i)*coord(i,1)*coord(i,2)
  rot(1,3)= rot(1,3) - masses(i)*coord(i,1)*coord(i,3)
  rot(2,3)= rot(2,3) - masses(i)*coord(i,2)*coord(i,3)
enddo
rot(2,1)=rot(1,2)
rot(3,1)=rot(1,3)
rot(3,2)=rot(2,3)

call inv3(rot) ! overwrite rot with its inverse

 do ip = 1, natoms
   ix = 3*(ip-1)
 do jp = 1, ip
   jx = 3*(jp-1)
   do ic = 1, 3
      JEND = 3
      IF (JP.EQ.IP) JEND = IC
   do jc = 1, JEND

        summ = 0.d0
        DO IA = 1, 3
        DO IB = 1, 3
           IF (TENSOR(IA,IB,IC).eq.0d0)cycle
           DO JA = 1, 3
           DO JB = 1, 3
              IF (TENSOR(JA,JB,JC).eq.0d0)cycle
              SUMM = SUMM+TENSOR(IA,IB,IC)*TENSOR(JA,JB,JC)*ROT(IA,JA)* &
                     COORD(IP,IB)*COORD(JP,JB)*sqrt(masses(ip)*masses(jp))
           enddo     
           enddo     
        enddo    
        enddo    

     i = ix + ic
     j = jx + jc
     prt(i,j)= summ  
     if (ic.eq.jc) prt(i,j) = prt(i,j) + sqrt(masses(ip)*masses(jp))/totmass
     if (i.ne.j)   prt(j,i) = prt(i,j)
   enddo ! end jc
   enddo ! end ic
 enddo ! end jp
 enddo ! end IP
end

subroutine inv3(m)
!
! subroutine to inverse a 3x3 matrix
!
implicit none
double precision :: m(3,3),w(3,3),z
w(1,1) = m(2,2)*m(3,3)-m(2,3)*m(3,2)
w(1,2) = m(1,3)*m(3,2)-m(1,2)*m(3,3)
w(1,3) = m(1,2)*m(2,3)-m(1,3)*m(2,2)
w(2,1) = m(2,3)*m(3,1)-m(2,1)*m(3,3)
w(2,2) = m(1,1)*m(3,3)-m(1,3)*m(3,1)
w(2,3) = m(1,3)*m(2,1)-m(1,1)*m(2,3)
w(3,1) = m(2,1)*m(3,2)-m(2,2)*m(3,1)
w(3,2) = m(1,2)*m(3,1)-m(1,1)*m(3,2)
w(3,3) = m(1,1)*m(2,2)-m(1,2)*m(2,1)
z = m(1,1)*w(1,1) + m(1,2)*w(2,1) + m(1,3)*w(3,1)
m(:,:) = w(:,:)/z
end subroutine inv3

subroutine transfc(coord,masses,natoms,n3tm,nf,fc,ibl,iba,ito,ilbe,nbl,nba,nto,&
                   nlbe,ntor,fcint,rmi,mtor)
!
! Subroutine to transform Hessian in Cartesian coordinates into Hessian in internal 
! coordinate and to calculate bar omega 
!
use constants
implicit none
integer :: natoms,n3tm,nf,ntor,nto,nbl,nba,nlbe
integer :: ibl(2,nf),iba(3,nf),ito(4,nf),ilbe(3,nf),indx(nf)
integer :: i,j,lwork,info,jmin,jmax,jcol
double precision :: coord(natoms,3)
double precision :: masses(natoms),fc(n3tm,n3tm)
double precision :: fcint(nf),rmi(ntor),t(3,3)
double precision :: fmat(nf,nf),u(n3tm,n3tm),bm(nf,n3tm),amat(n3tm,nf)
double precision :: gmat(nf,nf),bubt(nf,nf),bomega(ntor),w(ntor),mtor(ntor)
double precision :: tmp
double precision, allocatable :: work(:)

bm(:,:) = 0.d0

call bmat(nf,nbl,nba,nto,nlbe,ibl,iba,ito,ilbe,coord(:,1),coord(:,2),coord(:,3), &
          natoms,nf,bm,t)
u(:,:) = 0.0d0
do i = 1, natoms
   u(3*i-2,3*i-2) = 1.0d0/masses(i)
   u(3*i-1,3*i-1) = u(3*i-2,3*i-2)
   u(3*i  ,3*i  ) = u(3*i-2,3*i-2)
enddo
gmat(:,:)=matmul(bm(:,:),matmul( u(:,:),transpose(bm(:,:)) ) )
bubt(:,:)=gmat(:,:)

 call dgetrf(nf,nf,bubt,nf,indx,info)
 if (info.ne.0 ) then
    write(6,*) "dgetrf exit abnormally in transfc"
    stop
 endif
 lwork = -1

 allocate(work(1))
 call dgetri(nf,bubt,nf,indx,work,lwork,info)
 lwork=int(work(1))
 deallocate (work)
 allocate( work(lwork) )
 call dgetri(nf,bubt,nf,indx,work,lwork,info)
 if (info.ne.0 ) then
    write(6,*) "dgetri exit abnormally in transfc"
    stop
 endif

 amat(:,:)=matmul( (matmul(u(:,:),transpose(bm(:,:)))),bubt(:,:))
 fmat(:,:)=matmul(matmul(transpose(amat(:,:)),fc(:,:)),amat(:,:))

 tmp = 1.d0
 do i = 1, ntor 
    j = nf-ntor+i
    fcint(i) = fmat(j,j)
    bomega(i) = sqrt(fcint(i)/rmi(i))
    w(i) = 2.d0*fcint(i)/(mtor(i)**2)
!
 enddo
!   write(6,*) 'product of torsion ', tmp
!
 
 write(6,900)
900 format(/10X,'Uncoupled force constants (au), Pitzer moment of inertia (amu A^2),',/ &
          3X,'uncoupled torsional frequencies (cm-1), and uncoupled torsional barrier (kcal/mol)')
 write(6,910)
910 format(3X,80('-'))

jmin =1
jcol =5
do 
 jmax = min(jcol,ntor)
 write(6,'(3X,a,14x,10(i2,12X))')  'Torsion', (j, j =jmin ,jmax)
 write(6,'(3X,a,1X,10(1p,e13.5))') 'Force constant', (fcint(j),j=jmin,jmax)
 write(6,'(3X,a,8X,10(f8.3,5X))')  'Moment',(rmi(j)*au2amu,j=jmin,jmax)
 write(6,'(3X,a,6X,10(f8.3,5X))')  'Tor. Freq.', (bomega(j)*au2cm,j=jmin,jmax)
 write(6,'(3X,a,4X,10(f8.3,5X))')  'Rot. barrier', (w(j)*aukcal,j=jmin,jmax)
 write(6,'(3X,a,4X,10(f8.3,5X))')  'M_j,tau     ', (mtor(j),j=jmin,jmax)
 jmin = jmax + 1
 jcol = jcol + 4
 if(jmax == ntor) exit
 write(6,*)
enddo

IF(ALLOCATED(work)) DEALLOCATE(work) 
end subroutine transfc

subroutine shell(n,arr)
!
! This is a simple variant (more sophisticated versions exist) of the Shell sorting algorithm.
! It is a convenient choice here because we are sorting limited data that is already nearly well ordered.
! For large data sets another code should be chosen.
! 
implicit none
integer :: n, nn, lognb2, m, k, j, l, i
double precision, parameter :: ALN2I=1.d0/0.69314718d0,TINY=1.d-5 
double precision :: ARR(N), t

LOGNB2=INT(ALOG(FLOAT(N))*ALN2I+TINY)
m=n
do nn=1,LOGNB2
  m=m/2; k=n-m
  do j=1,k
    i=j
    10    continue
    l=i+m
    if(ARR(l).LT.ARR(i)) then
      t=ARR(i)
      ARR(i)=ARR(l)
      ARR(l)=t
      i=i-m
      if(i.GE.1) GOTO 10
    end if
  end do
end do
end subroutine shell


subroutine vibtor(coord,masses,totmass,fc,natoms,n3tm,nf,ibl,iba,ito,ilbe,nbl,nba,nto,nlbe,ntor, &
                  rmi, mtor, bh, dmat,detd)
!
! Subroutine to calculate torsional frequencies by internal coordinates
!
use constants
implicit none
integer :: natoms,n3tm,nbl,nba,nto,nlbe,ntor,nf,nintn,nontor
integer :: i,j,lwork,info,jmin,jmax,jcol
integer :: ibl(2,nf),iba(3,nf),ito(4,nf),ilbe(3,nf)
integer, allocatable :: indx(:)
double precision :: totmass,avem,norm,detd,prod_pcfreq
double precision :: coord(natoms,3),masses(natoms),com(3),t(3,3)
double precision :: fc(n3tm,n3tm)
double precision :: bfreq(ntor)
double precision :: rmi(ntor),mtor(ntor),dmat(ntor,ntor),bh(ntor),dmat2(ntor,ntor)
double precision :: dtmp(ntor,ntor),deigen(ntor)
character(len=1) :: space(n3tm)

double precision, allocatable :: pfreq(:),gftor(:,:),hw(:,:)
double precision, allocatable :: eigen(:),unity(:,:),bm(:,:),gmat(:,:),fmat(:,:),gf(:,:),u(:,:)
double precision, allocatable :: alambd(:),amat(:,:),bubt(:,:),binv(:,:)
double precision, allocatable :: work(:),wr(:),wi(:),vl(:),vr(:)
double precision, allocatable :: wrt(:),wit(:),vlt(:,:),vrt(:,:),gtor(:,:),ftor(:,:)
double precision, allocatable :: g12(:,:),g21(:,:),g11(:,:),g22(:,:)
double precision, allocatable :: gfnontor(:,:),gnontor(:,:),fnontor(:,:),hm(:,:),vb(:,:)
double precision, allocatable :: usvd(:,:),vtsvd(:,:),ssvd(:),bmtmp(:,:)
 
nintn =  nf 
nontor = nf - ntor 
space(:) = ' '
bfreq(:) = 0d0

allocate ( pfreq(n3tm),gftor(ntor,ntor),hw(ntor,ntor),gfnontor(nontor,nontor))
allocate ( fmat(nintn,nintn),gf(nintn,nintn),u(n3tm,n3tm),alambd(nintn),amat(n3tm,nintn))
allocate ( eigen(n3tm),bm(nintn,n3tm),gmat(nintn,nintn))
allocate ( bubt(nintn,nintn),binv(nintn,nintn))
allocate ( indx(nintn),wr(nintn),wi(nintn),vl(nintn),vr(nintn) )
allocate (wrt(ntor),wit(ntor),vlt(ntor,ntor),vrt(ntor,ntor))
allocate (gtor(ntor,ntor),ftor(ntor,ntor))
allocate (g11(nontor,nontor),g22(ntor,ntor),g12(nontor,ntor),g21(ntor,nontor))
allocate (gnontor(nontor,nontor),fnontor(nontor,nontor))
allocate (hm(nf,nf),vb(nf,nf))
allocate (usvd(nintn,nintn),vtsvd(n3tm,n3tm),ssvd(nintn),bmtmp(nintn,n3tm))

! Use GF method to do t torsional frequencies
!
! 1. Form B matrix
!
bm(:,:) = 0.d0
call bmat(nf,nbl,nba,nto,nlbe,ibl,iba,ito,ilbe,coord(:,1),coord(:,2),coord(:,3), &
         natoms,nintn,bm,t)
!
! 2 Form G matrix and F matrix
!
u(:,:) = 0.0d0
do i = 1, natoms
   u(3*i-2,3*i-2) = 1.0d0/masses(i)
   u(3*i-1,3*i-1) = u(3*i-2,3*i-2)
   u(3*i  ,3*i  ) = u(3*i-2,3*i-2)
enddo
gmat(:,:)=matmul(bm(:,:),matmul(u(:,:),transpose(bm(:,:))))
bubt(:,:)=gmat(:,:)
!
! Option 1: following Pulay and Fogarasi, mass is introduced in A matrix
! so that translation and rotation can be projected out
!
!   CONSTRUCT (BM*U*BMT)-1  AND A matrix
!
 call dgetrf(nintn,nintn,bubt,nintn,indx,info)
 if (info.ne.0 ) then
    write(6,*) "dgetrf exit abnormally"
    stop
 endif
 lwork = -1
 allocate(work(1))
 call dgetri(nintn,bubt,nintn,indx,work,lwork,info)
 lwork=int(work(1))
 deallocate (work)
 allocate( work(lwork) )
 call dgetri(nintn,bubt,nintn,indx,work,lwork,info)
 deallocate(work)
 if (info.ne.0 ) then
    write(6,*) "dgetri exit abnormally"
    stop
 endif
 amat(:,:)=matmul((matmul(u(:,:),transpose(bm(:,:)))),bubt(:,:))
 fmat(:,:)=matmul(matmul(transpose(amat(:,:)),fc(:,:)),amat(:,:))

 g11(:,:) = gmat(1:nontor,1:nontor)
 g22(:,:) = gmat(nontor+1:nf,nontor+1:nf)
 g12(:,:) = gmat(1:nontor,nontor+1:nf)
 g21(:,:) = gmat(nontor+1:nf,1:nontor)
 ftor(:,:)= fmat(nontor+1:nf,nontor+1:nf)

! Calculate (G11)-1 
 call dgetrf(nontor,nontor,g11,nontor,indx,info)
 if (info.ne.0 ) then
    write(6,*) "dgetrf exit abnormally"
    stop
 endif
 lwork = -1
 allocate(work(1))
 call dgetri(nontor,g11,nontor,indx,work,lwork,info)
 lwork=int(work(1))
 deallocate (work)
 allocate( work(lwork) )
 call dgetri(nontor,g11,nontor,indx,work,lwork,info)
 deallocate(work)
 if (info.ne.0 ) then
    write(6,*) "dgetri exit abnormally"
    stop
 endif
 
 gtor(:,:) = g22(:,:) - matmul(g21(:,:),matmul(g11(:,:),g12(:,:)  ) ) 
 gftor(:,:)=matmul(gtor(:,:),ftor(:,:))
 
 nintn=ntor
 lwork = -1
 allocate(work(1))
  call dgeev('N','N',ntor,gftor,ntor,WRT,WIT,VLT,ntor,VRT,ntor,WORK,LWORK,INFO)
 lwork=int(work(1))
 deallocate (work)
 allocate( work(lwork) )
 call dgeev('N','V',ntor,gftor,ntor,WRT,WIT,VLT,ntor,VRT,ntor,WORK,LWORK,INFO)
 deallocate (work)
 if(info.ne.0) then
    write(6,*) "dgeev exited abnormally in subroutine vibtor"
 endif

 prod_pcfreq=1d0
 space(:) = " "
 call shell(nintn,wrt) !sort the frequencies
 do i = 1, nintn
   bfreq(i)= sign(1.d0,wrt(i))*sqrt(abs(wrt(i)))
   pfreq(i) = bfreq(i)*au2cm 
   prod_pcfreq=prod_pcfreq*pfreq(i)
   if(pfreq(i).lt.0.d0) space(i)="i"
 enddo
    
 write(6,'(3X,a)')  "Primitively coupled torsional frequencies (cm-1) (Not used)"
 jcol = 6
 jmin = 1
 do 
  jmax = min(jcol,nintn)
  write(6,'(1x,6(f12.4,a))')  (abs(pfreq(j)),space(j),j=jmin,jmax) 
  jmin = jmax + 1
  jcol = jcol + 6
  if(jmax .eq. nintn) exit
 enddo
 write(6,'(3X,a,es12.4)') 'Product of primitively coupled torsional frequencies (in cm-1^ntor) ', &
       prod_pcfreq

!!inverse Gtor matrix
 g22(:,:) = gtor(:,:)
 call dgetrf(ntor,ntor,g22,ntor,indx,info)
 if (info.ne.0 ) then
    write(6,*) "dgetrf exit abnormally"
    stop
 endif
 lwork = -1
 allocate(work(1))
 call dgetri(ntor,g22,ntor,indx,work,lwork,info)
 lwork=int(work(1))
 deallocate (work)
 allocate( work(lwork) )
 call dgetri(ntor,g22,ntor,indx,work,lwork,info)
 deallocate(work)
 if (info.ne.0 ) then
    write(6,*) "dgetri exit abnormally"
    stop
 endif
 dmat(:,:) = g22(:,:)
!calculate determinant of D matrix
 lwork = -1
 allocate(work(1))
  call dgeev('N','N',ntor,g22,ntor,WRT,WIT,VLT,ntor,VRT,ntor,WORK,LWORK,INFO)
 lwork=int(work(1))
 deallocate (work)
 allocate( work(lwork) )
 call dgeev('N','V',ntor,g22,ntor,WRT,WIT,VLT,ntor,VRT,ntor,WORK,LWORK,INFO)
 deallocate (work)
 if(info.ne.0) then
    write(6,*) "dgeev exited abnormally in diagonalizing D matrix"
 endif
 deigen(:) = wrt(:)
 detd = 1d0
 do i = 1, ntor
   detd = detd*deigen(i)
   rmi(i) = dmat(i,i)
 enddo

  do i = 1, ntor
  do j = 1, ntor
    ftor(i,j) = fmat(i+nontor,j+nontor)/(mtor(i)*mtor(j))  
  enddo
  enddo
  hw(:,:) = ftor(:,:)

! solve the eigenvalues of hw matrix
 lwork = -1
 allocate(work(1))
  call dgeev('N','N',ntor,hw,ntor,WRT,WIT,VLT,ntor,VRT,ntor,WORK,LWORK,INFO)
 lwork=int(work(1))
 deallocate (work)
 allocate( work(lwork) )
 call dgeev('N','V',ntor,hw,ntor,WRT,WIT,VLT,ntor,VRT,ntor,WORK,LWORK,INFO)
 deallocate (work)
 if(info.ne.0) then
    write(6,*) "dgeev exited abnormally in subroutine vib"
 endif

 space(:) = " "
 call shell(nintn,wrt) !sort the eigenvalues
 do i = 1, nintn
   bh(i)= 2d0*wrt(i)
   if(bh(i) < 0d0) then
    write(6,*) 'Negative torsional barrier'
    stop
   endif
   if(wit(i).ne.0d0) write(6,*) 'Complex eigenvalue for Tor. barrier' 
 enddo

 write(6,'(3X,a)') "Effective coupled torsional barrier (kcal/mol)"
 jcol = 6
 jmin = 1
 do
  jmax = min(jcol,nintn)
  write(6,'(1x,6(f12.3,a))')  ((bh(j))*aukcal,space(j),j=jmin,jmax)
  jmin = jmax + 1
  jcol = jcol + 6
  if(jmax .eq. nintn) exit
 enddo

   
 write(6,'(/3X,a)') 'D matrix: '
 dmat = dmat*cau/a2bohr**2
 call matprint (6,dmat,ntor,ntor,ntor)
 write(6,'(3X,a)') 'Eigenvalues of D matrix:'
 jcol = 6
 jmin = 1
 do
  jmax = min(jcol,nintn)
  write(6,'(1x,6es14.4)')  (deigen(j)*cau/a2bohr**2,j=jmin,jmax)
  jmin = jmax + 1
  jcol = jcol + 6
  if(jmax .eq. nintn) exit
 enddo
 return
end subroutine vibtor

