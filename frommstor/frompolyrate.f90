!***********************************************************************
!     BMAT (From POLYRATE)
!     modified by Jingjing Zheng Sep. 21 2010
!              by Shaohong L. Li Sep.  9 2015: added out-of-plane bend
!              by Shaohong L. Li Oct.  ? 2015: added out-of-plane distance
!***********************************************************************
!
!SHLL 1509
SUBROUTINE bmat(NF,NBL,NBA,NTO,NOB,NOD,NLBE,IBL,IBA,ITO,IOB,IOD,ILBE &
     ,X,Y,Z,NATOM,NINT,BM,T)
!original
!SUBROUTINE bmat(NF,NBL,NBA,NTO,NLBE,IBL,IBA,ITO,ILBE,X,Y,Z,NATOM,NINT,BM,T)
!SHLL end
!
!     CALCULATES THE MATRIX DR/DX (B matrix)
!
implicit none
!GB nbl=#bond len;nba=#bond angle;nto=#torsion;nlbe=#linear bend;
!GB natom=#atom;nint=#internal coord
!GB nf should be some large # that defines the dim of array
integer :: irow, ij, nf, nbl, nba, nto, nlbe, natom, nint
integer :: i, j, ijk, im, k, l, ijkl
!GB IBL(n,i)=the n-th atom involved in i-th bond len 
integer :: IBL(2,NF),IBA(3,NF),ITO(4,NF),ILBE(3,NF)
double precision :: dist, angl, cose, sine
double precision :: BM(NINT,3*NATOM), X(NATOM),Y(NATOM),Z(NATOM)
double precision :: PX(3),PY(3),PZ(3),V(3),T(3,3),BMTE(3)
double precision :: rij, rij1, rij2, rij3, rjk1, rjk2, rjk3, rjk, pijk
double precision :: cijk, sijk, pxmg, pymg, pzmg, sjkl, sijkl
double precision :: rkl, rkl1, rkl2, rkl3, cjkl, cijkl, pjkl
!SHLL 1509
integer :: nob,IOB(4,NF),nod,IOD(4,NF)
double precision :: A243,S243,C243,E41(3),E42(3),E43(3) &
     ,R41,R42,R43,OOPB,COOPB,TOOPB,TMPV(3),BMV1(3),BMV2(3),BMV3(3) &
     ,BMV4(3),e12(3),e13(3),v14(3),r12,r13,theta,tmpv2(3) &
     ,tmpv3(3),rtmp
!SHLL end
!
IROW=0
!GB loop over # bond len
DO IJ=1,NBL
  IROW=IROW+1
!GB I,J=indices of two atoms of bond IJ
  I=IBL(1,IJ)
  J=IBL(2,IJ)
  RIJ1=X(J)-X(I)
  RIJ2=Y(J)-Y(I)
  RIJ3=Z(J)-Z(I)
!GB dist()=distance b/n atoms I J
  RIJ=DIST(I,J,X,Y,Z,NATOM)
!GB B mat for bond len = unit vector along bond
  BM(IROW,3*I-2)=-RIJ1/RIJ
  BM(IROW,3*I-1)=-RIJ2/RIJ
  BM(IROW,3*I)=-RIJ3/RIJ
  BM(IROW,3*J-2)=RIJ1/RIJ
  BM(IROW,3*J-1)=RIJ2/RIJ
  BM(IROW,3*J)=RIJ3/RIJ
enddo    
!
!GB Loop over bond angles
DO IJK=1,NBA
  IROW=IROW+1
!GB indices of three atoms
  I=IBA(1,IJK)
  J=IBA(2,IJK)
  K=IBA(3,IJK)
!GB vectors I->J, J->K
  RIJ1=X(J)-X(I)
  RIJ2=Y(J)-Y(I)
  RIJ3=Z(J)-Z(I)
  RJK1=X(K)-X(J)
  RJK2=Y(K)-Y(J)
  RJK3=Z(K)-Z(J)
!GB length of vectors
  RIJ=DIST(I,J,X,Y,Z,NATOM)
  RJK=DIST(J,K,X,Y,Z,NATOM)
!GB bond angle and its sine and cosine
  PIJK=ANGL(I,J,K,X,Y,Z,NATOM)
  CIJK=COSE(PIJK)
  SIJK=SINE(PIJK)
!GB B mat for bond angle - same as in my GB:Internal Coordinates.doc
  BM(IROW,3*I-2)=(-CIJK*RIJ1*RJK-RIJ*RJK1)/(SIJK*RIJ**2*RJK)
  BM(IROW,3*I-1)=(-CIJK*RIJ2*RJK-RIJ*RJK2)/(SIJK*RIJ**2*RJK)
  BM(IROW,3*I)=(-CIJK*RIJ3*RJK-RIJ*RJK3)/(SIJK*RIJ**2*RJK)
  BM(IROW,3*J-2)= (-RIJ*RIJ1*RJK+CIJK*RIJ1*RJK**2-CIJK*RIJ**2*RJK1+RIJ*RJK*RJK1) &
       /(SIJK*RIJ**2*RJK**2)
  BM(IROW,3*J-1)= (-RIJ*RIJ2*RJK+CIJK*RIJ2*RJK**2-CIJK*RIJ**2*RJK2+RIJ*RJK*RJK2) &
       /(SIJK*RIJ**2*RJK**2)
  BM(IROW,3*J)= (-RIJ*RIJ3*RJK+CIJK*RIJ3*RJK**2-CIJK*RIJ**2*RJK3+RIJ*RJK*RJK3)  &
       /(SIJK*RIJ**2*RJK**2)
  BM(IROW,3*K-2)=(RIJ1*RJK+CIJK*RIJ*RJK1)/(SIJK*RIJ*RJK**2)
  BM(IROW,3*K-1)=(RIJ2*RJK+CIJK*RIJ*RJK2)/(SIJK*RIJ*RJK**2)
  BM(IROW,3*K)=(RIJ3*RJK+CIJK*RIJ*RJK3)/(SIJK*RIJ*RJK**2)
enddo    
!
!
!     B-matrix for linear bend  - added 09xxYC96
!
DO IJK=1,NLBE
  IROW=IROW+1
  I=ILBE(1,IJK)
  J=ILBE(2,IJK)
  K=ILBE(3,IJK) 
  RIJ=DIST(I,J,X,Y,Z,NATOM)
  RJK=DIST(J,K,X,Y,Z,NATOM)
!
!  remember C-A-B ==> K-J-I : cases of z-axis
!
!  choice #1 : Z axis is assigned as rjk
!
!      PZ(1)=X(K)-X(J)
!      PZ(2)=Y(K)-Y(J)
!      PZ(3)=Z(K)-Z(J)
!
!  choice #2 : Z axis is assigned as rij
!
!      PZ(1)=X(J)-X(I)
!      PZ(2)=Y(J)-Y(I)
!      PZ(3)=Z(J)-Z(I)
!
!  choice #3 : Z axis is assigned as rik
!
  PZ(1)=X(K)-X(I)
  PZ(2)=Y(K)-Y(I)
  PZ(3)=Z(K)-Z(I)
!
  IF (ABS(PZ(3)/RJK).NE.1) THEN
!
!     set up local coordinates
!
    V(1)=1.0d0
    V(2)=1.0d0
    V(3)=1.0d0
!
!     find PY perpendicular to PZ and V
!
    CALL XPROD(V,PZ,PY)
!
!     find PX perpendicular to PY and PZ
!
    CALL XPROD(PY,PZ,PX)
!
!     normalize
!
    PXMG=0.0d0
    PYMG=0.0d0
    DO IM=1,3
      PXMG=PXMG+PX(IM)*PX(IM)
      PYMG=PYMG+PY(IM)*PY(IM)
    ENDDO
    PXMG=SQRT(PXMG)
    PYMG=SQRT(PYMG)
    PZMG=RJK
!
!   set up transpose of transformation matrix
!
    DO IM=1,3
      T(1,IM)=PX(IM)/PXMG
      T(2,IM)=PY(IM)/PYMG
      T(3,IM)=PZ(IM)/PZMG
    ENDDO
  ELSE
    DO IM=1,3
      IF(PZ(3)/PZMG.LT.0) THEN
         T(IM,IM)=-1.0d0
      ELSE
         T(IM,IM)=1.0d0
      ENDIF
    ENDDO
  ENDIF
!
!   the Ry mode ; Califano + Mathematica
!
!   in Califano the molecule in C-A-B
!   here is                     K-J-I
!
  BMTE(1)= 1.0d0/RIJ 
  BMTE(2)= -(RIJ+RJK)/(RIJ*RJK)
  BMTE(3)= 1.0d0/RJK
!
!   transform back to the original coordinates
!     B = B'T(transpose)
!
  BM(IROW,3*I-2)=BMTE(1)*T(1,1)
  BM(IROW,3*I-1)=BMTE(1)*T(1,2)
  BM(IROW,3*I)  =BMTE(1)*T(1,3)
  BM(IROW,3*J-2)=BMTE(2)*T(1,1)
  BM(IROW,3*J-1)=BMTE(2)*T(1,2)
  BM(IROW,3*J)  =BMTE(2)*T(1,3)
  BM(IROW,3*K-2)=BMTE(3)*T(1,1)
  BM(IROW,3*K-1)=BMTE(3)*T(1,2)
  BM(IROW,3*K)  =BMTE(3)*T(1,3)
!
!   the Rx mode 
!
  IROW=IROW+1
  BM(IROW,3*I-2)=BMTE(1)*T(2,1)
  BM(IROW,3*I-1)=BMTE(1)*T(2,2)
  BM(IROW,3*I)  =BMTE(1)*T(2,3)
  BM(IROW,3*J-2)=BMTE(2)*T(2,1)
  BM(IROW,3*J-1)=BMTE(2)*T(2,2)
  BM(IROW,3*J)  =BMTE(2)*T(2,3)
  BM(IROW,3*K-2)=BMTE(3)*T(2,1)
  BM(IROW,3*K-1)=BMTE(3)*T(2,2)
  BM(IROW,3*K)  =BMTE(3)*T(2,3)
enddo    
!
!  Torsion
!
!GB torsion
DO IJKL=1,NTO
  IROW=IROW+1
!GB indices of atoms
  I=ITO(1,IJKL)
  J=ITO(2,IJKL)
  K=ITO(3,IJKL)
  L=ITO(4,IJKL)
!GB vectors I->J, J->K, K->L and their length
  RIJ1=X(J)-X(I)
  RIJ2=Y(J)-Y(I)
  RIJ3=Z(J)-Z(I)
  RJK1=X(K)-X(J)
  RJK2=Y(K)-Y(J)
  RJK3=Z(K)-Z(J)
  RKL1=X(L)-X(K)
  RKL2=Y(L)-Y(K)
  RKL3=Z(L)-Z(K)
  RIJ=DIST(I,J,X,Y,Z,NATOM)
  RJK=DIST(J,K,X,Y,Z,NATOM)
  RKL=DIST(K,L,X,Y,Z,NATOM)
!GB two bond angles IJK and JKL and their sine and cosine
  PIJK=ANGL(I,J,K,X,Y,Z,NATOM)
  CIJK=COSE(PIJK)
  SIJK=SINE(PIJK)
  PJKL=ANGL(J,K,L,X,Y,Z,NATOM)
  CJKL=COSE(PJKL)
  SJKL=SINE(PJKL)
  CIJKL=((-RIJ2*RJK1+RIJ1*RJK2)*(-RJK2*RKL1+RJK1*RKL2)+ &
       (RIJ3*RJK1-RIJ1*RJK3)*(RJK3*RKL1-RJK1*RKL3)+     &
       (-RIJ3*RJK2+RIJ2*RJK3)*(-RJK3*RKL2+RJK2*RKL3))/  &
       (SIJK*SJKL*RIJ*RJK*RJK*RKL)
  SIJKL=((-RIJ3*RJK2+RIJ2*RJK3)*RKL1+(RIJ3*RJK1-RIJ1*RJK3)*RKL2+  &
       (-(RIJ2*RJK1)+RIJ1*RJK2)*RKL3)/(RIJ*RJK*RKL*SIJK*SJKL)
  BM(IROW,3*I-2)=SIJKL*RIJ1/(CIJKL*SIJK**2*RIJ**2)+ &
       CIJK*SIJKL*RJK1/(CIJKL*SIJK**2*RIJ*RJK)+     &
       (RJK3*RKL2-RJK2*RKL3)/(CIJKL*SIJK*SJKL*RIJ*RJK*RKL)
  BM(IROW,3*I-1)=SIJKL*RIJ2/(CIJKL*SIJK**2*RIJ**2)+ &
       CIJK*SIJKL*RJK2/(CIJKL*SIJK**2*RIJ*RJK)+     &
       (-RJK3*RKL1+RJK1*RKL3)/(CIJKL*SIJK*SJKL*RIJ*RJK*RKL)
  BM(IROW,3*I)=SIJKL*RIJ3/(CIJKL*SIJK**2*RIJ**2)+ &
       CIJK*SIJKL*RJK3/(CIJKL*SIJK**2*RIJ*RJK)+   &
       (RJK2*RKL1-RJK1*RKL2)/(CIJKL*SIJK*SJKL*RIJ*RJK*RKL)
  BM(IROW,3*J-2)=-(SIJKL*RIJ1/(CIJKL*SIJK**2*RIJ**2))+  &
       CIJK*SIJKL*(RIJ1-RJK1)/(CIJKL*SIJK**2*RIJ*RJK)-  &
       SIJKL*RJK1/(CIJKL*RJK**2)+SIJKL*RJK1/(CIJKL*SIJK**2*RJK**2)+  &
       SIJKL*RJK1/(CIJKL*SJKL**2*RJK**2)+  &
       CJKL*SIJKL*RKL1/(CIJKL*SJKL**2*RJK*RKL)+  &
       (-RIJ3*RKL2-RJK3*RKL2+RIJ2*RKL3+RJK2*RKL3)/  &
       (CIJKL*SIJK*SJKL*RIJ*RJK*RKL)
  BM(IROW,3*J-1)=-SIJKL*RIJ2/(CIJKL*SIJK**2*RIJ**2)+  &
       CIJK*SIJKL*(RIJ2-RJK2)/(CIJKL*SIJK**2*RIJ*RJK)- & 
       SIJKL*RJK2/(CIJKL*RJK**2)+SIJKL*RJK2/(CIJKL*SIJK**2*RJK**2)+  &
       SIJKL*RJK2/(CIJKL*SJKL**2*RJK**2)+ &
       CJKL*SIJKL*RKL2/(CIJKL*SJKL**2*RJK*RKL)+  &
       (RIJ3*RKL1+RJK3*RKL1-RIJ1*RKL3-RJK1*RKL3)/  & 
       (CIJKL*SIJK*SJKL*RIJ*RJK*RKL) 
  BM(IROW,3*J)=-SIJKL*RIJ3/(CIJKL*SIJK**2*RIJ**2)+   &
       CIJK*SIJKL*(RIJ3-RJK3)/(CIJKL*SIJK**2*RIJ*RJK)-  &
       SIJKL*RJK3/(CIJKL*RJK**2)+SIJKL*RJK3/(CIJKL*SIJK**2*RJK**2)+  &
       SIJKL*RJK3/(CIJKL*SJKL**2*RJK**2)+   &
       (-RIJ2*RKL1-RJK2*RKL1+RIJ1*RKL2+RJK1*RKL2)/  &
       (CIJKL*SIJK*SJKL*RIJ*RJK*RKL)+   &
       CJKL*SIJKL*RKL3/(CIJKL*SJKL**2*RJK*RKL)   
  BM(IROW,3*K-2)=-CIJK*SIJKL*RIJ1/(CIJKL*SIJK**2*RIJ*RJK)+  &
       SIJKL*RJK1/(CIJKL*RJK**2)-SIJKL*RJK1/(CIJKL*SIJK**2*RJK**2)-  &
       SIJKL*RJK1/(CIJKL*SJKL**2*RJK**2)+  &
       CJKL*SIJKL*(RJK1-RKL1)/(CIJKL*SJKL**2*RJK*RKL)+  &
       SIJKL*RKL1/(CIJKL*SJKL**2*RKL**2)+  &
       (RIJ3*RJK2-RIJ2*RJK3+RIJ3*RKL2-RIJ2*RKL3)/  &
       (CIJKL*SIJK*SJKL*RIJ*RJK*RKL)
  BM(IROW,3*K-1)=-CIJK*SIJKL*RIJ2/(CIJKL*SIJK**2*RIJ*RJK)+  &
       SIJKL*RJK2/(CIJKL*RJK**2)-SIJKL*RJK2/(CIJKL*SIJK**2*RJK**2)-  &
       SIJKL*RJK2/(CIJKL*SJKL**2*RJK**2)+  &
       CJKL*SIJKL*(RJK2-RKL2)/(CIJKL*SJKL**2*RJK*RKL)+  &
       SIJKL*RKL2/(CIJKL*SJKL**2*RKL**2)+  &
       (-RIJ3*RJK1+RIJ1*RJK3-RIJ3*RKL1+RIJ1*RKL3)/  &
       (CIJKL*SIJK*SJKL*RIJ*RJK*RKL)
  BM(IROW,3*K)=-CIJK*SIJKL*RIJ3/(CIJKL*SIJK**2*RIJ*RJK)+  &
       SIJKL*RJK3/(CIJKL*RJK**2)-SIJKL*RJK3/(CIJKL*SIJK**2*RJK**2)-  &
       SIJKL*RJK3/(CIJKL*SJKL**2*RJK**2)+  &
       (RIJ2*RJK1-RIJ1*RJK2+RIJ2*RKL1-RIJ1*RKL2)/  &
       (CIJKL*SIJK*SJKL*RIJ*RJK*RKL)+  &
       CJKL*SIJKL*(RJK3-RKL3)/(CIJKL*SJKL**2*RJK*RKL)+  &
       SIJKL*RKL3/(CIJKL*SJKL**2*RKL**2)
  BM(IROW,3*L-2)=-CJKL*SIJKL*RJK1/(CIJKL*SJKL**2*RJK*RKL)+ &
       (-RIJ3*RJK2+RIJ2*RJK3)/(CIJKL*SIJK*SJKL*RIJ*RJK*RKL)- &
       SIJKL*RKL1/(CIJKL*SJKL**2*RKL**2) 
  BM(IROW,3*L-1)=-CJKL*SIJKL*RJK2/(CIJKL*SJKL**2*RJK*RKL)+  &
       (RIJ3*RJK1-RIJ1*RJK3)/(CIJKL*SIJK*SJKL*RIJ*RJK*RKL)- &
       SIJKL*RKL2/(CIJKL*SJKL**2*RKL**2)
  BM(IROW,3*L)=(-RIJ2*RJK1+RIJ1*RJK2)/  &
       (CIJKL*SIJK*SJKL*RIJ*RJK*RKL)-   &
       CJKL*SIJKL*RJK3/(CIJKL*SJKL**2*RJK*RKL)-  &
       SIJKL*RKL3/(CIJKL*SJKL**2*RKL**2)
enddo    

!SHLL 1509
! out-of-plane bend
!
do IJKL=1, NOB
   IROW=IROW+1
   I=IOB(1,IJKL)
   J=IOB(2,IJKL)
   K=IOB(3,IJKL)
   L=IOB(4,IJKL)
   e41(1)=X(I)-X(L)
   e41(2)=Y(I)-Y(L)
   e41(3)=Z(I)-Z(L)
   e42(1)=X(j)-X(L)
   e42(2)=Y(j)-Y(L)
   e42(3)=Z(j)-Z(L)
   e43(1)=X(k)-X(L)
   e43(2)=Y(k)-Y(L)
   e43(3)=Z(k)-Z(L)
   R41=DIST(I,L,X,Y,Z,NATOM)
   R42=DIST(J,L,X,Y,Z,NATOM)
   R43=DIST(K,L,X,Y,Z,NATOM)
   e41=e41/R41
   e42=e42/R42
   e43=e43/r43
   a243=angl(j,l,k,x,y,z,natom)
   c243=cose(a243)
   s243=sine(a243)
   call xprod(e42,e43,tmpv)
   oopb=asin(dot_product(tmpv,e41)/s243)
   coopb=cose(oopb)
   toopb=tan(oopb)
   
   call xprod(e42,e43,tmpv)
   bmv1=(tmpv/coopb/s243-toopb*e41)/r41
   bm(irow,3*I-2)=bmv1(1)
   bm(irow,3*i-1)=bmv1(2)
   bm(irow,3*i)=bmv1(3)
   call xprod(e43,e41,tmpv)
   bmv2=(tmpv/coopb/s243-toopb/s243**2*(e42-e43*c243))/r42
   bm(irow,3*j-2)=bmv2(1)
   bm(irow,3*j-1)=bmv2(2)
   bm(irow,3*j)=bmv2(3)
   call xprod(e41,e42,tmpv)
   bmv3=(tmpv/coopb/s243-toopb/s243**2*(e43-e42*c243))/r43
   bm(irow,3*k-2)=bmv3(1)
   bm(irow,3*k-1)=bmv3(2)
   bm(irow,3*k)=bmv3(3)
   bmv4=-bmv1-bmv2-bmv3
   bm(irow,3*l-2)=bmv4(1)
   bm(irow,3*l-1)=bmv4(2)
   bm(irow,3*l)=bmv4(3)
end do
!SHLL 1510
! out-of-plane distance (aka pyramid height)
!
do ijkl=1, nod
   irow=irow+1
   I=iod(1,ijkl)
   J=iod(2,ijkl)
   K=iod(3,ijkl)
   L=iod(4,ijkl)
   e12(1)=X(j)-X(i)
   e12(2)=y(j)-y(i)
   e12(3)=z(j)-z(i)
   e13(1)=x(k)-x(i)
   e13(2)=y(k)-y(i)
   e13(3)=z(k)-z(i)
   v14(1)=x(l)-x(i)
   v14(2)=y(l)-y(i)
   v14(3)=z(l)-z(i)
   r12=dist(i,j,x,y,z,natom)
   r13=dist(i,k,x,y,z,natom)
   e12(:)=e12(:)/r12
   e13(:)=e13(:)/r13
   theta=angl(j, i, k,x,y,z,natom)
   call xprod(e12, e13, tmpv2)
   rtmp=dot_product(v14, tmpv2)

   call xprod(e13, v14, tmpv)
   tmpv3(:)=(cose(theta)*e12(:)-e13(:))/r12/sine(theta)
   bmv2 = (tmpv-rtmp*e12)/sine(theta)/r12 &
        - cose(theta)*rtmp/sine(theta)**2*tmpv3

   call xprod(v14, e12, tmpv)
   tmpv3(:)=(cose(theta)*e13(:)-e12(:))/r13/sine(theta)
   bmv3 = (tmpv-rtmp*e13)/sine(theta)/r13 &
        - cose(theta)*rtmp/sine(theta)**2*tmpv3

   call xprod(e12, e13, tmpv)
   bmv4 = tmpv/sine(theta)

   bmv1 = -bmv2-bmv3-bmv4

   if(rtmp/sine(theta)<0) then
      bmv1 = -bmv1
      bmv2 = -bmv2
      bmv3 = -bmv3
      bmv4 = -bmv4
   end if

   bm(irow,3*I-2)=bmv1(1)
   bm(irow,3*i-1)=bmv1(2)
   bm(irow,3*i)=bmv1(3)
   bm(irow,3*j-2)=bmv2(1)
   bm(irow,3*j-1)=bmv2(2)
   bm(irow,3*j)=bmv2(3)
   bm(irow,3*k-2)=bmv3(1)
   bm(irow,3*k-1)=bmv3(2)
   bm(irow,3*k)=bmv3(3)
   bm(irow,3*l-2)=bmv4(1)
   bm(irow,3*l-1)=bmv4(2)
   bm(irow,3*l)=bmv4(3)
end do
!SHLL end
END

!***********************************************************************
!     DIST
!***********************************************************************
!
DOUBLE PRECISION FUNCTION DIST(I,J,X,Y,Z,NATOM)
!
!     CALCULATES THE DISTANCE I-J
!     CALLED BY: BMAT
implicit none
integer :: i, j, natom
double precision :: X(NATOM),Y(NATOM),Z(NATOM)
DIST=SQRT((X(I)-X(J))**2+(Y(I)-Y(J))**2+(Z(I)-Z(J))**2)
END
!
!***********************************************************************
!     ANGL
!***********************************************************************
!
DOUBLE PRECISION FUNCTION ANGL(I,J,K,X,Y,Z,NATOM)
!
!     CALCULATES THE ANGLE I-J-K
!     CALLED BY:
!              BANGLE1,BMAT,BTENS,BTORSN,CONOUT,READINT,PTORS
!
implicit none
integer :: i, j, k, natom
double precision :: X(NATOM),Y(NATOM),Z(NATOM),dist
!
ANGL=((X(J)-X(I))*(X(J)-X(K))+(Y(J)-Y(I))*(Y(J)-Y(K))+  &
       (Z(J)-Z(I))*(Z(J)-Z(K)))/(DIST(I,J,X,Y,Z,NATOM)*DIST(J,K,X,Y,Z,NATOM))
IF (ABS(ANGL).GT.1d0) ANGL = 1.0d0
ANGL=ACOS(ANGL)
END
!
!***********************************************************************
!     SINE
!***********************************************************************
!
DOUBLE PRECISION FUNCTION sine(PHI)
!
!     CALLED BY: BMAT
!
implicit none
double precision, PARAMETER :: TINY=1.D-20 
double precision :: phi
SINE=SIN(PHI)
IF(ABS(SINE).LT.TINY) SINE=SIGN(TINY,SINE)
END
!
!***********************************************************************
!     COSE
!***********************************************************************
!
DOUBLE PRECISION FUNCTION cose(PHI)
!
!     CALLED BY: BMAT
!
implicit none
double precision, PARAMETER ::  TINY=1.D-20
double precision :: phi
COSE=COS(PHI)
IF(ABS(COSE).LT.TINY) COSE=SIGN(TINY,COSE)
END

!
!***********************************************************************
!     XPROD
!***********************************************************************
!
SUBROUTINE xprod(X,Y,Z)
!
!  FORMS THE CROSS PRODUCT OF VECTORS X AND Y, PLACING THE RESULT
!  IN VECTOR Z.  ALL VECTORS ARE THREE DIMENSIONAL.
!     CALLED BY: BMAT
!
implicit none
double precision :: X(3),Y(3),Z(3)
Z(1) = X(2)*Y(3) - X(3)*Y(2)
Z(2) = X(3)*Y(1) - X(1)*Y(3)
Z(3) = X(1)*Y(2) - X(2)*Y(1)
END
