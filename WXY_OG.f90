module xymodule
implicit none
integer, parameter :: rk=selected_real_kind(8)
integer, parameter :: ik=selected_int_kind(8)
real(kind=rk), parameter :: pi=4*atan(1.0)
integer :: L,N
contains

subroutine neighbors(ivic)
implicit none
integer(kind=ik) :: i,ir,il,iu,id
integer(kind=ik), intent(out), dimension(:,:) :: ivic

do i=1,N
 ir=i+1_ik
 if (mod(i,L)==0) ir=i+1_ik-L
 iu=i+L
 if (i>L*(L-1_ik)) iu=i+L-L**2
 il=i-1_ik
 if (mod(i,L)==1_ik) il=i-1_ik+L
 id=i-L
 if (i<=L) id=i-L+L**2
 ivic(i,1)=ir
 ivic(i,2)=iu
 ivic(i,3)=il
 ivic(i,4)=id
enddo
end subroutine neighbors

subroutine ham(v,u,en)
real(kind=rk) :: v,u,en
if ( (0 .le. v) .and. (v .le. 2*pi) .and. (0 .le. u) .and. (u .le. 2*pi) ) then
   en=-cos(v-u)
else
print*,"THESE ARE NOT ACCEPTABLE STATES: ",v,u
endif
end subroutine ham

subroutine wolffmove(iconf,ivic,temp,cid,cluster,dimclu)
implicit none
integer(kind=ik) :: dimclu
integer(kind=ik) :: low,high,isite
integer(kind=ik) :: i,j,k,ipsip
real(kind=rk) :: r,phi,test1,test2,temp,en1,en2,dot1,dot2
integer(kind=ik), dimension(:,:) :: ivic
integer(kind=ik), dimension(:) :: cid
real(kind=rk), dimension(:) :: iconf
real(kind=rk), dimension(:) :: cluster !theese will have the dimension of the lattice (N)
real(kind=rk), dimension(2) :: rv

dimclu=0_ik !initial dimension of the cluster is zero
low=0_ik !counts the number of sites of the cluster that have been selected
do j=1,N !initialize the cluster to zero, idea: the cluster cannot be bigger than N so you initialize a N-array with zeros
cid(j)=0_ik
cluster(j)=0_ik
enddo

call random_number(r) !pick a random site
isite=N*r+1_ik
call random_number(phi) !pick a random angle
phi=phi*2*pi
iconf(isite)=abs(modulo(2*phi-iconf(isite)+pi,2*pi))!flip the chosen site
cid(isite)=1_ik!signal that you have checked the isite site
cluster(low+1_ik)=isite!add the site to the cluster in the first position
high=1_ik !counts the number of sites in the cluster

do while(low.lt.high) !while the lower bound is less than the higher bound
 isite=cluster(low+1_ik)!take the starting site from the cluster
 do j=1,4
  ipsip=ivic(isite,j)!go through the neighbors of the selected spine
  if(cid(ipsip)==0_ik) then !check that the neighbour has not been asked yet
   dot1=cos(phi)*cos(iconf(isite))+sin(phi)*sin(iconf(isite))
   dot2=cos(phi)*cos(iconf(ipsip))+sin(phi)*sin(iconf(ipsip))
   test1=1-exp(2.*dot1*dot2/temp)
   call random_number(r)
   test2=r
   if(test2<test1) then !if th MC condition is satisfied
    cid(ipsip)=1_ik !add it to the list of asked sites
    cluster(high+1_ik)=ipsip!add the site to the cluster
    high=high+1_ik!raise the higher bound
    iconf(ipsip)=abs(modulo(2*phi-iconf(ipsip)+pi,2*pi))!rotate the added site
   endif
  endif
 enddo
 low=low+1_ik!raise the lower bound
enddo

dimclu=high!the dimension of the cluster is the number of sites which belong to it

end subroutine wolffmove

subroutine computeE(iconf,ivic,E,E2)
implicit none
integer(kind=ik) :: i,jj,k
real(kind=rk) :: E,E2
real(kind=rk) :: en
real(kind=rk), dimension(:) :: iconf
integer(kind=ik), dimension(:,:) :: ivic
E=0.0_rk
E2=0.0_rk
do i=1,size(iconf)
 do jj=1,2
  k=ivic(i,jj)
  call ham(iconf(i),iconf(k),en)
  E=E+en
 enddo
enddo
E=E/size(iconf)
E2=E**2
end subroutine computeE

subroutine computeM(iconf,M,M2,M4)
implicit none
integer(kind=ik) :: i
real(kind=rk) :: M,M2,M4,Mx,My
real(kind=rk), dimension(:) :: iconf
Mx=0.0_rk
My=0.0_rk
M=0.0_rk
M2=0.0_rk
M4=0.0_rk
do i=1,size(iconf)
 My=My+sin(iconf(i))
 Mx=Mx+cos(iconf(i))
enddo
M=sqrt(Mx**2+My**2)/size(iconf)
M2=M**2
M4=M**4
end subroutine computeM

subroutine computeS(iconf,ivic,S,S2)
implicit none
integer(kind=ik) :: i,k
real(kind=rk) :: S,S2
real(kind=rk), dimension(:) :: iconf
integer(kind=ik), dimension(:,:) :: ivic
S=0.0_rk
S2=0.0_rk
do i=1,size(iconf)
  k=ivic(i,1)
  S=S+sin(iconf(i)-iconf(k))
enddo
S=S/size(iconf)
S=S**2 !S=s**2
S2=S**2
end subroutine computeS

end module xymodule

program xy
use xymodule
implicit none
integer :: sizer
integer(kind=ik) :: Nrep,Niter,Neq,Nskip,Ninit,Ntemp,icount
integer(kind=ik) :: i,ii,j,jj,k,kk
integer(kind=ik) :: dimclu,dimclu2
real(kind=rk) r,t_i,t_f
real(kind=rk) acount,acc_rate
real(kind=rk) :: E,E2,sommaE,sommaE2,Emedio,E2medio,varE,errE
real(kind=rk) :: M,M2,M4,sommaM,sommaM2,sommaM4,Mmedio,M2medio,M4medio,varM,errM
real(kind=rk) :: S,sommaS,S2,sommaS2,varS,errS
real(kind=rk) :: Cv,Cv2,Cvmedio,Cv2medio
real(kind=rk) :: Chi,Chi2,Chimedio,Chi2medio
real(kind=rk) :: G,G2,Gmedio,G2medio
real(kind=rk) :: B,B2,Bmedio,B2medio
real(kind=rk) :: errEmedio,errMmedio,errCvmedio,errChimedio,errGmedio,errBmedio
real(kind=rk) :: dimclumedio,dimclu2medio,vardimclu,errdimclu
real(kind=rk) :: temp
integer, dimension(:), allocatable :: seme
integer(kind=ik), dimension(:,:), allocatable :: ivic
integer(kind=ik), dimension(:), allocatable :: cid
real(kind=rk), dimension(:), allocatable :: cluster,iconf
character (len = 15) :: conf,phys,err,input
logical :: exist

call random_seed(sizer)
allocate(seme(sizer))
do i = 1,sizer
 seme(i) = i
enddo
call random_seed(get=seme)
call random_seed(put=seme)

!--------------------------------------------!

!read(*,*) input
!open(40, file=input)
!read(40,*) temp
!print*,temp

L = 8
temp = 0. !5.
N=L**2
Neq=10000
Niter = 1e4
Nskip = 100
!temp = 2.
Nrep = 1
conf = "conf.txt"
phys = "phys.txt"
err = "err.txt"

inquire(file=conf, exist=exist)
  if (exist) then
    open(20, file=conf, status="old", position="append", action="write")
  else
    open(20, file=conf, status="new", action="write")
  end if
inquire(file=phys, exist=exist)
  if (exist) then
    open(21, file=phys, status="old", position="append", action="write")
  else
    open(21, file=phys, status="new", action="write")
  end if
inquire(file=err, exist=exist)
  if (exist) then
    open(22, file=err, status="old", position="append", action="write")
  else
    open(22, file=err, status="new", action="write")
  end if

allocate (iconf(N))
allocate (ivic(N,4))
allocate (cid(N))
allocate (cluster(N))
!--------------------------------------------!

call neighbors(ivic)

!print*,'------ XY MODEL ------'
!print*,'------ Temp = ', temp,'------'

Emedio=0.0_rk !energy
E2medio=0.0_rk
Mmedio=0.0_rk !magnetization
M2medio=0.0_rk
M4medio=0.0_rk
Cvmedio=0.0_rk !capacity
Cv2medio=0.0_rk
Chimedio=0.0_rk !susceptibility
Chi2medio=0.0_rk
Gmedio=0.0_rk !helicity modulus
G2medio=0.0_rk
Bmedio=0.0_rk !binder
B2medio=0.0_rk
dimclumedio=0.0_rk!cluster dimension

errEmedio=0.0_rk
errMmedio=0.0_rk
errCvmedio=0.0_rk
errChimedio=0.0_rk
errGmedio=0.0_rk
errBmedio=0.0_rk

do jj=1,Nrep !begin statistical loop
! print*,'------ Stat = ', jj,'------'

do k=1,N !initialize the ordered configuration
 call random_number(r)!pick a random angle
 r=r*2*pi
 iconf(k)=r
enddo

 sommaE=0.0_rk !energy
 sommaE2=0.0_rk

 sommaM=0.0_rk !magnetization
 sommaM2=0.0_rk
 sommaM4=0.0_rk

 sommaS=0.0_rk !helicity
 sommaS2=0.0_rk

 varE=0.0_rk
 errE=0.0_rk
 varM=0.0_rk
 errM=0.0_rk

 Cv=0.0_rk
 Cv2=0.0_rk
 Chi=0.0_rk
 Chi2=0.0_rk
 G=0.0_rk
 G2=0.0_rk
 B=0.0_rk
 B2=0.0_rk

 icount=0_ik

 do j=1,Neq !equilibration loop
  call wolffmove(iconf,ivic,temp,cid,cluster,dimclu)
  !if (mod(j,100)==1) then
  ! write(20,*) iconf '\n'
  !end if
 enddo

 do ii=1,Niter !steps loop
 ! if (mod(ii,50)==1_ik) print*,ii 
  do j=1,Nskip !decorrelation loop
   call wolffmove(iconf,ivic,temp,cid,cluster,dimclu)
  enddo
  call wolffmove(iconf,ivic,temp,cid,cluster,dimclu)
  dimclumedio=dimclumedio+dimclu
  dimclu2medio=dimclu2medio+dimclu**2
  icount=icount+1_ik
  call computeE(iconf,ivic,E,E2)
  sommaE=sommaE+E
  sommaE2=sommaE2+E2
  call computeM(iconf,M,M2,M4)
  sommaM=sommaM+M
  sommaM2=sommaM2+M2
  sommaM4=sommaM4+M4
  call computeS(iconf,ivic,S,S2)
  sommaS=sommaS+S
  sommaS2=sommaS2+S2
 enddo !ends steps loop

 varE=sommaE2/icount-(sommaE/icount)**2
 errE=sqrt(varE/Nrep)

 Cv=N*varE/temp**2 !this is Cv for one configuration
 Cv2=Cv**2

 varM=sommaM2/icount-(sommaM/icount)**2
 errM=sqrt(varM/Nrep)


 Chi=N*varM/temp !this is Chi for one configuration
 Chi2=Chi**2

 sommaE=sommaE/icount !this is the <E> for one simulation
 sommaE2=sommaE2/icount

 sommaM=sommaM/icount !this is the <M> for one simulation
 sommaM2=sommaM2/icount
 sommaM4=sommaM4/icount
 
 !print*, jj, 'energy', sommaE, errE
 !print*, jj, 'mag', sommaM, errM

 sommaS=sommaS/icount !this is the <s2> for one simulation
 sommaS2=sommaS2/icount

 G=-(sommaE/2.+N*sommaS/temp) !this is G for one simulation
 G2=G**2

 B=1.-sommaM4/(3.*sommaM2**2) !this is B for one simulation
 B2=B**2

 !the following add the found quantities to the statistical average
 Emedio=Emedio+sommaE
 E2medio=E2medio+sommaE**2

 Mmedio=Mmedio+sommaM
 M2medio=M2medio+sommaM**2
 M4medio=M4medio+sommaM**4

 Cvmedio=Cvmedio+Cv
 Cv2medio=Cv2medio+Cv2

 Chimedio=Chimedio+Chi
 Chi2medio=Chi2medio+Chi2

 Gmedio=Gmedio+G
 G2medio=G2medio+G2

 Bmedio=Bmedio+B
 B2medio=B2medio+B2

enddo !ends statistical loop

errEmedio=sqrt(E2medio/Nrep-(Emedio/Nrep)**2)
errMmedio=sqrt(M2medio/Nrep-(Mmedio/Nrep)**2)
errCvmedio=sqrt(Cv2medio/Nrep-(Cvmedio/Nrep)**2)
errChimedio=sqrt(Chi2medio/Nrep-(Chimedio/Nrep)**2)
errGmedio=sqrt(G2medio/Nrep-(Gmedio/Nrep)**2)
errBmedio=sqrt(B2medio/Nrep-(Bmedio/Nrep)**2)

Emedio=Emedio/Nrep
Mmedio=Mmedio/Nrep
Cvmedio=Cvmedio/Nrep
Chimedio=Chimedio/Nrep
Gmedio=Gmedio/Nrep
Bmedio=Bmedio/Nrep

dimclumedio=dimclumedio/(Nrep*Niter)
dimclu2medio=dimclu2medio/(Nrep*Niter)
vardimclu=dimclu2medio-dimclumedio**2
errdimclu=sqrt(vardimclu/(Nrep*Niter*N))

write(20,*) iconf
write(21,*) temp,Emedio,Mmedio,Cvmedio,Chimedio,Gmedio,Bmedio
write(22,*) temp,errEmedio,errMmedio,errCvmedio,errChimedio,errGmedio,errBmedio

write(6,*) 'Number of sites    =',N
write(6,*) 'Temperature        =',temp
write(6,*) 'Acceptance rate    =',acc_rate
write(6,*) 'Cluster dimension  =',dimclumedio
write(6,*)
write(6,*) 'Mean energy        =',Emedio,errEmedio
write(6,*) 'Magnetization      =',Mmedio,errMmedio
write(6,*) 'Specific heat      =',Cvmedio,errCvmedio
write(6,*) 'Susceptibility     =',Chimedio,errChimedio
write(6,*) 'Helicity modulus   =',Gmedio,errGmedio
write(6,*) 'Binder of M        =',Bmedio,errBmedio

close(20)
close(21)
close(22)

deallocate(seme)
deallocate(iconf)
deallocate(ivic)
deallocate(cid)
deallocate(cluster)

end program xy
