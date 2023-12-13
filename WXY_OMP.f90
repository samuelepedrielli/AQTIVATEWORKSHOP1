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
  ipsip=ivic(isite,j)!go through the neighbors of the selected site
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
integer(kind=ik) :: i,nn,k
real(kind=rk) :: E,E2
real(kind=rk) :: en
real(kind=rk), dimension(:) :: iconf
integer(kind=ik), dimension(:,:) :: ivic
E=0.0_rk
E2=0.0_rk
do i=1,N
 do nn=1,2
  k=ivic(i,nn)
  call ham(iconf(i),iconf(k),en)
  E=E+en
 enddo
enddo
E=E/N
E2=E**2
end subroutine computeE

subroutine computeM(iconf,M,M2)
implicit none
integer(kind=ik) :: i
real(kind=rk) :: M,M2,Mx,My
real(kind=rk), dimension(:) :: iconf
Mx=0.0_rk
My=0.0_rk
M=0.0_rk
M2=0.0_rk
do i=1,size(iconf)
 My=My+sin(iconf(i))
 Mx=Mx+cos(iconf(i))
enddo
M=sqrt(Mx**2+My**2)/size(iconf)
M2=M**2
end subroutine computeM

end module xymodule

program xy
use xymodule
use omp_lib
implicit none
integer :: sizer
integer(kind=ik) :: Nrep,Niter,Neq,Nskip,Ninit,Ntemp,icount,t1,t2,rate
integer(kind=ik) :: i,ii,j,jj,k,kk,isite,ipsip,nn
integer(kind=ik) :: dimclu,dimclu2
real(kind=rk) r,t_i,t_f
real(kind=rk) acount,acc_rate
real(kind=rk) :: E,E2,sommaE,sommaE2,Emedio,E2medio,varE,errE,somma1E,somma1E2,en,Mx,My
real(kind=rk) :: M,M2,M4,sommaM,sommaM2,sommaM4,Mmedio,M2medio,M4medio,varM,errM,somma1M,somma1M2,somma1M4
real(kind=rk) :: S,sommaS,S2,sommaS2,varS,errS,somma1S,somma1S2
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
real(kind=rk), dimension(:), allocatable :: cluster,iconf,initial
character (len = 15) :: conf,phys,err,input
logical :: exist

!--------------------------------------------!
!read(*,*) input
!open(40, file=input)
!read(40,*) temp
!print*,temp

L=16
temp=.0
N=L**2
Neq=1e4
Niter=Neq
Nskip=100
Nrep=10
phys="phys.txt"
err="err.txt"

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

!--------------------------------------------!
call system_clock(count_rate=rate)
call random_seed(sizer)
allocate(seme(sizer))

allocate (iconf(N))
allocate (ivic(N,4))
allocate (cid(N))
allocate (cluster(N))


call neighbors(ivic)

!print*,'------ XY MODEL ------'
!print*,'------ Temp = ', temp,'------'

Emedio=0.0_rk !energy
E2medio=0.0_rk
Mmedio=0.0_rk !magnetization
M2medio=0.0_rk
dimclumedio=0.0_rk!cluster dimension
errEmedio=0.0_rk
errMmedio=0.0_rk

call system_clock(t1)

!$OMP PARALLEL DO PRIVATE(jj,iconf,cid,cluster,seme,icount,sommaE,sommaM,E,en)
do jj=1,Nrep !begin statistical loop
! print*,'------ Stat = ', jj,'------'

 do i = 1,sizer
  seme(i) = i+omp_get_thread_num()+1
 enddo
 call random_seed(get=seme)
 call random_seed(put=seme)

 do k=1,N !initialize the disordered configuration
  call random_number(r)!pick a random angle
  r=r*2*pi
  iconf(k)=r
 enddo
 

 sommaE=0.0_rk !energy
 sommaM=0.0_rk !magnetization
 icount=0

 do j=1,Neq !equilibration loop
  call wolffmove(iconf,ivic,temp,cid,cluster,dimclu)
 enddo

 do ii=1,Niter !steps loop
  E=0.0_rk
  Mx=0.0_rk
  My=0.0_rk
  M=0.0_rk
  !if (mod(ii,50)==1_ik) print*,ii
  do j=1,Nskip !decorrelation loop
   call wolffmove(iconf,ivic,temp,cid,cluster,dimclu)
  enddo
  call wolffmove(iconf,ivic,temp,cid,cluster,dimclu)
  dimclumedio=dimclumedio+dimclu
  icount=icount+1
    ! COMPUTE ENERGY AND MAGNETIZATION
   do isite=1,N
    My=My+sin(iconf(isite))
    Mx=Mx+cos(iconf(isite))
    M=sqrt(Mx**2+My**2)/N
    sommaM=sommaM+M
    do nn=1,2
     en=0.0_rk
     ipsip=ivic(isite,nn)
     call ham(iconf(isite),iconf(ipsip),en)
     E=E+en
    enddo
   enddo
   E=E/N
   sommaE=sommaE+E
 enddo !ends steps loop

 sommaE=sommaE/icount !this is the <E> for one simulation
 sommaM=sommaM/icount !this is the <M> for one simulation

 Emedio = Emedio + sommaE
 E2medio = E2medio + sommaE**2
 Mmedio = Mmedio + sommaM
 M2medio = M2medio + sommaM**2
print*, jj, sommaE

enddo !ends statistical loop
!$OMP END PARALLEL DO

call system_clock(t2)

errEmedio=sqrt(E2medio/Nrep-(Emedio/Nrep)**2)
errMmedio=sqrt(M2medio/Nrep-(Mmedio/Nrep)**2)

Emedio=Emedio/Nrep
Mmedio=Mmedio/Nrep

dimclumedio=dimclumedio/(Nrep*Niter)

write(21,*) temp,Emedio,Mmedio !,Cvmedio,Chimedio,Gmedio,Bmedio
write(22,*) temp,errEmedio,errMmedio !,errCvmedio,errChimedio,errGmedio,errBmedio

write(6,*) 'Number of sites    =',N
write(6,*) 'Temperature        =',temp
write(6,*) 'Acceptance rate    =',acc_rate
write(6,*) 'Cluster dimension  =',dimclumedio
write(6,*)
write(6,*) 'Mean energy        =',Emedio,errEmedio
write(6,*) 'Magnetization      =',Mmedio,errMmedio
write(6,*) 'Execution time [s] =',real(t2-t1)/real(rate)

close(21)
close(22)

deallocate(seme)
deallocate(iconf)
deallocate(ivic)
deallocate(cid)
deallocate(cluster)

end program xy
