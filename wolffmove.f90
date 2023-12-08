program wolffmove
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

end program wolffmove
