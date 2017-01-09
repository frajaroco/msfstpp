C     Francisco J. Rodriguez-Cortes, November 2016
C
C     This code provides an edge-corrected non-parametric kernel based
C     estimator of the standardized spatial variance mark function. 
C

      subroutine Vmrcoreinh(x,y,txy,n,s,ns,slambda,ks,hs,wrs,wts,
     +     wbi,wbimod,wss,edg,emr,Vmr)
     
      implicit real*8(a-h,o-z)

      integer i,j,iu,n,ns,ks,edg
      double precision inhwij,inhvij,hs,kerns,Vmrminh,Vmrninh,Vmr
      double precision hij,mij,xi,yi,ti,two,wrs,wts,wbi,x,y,txy
      double precision wbimod,wss,slambda,emr
      dimension x(n),y(n),txy(n),s(ns),Vmrminh(ns),Vmrninh(ns)
      dimension wrs(n,n),wts(n,n),wbi(n,ns),wbimod(n,ns),wss(ns)
      dimension ks(3),edg(6),slambda(n),emr(ns),Vmr(ns)
       
       Vmrminh=0d0
       Vmrninh=0d0
      
          two=2d0

      do iu=1,ns
      do i=1,n
         xi=x(i)
         yi=y(i)
         ti=txy(i)      
      do j=1,n
      if (j.ne.i) then
      hij=sqrt(((xi-x(j))**two)+((yi-y(j))**two))
      mij=(ti-emr(iu))**two
      if (ks(1).eq.1) then
      kerns=boxkernel((s(iu)-hij)/hs,hs)
      else if (ks(2).eq.1) then
      kerns=ekernel((s(iu)-hij)/hs,hs)
      else if (ks(3).eq.1) then
      kerns=qkernel((s(iu)-hij)/hs,hs)
      end if
      if (kerns.ne.0d0) then
C    none
      if (edg(1).eq.1) then
       inhwij=(mij*kerns)/(slambda(i)*slambda(j))
       inhvij=kerns/(slambda(i)*slambda(j))
       Vmrminh(iu)=Vmrminh(iu)+inhwij
       Vmrninh(iu)=Vmrninh(iu)+inhvij  
      end if                  
C    isotropic
      if (edg(2).eq.1) then                  
      inhwij=(mij*kerns*wrs(i,j))/(slambda(i)*slambda(j))
      inhvij=(kerns*wrs(i,j))/(slambda(i)*slambda(j))
      Vmrminh(iu)=Vmrminh(iu)+inhwij
      Vmrninh(iu)=Vmrninh(iu)+inhvij
      end if
C    border
      if (edg(3).eq.1) then                  
      inhwij=(mij*kerns*wbi(i,iu))/(slambda(i)*slambda(j))
      inhvij=(kerns*wbi(i,iu))/(slambda(i)*slambda(j))
      Vmrminh(iu)=Vmrminh(iu)+inhwij
      Vmrninh(iu)=Vmrninh(iu)+inhvij
      end if
C    modified.border
      if (edg(4).eq.1) then
      inhwij=(mij*kerns*wbimod(i,iu))/(slambda(i)*slambda(j))
      inhvij=(kerns*wbimod(i,iu))/(slambda(i)*slambda(j))
      Vmrminh(iu)=Vmrminh(iu)+inhwij
      Vmrninh(iu)=Vmrninh(iu)+inhvij 
      end if                  
C    translate
      if (edg(5).eq.1) then
      inhwij=(mij*kerns*wts(i,j))/(slambda(i)*slambda(j))
      inhvij=(kerns*wts(i,j))/(slambda(i)*slambda(j))
      Vmrminh(iu)=Vmrminh(iu)+inhwij
      Vmrninh(iu)=Vmrninh(iu)+inhvij
      end if
C    setcovf         
      if (edg(6).eq.1) then
      inhwij=(mij*kerns*wss(iu))/(slambda(i)*slambda(j))
      inhvij=(kerns*wss(iu))/(slambda(i)*slambda(j))
      Vmrminh(iu)=Vmrminh(iu)+inhwij
      Vmrninh(iu)=Vmrninh(iu)+inhvij
      end if
      end if
      end if
       end do
       end do
       Vmr(iu)=Vmrminh(iu)/Vmrninh(iu)
       end do
      
        return
        
        end  
