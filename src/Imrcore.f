C     Francisco J. Rodriguez-Cortes, November 2016
C
C     This code provides a non-parametric kernel based estimator of the
C     spatial Moran's I-statistics.
C

       subroutine Imrcore(x,y,txy,n,s,ns,ks,hs,mumr,Imr)

       implicit real*8(a-h,o-z)

       integer i,j,iu,n,ns,ks
       double precision wij,vij,hs,kerns,Irm,Irn,Imr,x,y,txy
       double precision hij,mij,xi,yi,ti,two,mumr
       dimension x(n),y(n),txy(n),s(ns),Irm(ns),Irn(ns),Imr(ns),ks(3)

       Irm=0d0
       Irn=0d0

          two=2d0

       do iu=1,ns
        do i=1,n
         xi=x(i)
         yi=y(i)
         ti=txy(i)
          do j=1,n
           if (j.ne.i) then
            hij=sqrt(((xi-x(j))**two)+((yi-y(j))**two))
            mij=(ti-mumr)*(txy(j)-mumr)
              if (ks(1).eq.1) then
               kerns=boxkernel((s(iu)-hij)/hs,hs)
                else if (ks(2).eq.1) then
                 kerns=ekernel((s(iu)-hij)/hs,hs)
                  else if (ks(3).eq.1) then
                   kerns=qkernel((s(iu)-hij)/hs,hs)
              end if
             if (kerns.ne.0d0) then
                    wij=mij*kerns
                    vij=kerns
                    Irm(iu)=Irm(iu)+wij
                    Irn(iu)=Irn(iu)+vij
             end if
           end if
          end do
          end do
            Imr(iu)=Irm(iu)/Irn(iu)
          end do

        return

        end
