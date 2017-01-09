C     Francisco J. Rodriguez-Cortes, November 2016
C
C     This code provides a non-parametric kernel based estimator of the
C     spatial variance mark function.
C

       subroutine Vmrcore(x,y,txy,n,s,ns,ks,hs,emr,Vmr)

       implicit real*8(a-h,o-z)

       integer i,j,iu,n,ns,ks
       double precision wij,vij,hs,kerns,Vrm,Vrn,Vmr,x,y,txy
       double precision hij,mij,xi,yi,ti,two,emr
       dimension x(n),y(n),txy(n),s(ns),Vrm(ns),Vrn(ns)
       dimension Vmr(ns),ks(3),emr(ns)
       
       Vrm=0d0
       Vrn=0d0

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
                    wij=mij*kerns
                    vij=kerns
                    Vrm(iu)=Vrm(iu)+wij
                    Vrn(iu)=Vrn(iu)+vij
             end if
           end if
          end do
          end do
            Vmr(iu)=Vrm(iu)/Vrn(iu)
          end do

        return

        end
