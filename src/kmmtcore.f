C     Francisco J. Rodriguez-Cortes, November 2016
C
C     This code provides a non-parametric kernel based estimator of the
C     temporal mark correlation function.
C

       subroutine kmmtcore(xy,txy,n,t,nt,kt,ht,nxw,kmmt)

       implicit real*8(a-h,o-z)

       integer i,j,iv,n,nt,kt
       double precision wij,vij,ht,kernt,ktm,ktn,kmmt,xy,txy,nxw
       double precision tij,mij,xyi1,xyi2,ti
       dimension xy(n,2),txy(n),t(nt),ktm(nt),ktn(nt),kmmt(nt),kt(3)

       ktm=0d0
       ktn=0d0

       do iv=1,nt
        do i=1,n
         xyi1 = xy(i,1)
         xyi2 = xy(i,2)
         ti=txy(i)
          do j=1,n
           if (j.ne.i) then
            tij=dabs(ti-txy(j))
            mij = (xyi1*xy(j,1) + xyi2*xy(j,2))/nxw
              if (kt(1).eq.1) then
               kernt=boxkernel((t(iv)-tij)/ht,ht)
                else if (kt(2).eq.1) then
                 kernt=ekernel((t(iv)-tij)/ht,ht)
                  else if (kt(3).eq.1) then
                   kernt=qkernel((t(iv)-tij)/ht,ht)
              end if
             if (kernt.ne.0d0) then
                    wij=mij*kernt
                    vij=kernt
                    ktm(iv)=ktm(iv)+wij
                    ktn(iv)=ktn(iv)+vij
             end if
           end if
          end do
          end do
            kmmt(iv)=ktm(iv)/ktn(iv)
          end do

        return

        end
