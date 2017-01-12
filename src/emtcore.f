C     Francisco J. Rodriguez-Cortes, November 2016
C
C     This code provides a non-parametric kernel based estimator of the
C     temporal mean function.
C

       subroutine emtcore(x,y,txy,n,t,nt,kt,ht,emt)

       implicit real*8(a-h,o-z)

       integer i,j,iv,n,nt,kt
       double precision wij,vij,ht,kernt,etm,etn,emt,x,y,txy
       double precision tij,mij,xi,yi,ti
       dimension x(n),y(n),txy(n),t(nt),etm(nt),etn(nt),emt(nt),kt(3)

       etm=0d0
       etn=0d0

          two=2d0

       do iv=1,nt
        do i=1,n
         xi=x(i)
         yi=y(i)
         ti=txy(i)
          do j=1,n
           if (j.ne.i) then
            tij=abs(ti-txy(j))
            mij=sqrt(((xi-x(j))**two)+((yi-y(j))**two))
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
                    etm(iv)=etm(iv)+wij
                    etn(iv)=etn(iv)+vij
             end if
           end if
          end do
          end do
            emt(iv)=etm(iv)/etn(iv)
          end do

        return

        end
