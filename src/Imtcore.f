C     Francisco J. Rodriguez-Cortes, November 2016
C
C     This code provides a non-parametric kernel based estimator of the
C     temporal Moran's I-statistics.
C

       subroutine Imtcore(snorm,txy,n,t,nt,kt,ht,mumt,Imt)

       implicit real*8(a-h,o-z)

       integer i,j,iv,n,nt,kt
       double precision wij,vij,ht,kernt,Itm,Itn,Imt,snorm,txy
       double precision tij,mij,snormi,ti,mumt
       dimension snorm(n),txy(n),t(nt),Itm(nt),Itn(nt),Imt(nt),kt(3)

       Itm=0d0
       Itn=0d0

       do iv=1,nt
        do i=1,n
         snormi=snorm(i)
         ti=txy(i)
          do j=1,n
           if (j.ne.i) then
            tij=abs(ti-txy(j))
            mij=(snormi-mumt)*(snorm(j)-mumt)
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
                    Itm(iv)=Itm(iv)+wij
                    Itn(iv)=Itn(iv)+vij
             end if
           end if
          end do
          end do
            Imt(iv)=Itm(iv)/Itn(iv)
          end do

        return

        end
