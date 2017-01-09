C     Francisco J. Rodriguez-Cortes, November 2016
C
C     This code provides a non-parametric kernel based estimator of the
C     temporal variance mark function.
C
       subroutine Vmtcore(snorm,txy,n,t,nt,kt,ht,emt,Vmt)

       implicit real*8(a-h,o-z)

       integer i,j,iv,n,nt,kt
       double precision wij,vij,ht,kernt,Vtm,Vtn,Vmt,snorm,txy
       double precision tij,mij,snormi,ti,emt,two
       dimension snorm(n),txy(n),t(nt),Vtm(nt),Vtn(nt)
       dimension Vmt(nt),kt(3),emt(nt)
       

       Vtm=0d0
       Vtn=0d0

          two=2d0
          
       do iv=1,nt
        do i=1,n
         snormi=snorm(i)
         ti=txy(i)
          do j=1,n
           if (j.ne.i) then
            tij=abs(ti-txy(j))
            mij=(snormi-emt(iv))**two
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
                    Vtm(iv)=Vtm(iv)+wij
                    Vtn(iv)=Vtn(iv)+vij
             end if
           end if
          end do
          end do
            Vmt(iv)=Vtm(iv)/Vtn(iv)
          end do

        return

        end
