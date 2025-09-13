C     Francisco J. Rodriguez-Cortes, November 2016
C
C     This code provides a non-parametric kernel based estimator of the
C     temporal t-mark function.
C
subroutine kmtcore(xy, txy, n, t, nt, kt, ht, nxw, kmt)
  implicit none

  !-- argumentos
  integer, intent(in) :: n        ! número de puntos (filas de xy)
  integer, intent(in) :: nt       ! número de puntos en vector t
  integer, intent(in) :: kt(3)    ! indicador de kernel (kt(1) o kt(2) o kt(3) == 1)
  double precision, intent(in) :: xy(n,2)   ! matriz n x 2
  double precision, intent(in) :: txy(n)    ! tiempos asociados a cada fila de xy
  double precision, intent(in) :: t(nt)     ! vector de tiempos donde evaluar
  double precision, intent(in) :: ht        ! ancho de banda temporal
  double precision, intent(in) :: nxw       ! factor de normalización (no debe ser 0)
  double precision, intent(out) :: kmt(nt)  ! resultado: kmt para cada t(iv)

  !-- variables locales
  integer :: iv, i, j
  double precision :: ktm(nt), ktn(nt)
  double precision :: wij, vij, kernt, tij, mij
  double precision :: xyi(2), ti

  !-- inicializar
  ktm = 0d0
  ktn = 0d0
  kmt = 0d0

  !-- bucle sobre tiempos de evaluación
  do iv = 1, nt
     ! para cada punto i
     do i = 1, n
        xyi = xy(i, :)     ! subarray: fila i (dos componentes)
        ti = txy(i)
        do j = 1, n
           if (j /= i) then
              tij = dabs(ti - txy(j))
              ! producto escalar entre vectores de dimensión 2
              mij = dot_product(xyi, xy(j, :)) / nxw

              ! seleccionar kernel temporal (asegúrate de que las funciones existen)
              if (kt(1) == 1) then
                 kernt = boxkernel((t(iv) - tij)/ht, ht)
              else if (kt(2) == 1) then
                 kernt = ekernel((t(iv) - tij)/ht, ht)
              else if (kt(3) == 1) then
                 kernt = qkernel((t(iv) - tij)/ht, ht)
              else
                 kernt = 0d0
              end if

              if (kernt /= 0d0) then
                 wij = mij * kernt
                 vij = kernt
                 ktm(iv) = ktm(iv) + wij
                 ktn(iv) = ktn(iv) + vij
              end if
           end if
        end do
     end do

     ! evitar división por cero
     if (ktn(iv) == 0d0) then
        kmt(iv) = 0d0
     else
        kmt(iv) = ktm(iv) / ktn(iv)
     end if
  end do

  return
end subroutine kmtcore
