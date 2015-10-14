SUBROUTINE studnt (t, doff, fn_val)

! N.B. Argument IFAULT has been removed.
 
! Code converted using TO_F90 by Alan Miller
! Date: 2002-01-02  Time: 21:19:45

!     ALGORITHM AS 27  APPL. STATIST. VOL.19, NO.1

!     Calculate the upper tail area under Student's t-distribution

!     Translated from Algol by Alan Miller

IMPLICIT NONE
double precision t
double precision doff
double precision fn_val

!     Local variables

REAL     :: v, x, tt
LOGICAL  :: pos
REAL, PARAMETER  :: four = 4.0, one = 1.0, zero = 0.0, half = 0.5
REAL, PARAMETER  :: a1 = 0.09979441, a2 = -0.581821, a3 = 1.390993,  &
                    a4 = -1.222452, a5 = 2.151185
REAL, PARAMETER  :: b1 = 5.537409, b2 = 11.42343
REAL, PARAMETER  :: c1 = 0.04431742, c2 = -0.2206018, c3 = -0.03317253,  &
                    c4 = 5.679969, c5 = -12.96519
REAL, PARAMETER  :: d1 = 5.166733, d2 = 13.49862
REAL, PARAMETER  :: e1 = 0.009694901, e2 = -0.1408854, e3 = 1.88993,  &
                    e4 = -12.75532, e5 = 25.77532
REAL, PARAMETER  :: f1 = 4.233736, f2 = 14.3963
REAL, PARAMETER  :: g1 = -9.187228E-5, g2 = 0.03789901, g3 = -1.280346,  &
                    g4 = 9.249528, g5 = -19.08115
REAL, PARAMETER  :: h1 = 2.777816, h2 = 16.46132
REAL, PARAMETER  :: i1 = 5.79602E-4, i2 = -0.02763334, i3 = 0.4517029,  &
                    i4 = -2.657697, i5 = 5.127212
REAL, PARAMETER  :: j1 = 0.5657187, j2 = 21.83269

!     Check that number of degrees of freedom > 4.

!IF (doff <= four) THEN
  !WRITE(*, *) '** Error in AS27 - degrees of freedom <= 4  **'
!        fn_val = -0.5
!  RETURN
!END IF

!     Evaluate series.


v = one / doff
pos = (t >= zero)
tt = ABS(t)
x = half*(one +   &
    tt*(((a1 + v*(a2 + v*(a3 + v*(a4 + v*a5)))) / (one - v*(b1 - v*b2))) +  &
    tt*(((c1 + v*(c2 + v*(c3 + v*(c4 + v*c5)))) / (one - v*(d1 - v*d2))) +  &
    tt*(((e1 + v*(e2 + v*(e3 + v*(e4 + v*e5)))) / (one - v*(f1 - v*f2))) +  &
    tt*(((g1 + v*(g2 + v*(g3 + v*(g4 + v*g5)))) / (one - v*(h1 - v*h2))) +  &
    tt*((i1 + v*(i2 + v*(i3 + v*(i4 + v*i5)))) / (one - v*(j1 - v*j2))) ))))) ** (-8)
IF (pos) THEN
  fn_val = x
ELSE
  !fn_val = one - x
        fn_val = x
END IF

RETURN
END SUBROUTINE studnt