PROGRAM Diferential_Equations
    IMPLICIT NONE

    REAL(4), DIMENSION(100) :: xlist, tlist
    REAL(4) :: t0, x0, dt, dx_dt
    INTEGER :: N, i, ios
    REAL(4), EXTERNAL :: f

    t0 = 0.0
    x0 = 1.0
    dt = 0.1
    N = 100


    CALL Euler_method(t0, x0, dt, N, tlist, xlist)

   

     OPEN(UNIT=100, IOSTAT=ios, FILE='euler.txt', STATUS='new', ACTION='write')

     IF (ios /= 0) THEN
       PRINT*, 'ERROR: Unable to open file'
       STOP
     ELSE
       WRITE(100, *) 't ', '                   ', 'x '
       DO i = 1, N 
        WRITE(100, *) tlist(i), xlist(i)                          
       END DO                       
     CLOSE(100)
     END IF




    CALL Midpoint(t0, x0, dt, N, tlist, xlist)


     OPEN(UNIT=101, IOSTAT=ios, FILE='Midpoint.txt', STATUS='replace', ACTION='write')

     IF (ios /= 0) THEN
       PRINT*, 'ERROR: Unable to open file'
       STOP
     ELSE
      WRITE(101, *) 't ', '                   ', 'x '
       DO i = 1, N 
        WRITE(101, *) tlist(i), xlist(i)                             
       END DO                       
     CLOSE(101)
     END IF

   

END PROGRAM Diferential_Equations

REAL(4) FUNCTION f(x) RESULT (dx_dt)
    REAL(4) :: x
    dx_dt = -x
END FUNCTION f

SUBROUTINE Euler_method(t0, x0, dt, N, tlist, xlist)
    IMPLICIT NONE
    REAL(4), INTENT(IN) :: t0, x0, dt
    REAL(4), DIMENSION(100), INTENT(OUT) :: tlist, xlist
    INTEGER(4) :: N, i
    REAL :: f, fx, t, x

    x = x0
    t = t0
    tlist(1) = t
    xlist(1) = x

    DO i = 1, N
      fx = f(x)
      x = x + fx * dt
      t = t + dt
      tlist(i + 1) = t
      xlist(i + 1) = x
    END DO
  
 
END SUBROUTINE Euler_method


SUBROUTINE Midpoint(t0, x0, dt, N, tlist, xlist)
    IMPLICIT NONE
    REAL(4), INTENT(IN) :: t0, x0, dt
    REAL(4), DIMENSION(100), INTENT(OUT) :: tlist, xlist
    INTEGER(4) :: N, i
    REAL :: f, fx, t, x

    x = x0
    t = t0
    tlist(1) = t
    xlist(1) = x

    DO i = 1, N
      fx = f(x)
      x = x + fx * (dt / 2.0)
      t = t + dt
      tlist(i + 1) = t
      xlist(i + 1) = x
    END DO
  
  
END SUBROUTINE Midpoint