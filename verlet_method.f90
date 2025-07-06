PROGRAM Diferential_EquationsII
    IMPLICIT NONE

    REAL(4), DIMENSION(300) :: tlist, x1list, y1list, x2list, y2list
    INTEGER :: N, i, ios, tf
    REAL(4), EXTERNAL :: f, v
    REAL, DIMENSION(2) :: x0, y0
    REAL :: t, dt, t0

    t0 = 0.0
    tf = 30
    x0 = (/2.0, 0.0/)
    dt = 0.1
    y0 = (/0.0, 0.5/)
    N = (tf - t0) / dt


    CALL Euler_method(y0, t0, x0, dt, N, tlist, x1list, x2list, y1list, y2list)

   

     OPEN(UNIT=100, IOSTAT=ios, FILE='EulerBonus.txt', STATUS='new', ACTION='write')

     IF (ios /= 0) THEN
       PRINT*, 'ERROR: Unable to open file'
       STOP
     ELSE
       WRITE(100, *) 't ', '                      ', 'x1 ', '           ', 'x2 ', '            ', 'y1 ', '              ', 'y2 '
       DO i = 1, N 
        WRITE(100, *) tlist(i), x1list(i), x2list(i), y1list(i), y2list(i) 
       END DO                       
     CLOSE(100)
     END IF


  





    CALL Verlet_method(y0, t0, x0, dt, N, tlist, x1list, x2list, y1list, y2list)


     OPEN(UNIT=102, IOSTAT=ios, FILE='VerletBonus.txt', STATUS='new', ACTION='write')

     IF (ios /= 0) THEN
       PRINT*, 'ERROR: Unable to open file'
       STOP
     ELSE
      WRITE(102, *) 't ', '            ', 'x1 ', '                ', 'x2 ', '              ', 'y1 ', '                ', 'y2 '
       DO i = 1, N 
        WRITE(102, *) tlist(i), x1list(i), x2list(i), y1list(i), y2list(i)                                
       END DO                       
     CLOSE(102)
     END IF



   

END PROGRAM Diferential_EquationsII

!REAL(4) FUNCTION f(x)
!   REAL(4), DIMENSION(2) :: x
!END FUNCTION f
!

!REAL(4) FUNCTION v(x) 
!    REAL(4), DIMENSION(2) :: x, dy_dt
!    INTEGER :: i
!    REAL :: sum_x
!    DO i = 1, 2
!     !dy_dt(i) =  (-x(i) / ((sqrt(sum((x(i)) ** 2))) ** 3))
!     sum_x = sum_x + x(i)
!     dy_dt(i) = -(x(i)/(sqrt(sum_x ** 2) )** 3)
!   
!    END DO
!END FUNCTION v




REAL(4) FUNCTION v(a, b) 
    REAL(4) :: a, b, dy_dt
     
    dy_dt = -a / (sqrt((a**2 + b**2))) ** 3
   
END FUNCTION v






SUBROUTINE Verlet_method(y0, t0, x0, dt, N, tlist, x1list, x2list, y1list, y2list)
    IMPLICIT NONE
    REAL(4), INTENT(IN) :: t0, dt
    REAL, DIMENSION(2), INTENT(IN) :: x0, y0
    REAL ::  x1, x2, y1, y2
    REAL(4), DIMENSION(300), INTENT(OUT) :: tlist, x1list, y1list, x2list, y2list
    INTEGER(4) :: N, i, j
    REAL :: fy1, fy2
    REAL :: t
    REAL, external :: v 
    

    x1 = x0(1)
    x2 = x0(2)
    y1 = y0(1)
    y2 = y0(2)
    t = t0
    
    
    tlist(1) = t
    x1list(1) = x0(1)
    x2list(1) = x0(2)
    y1list(1) = y0(1)
    y2list(1) = y0(2)

    DO i = 1, N
    
    
    
      fy1 = v(x1, x2)
      y1 = y1 + fy1 * (dt / 2.0)
      x1 = x1 + y1 * dt 
      fy1 = v(x1, x2)
      y1 = y1 + fy1 * (dt / 2.0)
  
  
      
      fy2 = v(x2, x1)
      y2 = y2 + fy2 * (dt / 2.0)
      x2 = x2 + y2 * dt 
      fy2 = v(x2, x1)
      y2 = y2 + fy2 * (dt / 2.0)
      
      
      
      
      t = t + dt
      tlist(i + 1) = t
      x1list(i + 1) = x1
      x2list(i + 1) = x2
      y1list(i + 1) = y1
      y2list(i + 1) = y2
       
    END DO
  
 
END SUBROUTINE Verlet_method




SUBROUTINE Euler_method(y0, t0, x0, dt, N, tlist, x1list, x2list, y1list, y2list)
    IMPLICIT NONE
    REAL(4), INTENT(IN) :: t0, dt
    REAL, DIMENSION(2), INTENT(IN) :: x0, y0
    REAL ::  x1, x2, y1, y2
    REAL(4), DIMENSION(300), INTENT(OUT) :: tlist, x1list, y1list, x2list, y2list
    INTEGER(4) :: N, i, j
    REAL :: fy1, fy2
    REAL :: t
    REAL, external :: f, v 



    x1 = x0(1)
    x2 = x0(2)
    y1 = y0(1)
    y2 = y0(2)
    t = t0
    

    tlist(1) = t   
    x1list(1) = x0(1)
    x2list(1) = x0(2)
    y1list(1) = y0(1)
    y2list(1) = y0(2)

    DO i = 1, N
    
      fy1 = v(x1, x2)
      fy2 = v(x2, x1)
      x1 = x1 + y1 * dt
      y1 = y1 + fy1 * dt
      
      
      x2 = x2 + y2 * dt
      y2 = y2 + fy2 * dt
      
      
      t = t + dt
      tlist(i + 1) = t
      x1list(i + 1) = x1
      x2list(i + 1) = x2
      y1list(i + 1) = y1
      y2list(i + 1) = y2
    END DO
  
 
END SUBROUTINE Euler_method

