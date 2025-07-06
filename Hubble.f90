PROGRAM Hubbles_Constant
  IMPLICIT NONE
  INTEGER :: i, ios
  REAL(8), DIMENSION(15) :: Speed, Distance
  REAL(8) :: sum_x, sum_y, sum_xsqr, sum_xy
  REAL(8) :: Beta_1, Beta_2, SE_Beta_1, SE_Beta_2, Variance, R_squared
  REAL(8) :: x_avg, y_avg, S_res, S_tot
  INTEGER :: N
  REAL(8) :: Matrix(2, 2), Determinant, Determinant_b1, Determinant_b2

  
  REAL(8) :: det

  N = 15

  OPEN(UNIT=100, IOSTAT=ios, FILE='hubble_data.txt', STATUS='old', ACTION='READ')
  IF (ios /= 0) THEN
    PRINT*, 'Error opening the file.'
    STOP
  END IF
  
  READ(100, *)

  
  DO i = 1, N
    READ(100, *) Speed(i), Distance(i)
  END DO
  CLOSE(100)

  
  sum_x = 0.0
  sum_y = 0.0
  sum_xsqr = 0.0
  sum_xy = 0.0

  
  DO i = 1, N
    sum_x = sum_x + Distance(i)
    sum_y = sum_y + Speed(i)
    sum_xsqr = sum_xsqr + Distance(i)**2
    sum_xy = sum_xy + Distance(i) * Speed(i)
  END DO

  
  x_avg = sum_x / N
  y_avg = sum_y / N

  !Here I filled to matrix to apply cramer's rule
  Matrix(1, 1) = N
  Matrix(1, 2) = sum_x
  Matrix(2, 1) = sum_x
  Matrix(2, 2) = sum_xsqr

  
  Determinant = det(Matrix)

  IF (ABS(Determinant) < 1.0D-10) THEN
    PRINT*, 'The system does not have a unique solution (determinant is zero).'
    STOP
  END IF

  
  Matrix(1, 1) = sum_y
  Matrix(2, 1) = sum_xy
  Determinant_b1 = det(Matrix) / Determinant

  Matrix(1, 1) = N
  Matrix(2, 1) = sum_x
  Matrix(1, 2) = sum_y
  Matrix(2, 2) = sum_xy
  Determinant_b2 = det(Matrix) / Determinant

  
  Beta_1 = Determinant_b1
  Beta_2 = Determinant_b2

  
  S_res = 0.0
  S_tot = 0.0

  DO i = 1, N
    S_res = S_res + (Speed(i) - (Beta_1 + Beta_2 * Distance(i)))**2
    S_tot = S_tot + (Speed(i) - y_avg)**2
  END DO

  
  Variance = S_res / (N - 2)

  
  SE_Beta_2 = SQRT(Variance / SUM((Distance - x_avg)**2))
  SE_Beta_1 = SQRT(Variance * (1.0 / N + x_avg**2 / SUM((Distance - x_avg)**2)))

  
  R_squared = 1.0 - (S_res / S_tot)

  
  OPEN(UNIT=101, FILE='fit.txt', STATUS='REPLACE', ACTION='WRITE')
  WRITE(101,*) 'Beta_1: ', Beta_1
  WRITE(101,*) 'Beta_2: ', Beta_2
  WRITE(101,*) 'SE_Beta_1: ', SE_Beta_1
  WRITE(101,*) 'SE_Beta_2: ', SE_Beta_2
  WRITE(101,*) 'R_squared: ', R_squared
  CLOSE(101)

END PROGRAM Hubbles_Constant



REAL(8) FUNCTION det(matrix)
  IMPLICIT NONE
  REAL(8), DIMENSION(2, 2), INTENT(IN) :: matrix
  det = matrix(1, 1) * matrix(2, 2) - matrix(1, 2) * matrix(2, 1)
END FUNCTION det
