PROGRAM Black_Body
  USE Module_Ass10
  IMPLICIT NONE
  REAL(4), DIMENSION(1612) :: Intensity, Wavelenght, Residual   
  REAL(4), DIMENSION(1612, 2) :: Jacobian      
  REAL(4), DIMENSION(2, 2) :: Matrix, Inverse_Matrix, Var_Beta  
  REAL(4), DIMENSION(2) :: Delta_Beta, Betas, Auxiliar_Matrix, SE_Beta  
  INTEGER(4) :: ios, i, N, k, maxIter          
  REAL(4) :: Beta_1, Beta_2, Determinant, Alpha, Tolerance, S_Residual, Temperature 



  MaxIter = 1000  
  Alpha = 0.1         
  Tolerance = 1.0E-6      
  N = 1612        


  
  OPEN(UNIT=100, IOSTAT=ios, FILE='sun_data.txt', STATUS='old', ACTION='read')
  IF (ios /= 0) THEN
    PRINT*, 'ERROR: Unable to open file'
    STOP
  ELSE
    READ(100,*)  
    READ(100,*)  
    DO i = 1, N
      READ(100, *) Wavelenght(i), Intensity(i)  
    END DO
    CLOSE(100)  
  END IF



  
  ! Initial guess
  Betas = (/1.0, 1.0/)

  ! Gauss-Newton code
  DO k = 1, MaxIter
  
    DO i = 1, N
      Residual(i) = Intensity(i) - B(Wavelenght(i), Betas(1), Betas(2))  
      Jacobian(i, 1) = Partial_Derivate_Beta1(Wavelenght(i), Betas(2))                    
      Jacobian(i, 2) = Partial_Derivate_Beta2(Wavelenght(i), Betas(1), Betas(2))    
    END DO

    !  (Matrix = Jacobian^T * Jacobian)
    Matrix = MATMUL(TRANSPOSE(Jacobian), Jacobian)
        
    Determinant = (Matrix(1,1) * Matrix(2,2)) - (Matrix(1,2) * Matrix(2,1))
   !(Inverse_Matrix = M^-1)
    Inverse_Matrix(1,1) = Matrix(2,2) / Determinant
    Inverse_Matrix(1,2) = -Matrix(2,1) / Determinant
    Inverse_Matrix(2,1) = -Matrix(1,2) / Determinant
    Inverse_Matrix(2,2) = Matrix(1,1) / Determinant

    !J^T*r = Auxiliar_Matrix


    Auxiliar_Matrix = MATMUL(TRANSPOSE(Jacobian), residual)
    Delta_Beta = MATMUL(Inverse_Matrix, Auxiliar_Matrix)  


    IF (SUM(Delta_Beta**2) <= Tolerance) EXIT

 
    Betas = Betas + Alpha * Delta_Beta


  END DO 






  S_Residual = SUM(Residual**2) 
  Var_Beta = (S_Residual / (N - 2)) * Inverse_Matrix
  SE_Beta(1) = SQRT(Var_Beta(1,1))  
  SE_Beta(2) = SQRT(Var_Beta(2,2))  



  Temperature = 14387.77 / Betas(2)  








OPEN(UNIT=101, IOSTAT=ios, FILE='fit_Ass10.txt', STATUS='new', ACTION='write')

IF (ios /= 0) THEN
  PRINT*, 'ERROR: Unable to open file'
  STOP
ELSE
  WRITE(101, *) 'Beta_1: ', Betas(1)                      
  WRITE(101, *) 'Beta_2: ', Betas(2)                 
  WRITE(101, *) 'Standard Error of Beta_1: ', SE_Beta(1)        
  WRITE(101, *) 'Standard Error of Beta_2: ', SE_Beta(2)        
  WRITE(101, *) 'Estimated Surface Temperature of the Sun:', Temperature           

  CLOSE(101)
  END IF

END PROGRAM Black_Body
