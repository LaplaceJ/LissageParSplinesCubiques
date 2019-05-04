       SUBROUTINE CHOLESKY_SOLVE(M, N, A, LDA, B ,  RES , NMAX)
* Les paramètres formels (LDA = Leading Dimension of A)
      INTEGER M, N, LDA, i , NMAX
      DOUBLE PRECISION A(LDA,*), B(*)
* Variables locales
      DOUBLE PRECISION ATA(N,N), ATB(N) , RES(NMAX)
      INTEGER INFO
*
      !WRITE (*,*) 'Résolution par la méthode du Cdt Cholesky'
* ATA = Transpose (A) . A   (DGEMM est une BLAS pour le pdt de matrices)
      CALL DGEMM ('Transpose', 'No Transpose', 
     $             N, N, M, 1D0, A, LDA, A, LDA, 0D0, ATA, N)
      !CALL DPRINT_MPL ('A**T . A', N, N, ATA, N)
* ATB = Transpose (A) . b   (DGEMV est une BLAS pour le pdt matrice . vecteur)
      CALL DGEMV ('Transpose', M, N, 1D0, A, LDA, B, 1, 0D0, ATB, 1)
      !CALL DPRINT_MPL ('A**T . b', N, 1, ATB, N)
* Factorisation de Cholesky : ATA = L . Transpose (L)
* DPOTRF est une fonction LAPACK pour la méthode de Cholesky
      CALL DPOTRF ('Lower', N, ATA, N, INFO)
      IF (INFO .NE. 0) STOP 'Erreur DPOTRF'
      !CALL DPRINT_MPL ('Après DPOTRF, A', N, N, ATA, N)
* Résolution L . y = ATB avec résultat dans ATB.
* DTRSV est une BLAS pour la substitution avant/arrière
      CALL DTRSV ('Lower', 'No Transpose', 'Not Unit Triangular',
     $              N, ATA, N, ATB, 1)
* Résolution Transpose (L) . x = y (résultat dans ATB)
      CALL DTRSV ('Lower', 'Transpose', 'Not Unit Triangular',
     $              N, ATA, N, ATB, 1)
      
      DO i = 1 , N + 1 
        RES(i) = ATB(i) 
      ENDDO
* CALL DPRINT_MPL ('solution', N, 1, ATB, N)
      END SUBROUTINE
