* fonction factorisant une matrice sous la 
* forme de cholesky.  A = btb
*
* Cette fonction calcul la matrice b diagonal supérieur 
* ou inférieur.
*
*
*
*
*
* Entier : 
* UPLO matrice supérieur(1) ou inférieur(0)
* N nombre de colonnes de la matrice 
* LDA La dimension principale de A
* i, j, k indices de boucle
*
* DOUBLE PRECISION
* A matrice à factoriser
* tmp variable tampon

* write(*,*)  2 * (VECT_H(i)  + VECT_H(i + 1 ) )

      SUBROUTINE MAKE_Msigma( Msigma , N_MS, LDMS, VECT_S , N_S)
      IMPLICIT NONE
      
      INTEGER LDMS , N_MS  ,N_S , i 
      DOUBLE PRECISION Msigma(LDMS,N_MS) , VECT_S(N_S)


     
      DO i = 1 , N_MS 
        Msigma(i,i ) =   VECT_S(i) 
      ENDDO
      
  
      END SUBROUTINE
