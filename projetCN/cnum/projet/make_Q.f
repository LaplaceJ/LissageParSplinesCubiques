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
      SUBROUTINE MAKE_Q ( Q , N_Q, LDQ, VECT_H , N_H)
      IMPLICIT NONE
      
      INTEGER LDQ , N_Q  ,N_H , i 
      DOUBLE PRECISION Q(LDQ,N_Q) , VECT_H(N_H)


     
      DO i = 1 , N_Q 
* r(i-1)
        Q(i,i ) =  3 / VECT_H(i)
* ri
        Q(i+ 2 ,i ) = 3 / VECT_H(i + 1 )
* fi 
        Q(i+ 1 ,i  ) =   - (Q(i,i ) + Q(i + 2 ,i )  ) 
      ENDDO
      
  
      END SUBROUTINE
