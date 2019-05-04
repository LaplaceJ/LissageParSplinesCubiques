* Cette fonction calcul la matrice R
*
*
*
*
*
* Entier : 
* N_R nombre de ligne de la matrice 
* N_H nombre ligne du vecteur H 
* LDR La dimension principale de R
* i indice de boucle
*
* DOUBLE PRECISION
* R matrice à calculer 
* VECT_H vecteur qui permet de calculer R

      SUBROUTINE MAKE_R ( R , N_R, LDR, VECT_H , N_H)
      IMPLICIT NONE
      
      INTEGER LDR , N_R  ,N_H , i 
      DOUBLE PRECISION R(LDR,N_R) , VECT_H(N_H)
* p      
      DO i = 1 , N_R 
        R(i,i) =  2 * (VECT_H(i)  + VECT_H(i + 1 ) )
      ENDDO
* h      
      DO i = 1 , N_R  - 1 
      R(i,i + 1 ) = VECT_H(i + 1)
      R(i + 1 ,i) = VECT_H(i + 1)
      ENDDO
      
  
      END SUBROUTINE
