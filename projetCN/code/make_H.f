*
* Cette fonction calcul le vecteur  H
*
*
*
*
*
* Entier : 
* N_H taille du vecteur H
* N_X taille du vecteur X ( n_h + 1 ) 
* i indice de boucle
*
* DOUBLE PRECISION
* VECT_H vecteur de ???? 
* VECT_X vecteur des abscisse 


      SUBROUTINE MAKE_H (VECT_H , N_H , VECT_X , N_X )
      IMPLICIT NONE
      
      INTEGER N_H , N_X , i 
      DOUBLE PRECISION VECT_H(N_H), VECT_X(N_X) 
      
      
      DO i = 1 , N_H 
       !calcul de VECT_H 
       VECT_H(i) = VECT_X(i + 1 ) - VECT_X(i) 
      ENDDO
      
      END SUBROUTINE
