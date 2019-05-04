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

      SUBROUTINE make_muqt( Q , NQT, MQT, LDQT,  
     $ LAMBDA, RES , NR, LDR )
      IMPLICIT NONE
      
      INTEGER LDQT , NQT , MQT , i , j ,  NR , LDR
      DOUBLE PRECISION Q(LDQT,NQT) , MU , 
     $ RES(LDR , NR) , LAMBDA 
      
      MU = (2 - 2 * LAMBDA) / (3 * LAMBDA) 
     
      DO i = 1 , NQT 
        DO j = 1 , MQT 
        RES(i , j ) =   Q(i,j) * MU 
        ENDDO
      ENDDO
      
  
      END SUBROUTINE
