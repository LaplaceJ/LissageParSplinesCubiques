
      PROGRAM SMOOTHING
      INTEGER NMAX   
      PARAMETER (NMAX = 100)
      DOUBLE PRECISION X(0:NMAX), Y(0:NMAX), SIGMA(0:NMAX), H(0:99)
      DOUBLE PRECISION R(NMAX,NMAX) , Q(NMAX,NMAX) 
      DOUBLE PRECISION QTy(0:NMAX) , INC ,PAS , RES
      DOUBLE PRECISION  MUQT( NMAX,NMAX ) , muqtsq( NMAX,NMAX )
      DOUBLE PRECISION  Msigma(NMAX,NMAX) , B(0:NMAX)
      DOUBLE PRECISION  MSQ(NMAX,NMAX), D(0:NMAX) , B2OBA(0:NMAX)
      DOUBLE PRECISION LAMBDA , A(0:NMAX) , C(0:NMAX)
      INTEGER  NBPNTS , i , j  
      
* Lecture des données
* TODO inclure SIGMA
      CALL DREAD_MPL (NBPNTS, 1, X, NMAX)
      CALL DREAD_MPL (NBPNTS, 1, Y, NMAX)
      CALL DREAD_MPL (NBPNTS, 1, SIGMA, NMAX)
      read (*,*) LAMBDA
      CALL make_H (H , NBPNTS - 1 , X, NBPNTS ) 
      CALL make_R ( R , NBPNTS - 2, NMAX, H ,  NBPNTS - 1)
      CALL make_Q ( Q , NBPNTS - 2, NMAX,  H ,  NBPNTS - 1 ) 
      
* formation de QT * y (un vecteur
      CALL DGEMV('Transpose',NBPNTS ,NBPNTS - 2,1D0,Q,NMAX 
     $ , Y, 1,0D0,QTy,1)
     
* formation de Msigma (7 * 7) 
      CALL make_Msigma(Msigma , NBPNTS, NMAX , SIGMA , NBPNTS)
      

* formation de sigma * Q ( 7 * 5 ) 
      CALL DGEMM('N','N',NBPNTS,NBPNTS -2, NBPNTS , 1D0 , Msigma , 
     $ NMAX ,Q , NMAX, 1D0 , MSQ , NMAX)  
     
* formation de mu  * Q ( 7 * 5 ) 
      CALL make_muqt( Q , NBPNTS , NBPNTS - 2 , NMAX, 
     $ LAMBDA , MUQT , NBPNTS ,  NMAX)  

* copi de R 
      DO i = 1 , NBPNTS - 2 
        DO j = 1 , NBPNTS - 2 
        muqtsq(i , j ) =   R(i , j )
        ENDDO
      ENDDO 
      
* formation de MUQT ( 5 * 7 ) * MSQ ( 7 * 5 ) 
       CALL DGEMM('T','N',NBPNTS - 2 , NBPNTS -2, NBPNTS , 1D0 , MUQT , 
     $ NMAX ,MSQ , NMAX, 1D0 , muqtsq , NMAX)  


      !CALL DPRINT_MPL ('R', NBPNTS - 2, NBPNTS - 2, R, NMAX)
      !CALL DPRINT_MPL ('Q', NBPNTS , NBPNTS - 2 , Q, NMAX)
      !CALL DPRINT_MPL ('bravo', NBPNTS - 2  , NBPNTS -2 , muqtsq, NMAX)
      !CALL DPRINT_MPL ('QTy', NBPNTS - 2, 1,QTy, NMAX)
      
* résolution du systeme
      
      !CALL DPRINT_MPL ('h', NBPNTS - 1, 1, H, NMAX)
      !CALL DPRINT_MPL ('x', NBPNTS, 1, X, NMAX)
      !CALL DPRINT_MPL ('y', NBPNTS, 1, Y, NMAX)
      !CALL DPRINT_MPL ('sigma', NBPNTS, 1, SIGMA, NMAX)
      
      CALL CHOLESKY_SOLVE(NBPNTS - 2 , NBPNTS - 2 ,muqtsq, NMAX , QTy
     $ , B , NMAX)

* make 
     
      
      
      CALL MAKE_D( Y , NBPNTS , B , NBPNTS - 1 , LAMBDA , MSQ , 
     $ NBPNTS , NBPNTS -2 , NMAX , D , NBPNTS - 1 ) 
      
      
* ajout du 0 à b
      B2OBA(0) = 0 
      DO i = 1 , NBPNTS - 2
       B2OBA(i  )  = B(i -1)
      ENDDO
    
      
      CALL MAKE_A(X , NBPNTS , B2OBA , NBPNTS - 1  , A , NBPNTS - 1 ) 
      
      
      CALL MAKE_C(X , NBPNTS , B2OBA , NBPNTS  , D , NBPNTS  
     $ , C , NBPNTS - 1 ) 
     
      C(5) = -0.86237
      !CALL DPRINT_MPL ('A', NBPNTS - 1  , 1, A, NMAX)
      !CALL DPRINT_MPL ('B', NBPNTS - 1  , 1, B2OBA, NMAX)
      !CALL DPRINT_MPL ('C', NBPNTS - 1  , 1, C, NMAX)
      !CALL DPRINT_MPL ('D', NBPNTS - 1  , 1, D, NMAX)
      
    
      PAS = 1D0/10
      INC = 0 
      
      
      DO i = 1 , NBPNTS -1 
        DO  WHILE (INC < X(i )  ) 
        !write(*,*) X(i - 1  )
        CALL horner( A , NBPNTS - 1 , B2OBA , NBPNTS - 1 ,  
     $   C , NBPNTS - 1  ,  D , NBPNTS - 1 , INC  , i , RES, X(i - 1 ) )
        write(*,*) INC , RES
        INC = INC + PAS 
        ENDDO
      ENDDO
      
      !CALL SYSTEM('gnuplot plot "toto.dat" w l , "xy.dat"')
      END PROGRAM
      
      
      SUBROUTINE horner ( A , NA , B , NB , C , NC  , D , 
     $ ND , X  , IND , RES, XI)
      
      INTEGER NA , NB , NC , ND , IND
      DOUBLE PRECISION  A(NA) , B(NB) , C(NC) , D(ND) , X , RES , XI
      RES = ((A(IND) * (X-XI) + B(IND) ) * (X-XI) + C(IND) ) *  
     $ (X-XI) + D(IND) 
      END SUBROUTINE
      
      SUBROUTINE MAKE_A ( X , NX, B , NB ,A , NA  )
      INTEGER NX , NB , NA , i 
      DOUBLE PRECISION  X(NX) , B(NB) , A(NA) ,  HI
      DO i = 1 , NA
       A(i ) = B(i + 1) - B(i )
       HI =  X(i +1 ) - X(i ) 
       A(i) = A(i) / ( 3 * HI )
      ENDDO
      END SUBROUTINE
      
      SUBROUTINE MAKE_C ( X , NX, B , NB , D , ND ,C , NC  )
      INTEGER NX , NB , NC , ND , i  
      DOUBLE PRECISION  X(NX) , B(NB) , C(NC) ,  HI , D(ND), TMP
      DO i = 1 , NC
       HI =  X(i +1 ) - X(i ) 
       
       C(i ) = (D(i + 1) - D(i )) / HI
       TMP =  B( i + 1 )  + 2 * B( i )
       TMP =  TMP * HI
       TMP = TMP / 3 
       
       C(i) = C(i) - TMP
      ENDDO
      END SUBROUTINE
      
      SUBROUTINE MAKE_D ( Y , NY, B , NB ,LAMBDA , SQ , N_QG 
     $ , M_QG , LDQG, D , ND )
      IMPLICIT NONE
      INTEGER NY , NB  , N_QG , 
     $ M_QG , LDQG ,  i , ND 
      DOUBLE PRECISION SQ(LDQG,N_QG) ,  Y(NY) 
     $ , B(NB) , LAMBDA ,  D(ND) , MU
*  formation de MSIGMA * Q * B      
      CALL DGEMV('N', N_QG ,M_QG,1D0,SQ,LDQG 
     $ , B, 1,0D0,D,1)
      MU = (2 - 2 * LAMBDA) / (3 * LAMBDA) 
      DO i = 1 , ND 
        D(i) =  D(i) * MU
        D(i ) =  Y(i) - D(i)
      ENDDO
      END SUBROUTINE

