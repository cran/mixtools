      SUBROUTINE MUDEPTH(N,T,D,MPT,X,COUNT,SDEP)
C
      INTEGER N,D,T,COUNT(T)
C
      REAL MPT(T,D),X(N,D),SDEP(T)
      INTEGER I,J,K,M,L
      REAL D1,D2,D3,D4,D5
C
      DO 200 L=1,T
C
      COUNT(L)=0
      SDEP(L)=0.0
C  
C     COUNT HOW MANY CIRCLES WILL COVER THE POINT
C
      M=N-1
      DO 10 I=1,M
         DO 50 J=(I+1),N
            D1=0.0
            D2=0.0
            D3=0.0
            DO 100 K=1,D
               D1=D1+(X(I,K)-MPT(L,K))**2
               D2=D2+(X(J,K)-MPT(L,K))**2
               D3=D3+(X(I,K)-X(J,K))**2
100   CONTINUE
               D4=D1+D2-D3
C
               IF (D4.LE.0.0) THEN
                  COUNT(L)=COUNT(L)+1
               ENDIF
C
50    CONTINUE
10    CONTINUE
      D5=(N*(N-1)/8)**(.5)
      SDEP(L)=(COUNT(L)-N*(N-1)/4)/D5
C
200   CONTINUE
C
      RETURN
      END
C
C
C
