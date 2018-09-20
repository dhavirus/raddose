C-----------------------------------------------------------------------        
      SUBROUTINE UPCASE(LINE)
C     CONVERT LOWER-CASE CHARACTERS TO UPPER CASE:                              
      CHARACTER LINE*(*)
      DO 100 N=1,LEN(LINE)
        IC=ICHAR(LINE(N:N))
        IF(IC.GE.97.AND.IC.LE.122) LINE(N:N)=CHAR(IC-32)
 100        CONTINUE
      END

