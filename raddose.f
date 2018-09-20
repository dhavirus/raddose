C
C       RADDOSE
C	
C 	A program to calculate the absorbed X-ray dose of macromolular
C       crystals in X-ray beams
C
C       James Murray, Elspeth Garman and Raimond Ravelli
C       using the mucal.f routine from Pathikrit Bandyopadhyay at 
C 	the University of Notre Dame


	PROGRAM RADDOSE
C	==============
	IMPLICIT NONE

C       error variables
	
	INTEGER IERR, LFLAG, LENSTR
	CHARACTER PLOTF*10

C       things for parser
	INTEGER MAXTOK
	PARAMETER (MAXTOK=200)
	CHARACTER KEY*4,CVALUE(MAXTOK*4)*4,LINE*256
	INTEGER IBEG(MAXTOK), IEND(MAXTOK), ITYP(MAXTOK), IDEC(MAXTOK), NTOK
	REAL FVALUE(MAXTOK)
	LOGICAL LEND, ANISFLG
	INTEGER WEBO

C       Miscellaneous variables
	INTEGER N,K,I,IUNOUT
	CHARACTER NAME*2
	REAL V, X

C       MATHEMATICAL AND PHYSICAL CONSTANTS
	REAL PI
	PARAMETER (PI=3.1415926536)
	REAL NAVO,AMU, PDENSI
	REAL RDENSI, DDENSI
	PARAMETER (NAVO=6.022e+23)
C       Atomic mass unit 1.66*10^-24 g
	PARAMETER (AMU=1.66E-24)
C       Protein density 1.35 g / ml
	PARAMETER (PDENSI=1.35)
C       RNA density  g / ml 
	PARAMETER (RDENSI=2.00)
C       DNA density g / ml 
	PARAMETER (DDENSI=2.00)
C       Mass in Da of "average" amino acid
	REAL AAMASS
C       Mass in Da of "average" nucleotide RNA,DNA
	REAL NPMASS, DNMASS
	PARAMETER (AAMASS=110.0)
C       Mass in Da of "average" deoxynucleotide DNA
	PARAMETER (DNMASS=312.0)
C       Mass in Da of "average" nucleotide RNA
	PARAMETER (NPMASS=321.0)
C Henderson limit default 2E7 Gray
	REAL HENDLI 

C	Factor between keV and A
	REAL KEVA
	PARAMETER (KEVA=12.398)
C	Copper K alpha wavelength
 	REAL CUKA 
	PARAMETER (CUKA=1.541838)

C Convective heat transfer coefficient W m**-2 K**-1 ~320
	REAL CONVHE
C Heat Capacity J/kg/K ~500
	REAL CP
C Heat transfer things
	REAL DELTAT,TSYS,PABS
	
C	No of elements
	INTEGER ELEMEN
	PARAMETER (ELEMEN=95)

C	Beam and experimental parameters	
C	Beam energy in keV or wavelength in A
	REAL EN
	REAL WAVE
	REAL PHOSEC,CFLUX
	REAL XBEA, YBEA, SIGX,SIGY
	REAL XFWHM, YFWHM,CORREF
	REAL XCRY, YCRY, ZCRY
	INTEGER IMAGES
	REAL EXPO
	REAL LRANGE, URANGE, GSTEP

C       Crystal parameters
	REAL CELL(6)
	REAL PROT(ELEMEN), SOLV(ELEMEN), SOLVMM(ELEMEN), TATM(ELEMEN)
	REAL SOLVEN, VOLUME, ULT, PATM, NWAT, NDRY
	INTEGER NRES, NMON, Z, NELEM
	INTEGER NRNA, NDNA
	REAL MASS,DENS,MU, FRAC, FRACC
	REAL CONC

C	Crystal absorption parameters, cross-sections in barns
	REAL XSEC_TO(ELEMEN), XSEC_PE(ELEMEN),XSEC_IN(ELEMEN), XSEC_CO(ELEMEN)	
	REAL DENSI, DENSIG, GRYABS, HENDT,PCABS

C	Total cross-sections for the whole crystal
	REAL XSCR_TO, XSCR_PE, XSCR_CO, XSCR_IN
C	Cross-sections for protein
	REAL XSPR_TO, XSPR_PE, XSPR_CO, XSPR_IN
C	Cross-sections for solvent
	
	
	REAL CVOL,BEAMF,PCELLH,PCELLS,IMDOSE,DIDOSE
	REAL CELVOL

C	Parameters for mucal
	REAL XSEC(10), EDGES(9), FLY(4)
	INTEGER ER
	CHARACTER UNIT
	LOGICAL ERF
	CHARACTER*2 MANE

	EXTERNAL  CCPERR,CCPFYP,CCPRCS,MEMOPARSE, NAME2Z, Z2NAME,
     +            MUCAL,CCPDPN,QCLOSE,LENSTR

	COMMON /SPLINE/ SPLIF,ESPLI,ZSPLI,SPLI
C	COMMON /MUCAL/ 
	REAL SPLI(2,1000)
	CHARACTER SPLIF*132
	CHARACTER ESPLI*2
	INTEGER ZSPLI



C===================================================================
C       INITIALIZATIONS
C==================================================================
	CALL CCPFYP
	CALL CCP4H_INIT()
	CALL CCPRCS(6,'RADDOSE','$DATE: 2003/01/29 12:30:24  $')
	
	IERR = 0
	LFLAG = 0

	IUNOUT = 0 

	NAME = '  '

C       Program parameters
C       WEBO - WEBOUTPUT 1 for cgi script version
	WEBO = 0

C Parameters of the crystal
	PATM=0.0
	NRES = 0
	NMON = 0 
	CVOL = -1.0
	SOLVEN = -1.0
	ULT = 0.0
	MASS = 0 
	DENSI = 0
	MU = 0.0
	NWAT = 0.0
	NDRY = 0.0
	DO 10 I=1,ELEMEN
	   PROT(I)=0.0
	   SOLV(I)=0.0
	   SOLVMM(I)=0.0
	   TATM(I)=0.0
 10	CONTINUE

C SPLINOR 
	ESPLI='XX' 
	ZSPLI=999
C Experimental parameters
C	BEAM energy
	EN = -1.0
C 	Copper K alpha default
	WAVE = CUKA 
	XBEA = 0.1
	YBEA = 0.1
	XCRY = 0.1
	YCRY = 0.1
	ZCRY = 0.1
	XFWHM = -1.0
	YFWHM = -1.0
	CORREF = 1.0
	PHOSEC = 1.0E9
	IMAGES = 90
	EXPO = 10.0
	HENDT = 0.0
	HENDLI = 2E7

	LRANGE = -1.0
	URANGE = 20.0
	GSTEP = 0.01
	
C       CP = 500 J / kg
	CP = 500
C       CONVHE = 100 W m-2 K-1
	CONVHE = 320
	PABS = 0.0

C	mucal stuff
	UNIT = 'B'
	ERF = .TRUE.

C====================================================================
C       Read Command input
C====================================================================
 100	line=' ' 
	NTOK=MAXTOK
C       call parser in parser.f to parse command line(s) 
	CALL PARSER(
     +  KEY,LINE,IBEG,IEND,ITYP,FVALUE,
     +	CVALUE,IDEC,NTOK,LEND,.TRUE.)
C       Check for end-of-file lend ==.true.
	if(LEND) go to 200
C       Check for blank line, if so try again
	if(NTOK.eq.0) GO TO 100

C Start parsing input keywords

C..................................................................
C end of input tokens
C END
C..................................................................
	IF (KEY .eq. 'END') THEN 
	   GO TO 200
C..................................................................
C optional REMARK echo to standard out
C REMARK
	ELSEIF (KEY .EQ. 'REMA') THEN
	   WRITE(*, '(a)')  line(1:lenstr(line)) 

C...................................................................
C Get unit cell parameters for volume calculation if not given already
C CELL 
	ELSEIF (KEY .EQ. 'CELL') THEN 
	   IF (NTOK .LT. 4) THEN
	      WRITE(6, '(a)') '** Too few numbers on CELL command **'
	      IERR = IERR + 1
	      GO TO 100
	   ENDIF
	   CELL(4) = 90
	   CELL(5) = 90
	   CELL(6) = 90
	   DO 110, I=1,6
	      CALL GTTREA(I+1,CELL(I),LFLAG,NTOK,ITYP,FVALUE)
	      IF (LFLAG .NE. 0 .AND. I .LE. 3) GO TO 199
 110	   CONTINUE
	   CVOL = -1.0
	   
C..................................................................
C VOLUME
C..................................................................
	ELSEIF (key .EQ. 'VOLU') THEN
	   n = 2
	   CALL GTTREA(N,CVOL,LFLAG,NTOK,ITYP,FVALUE)
	   IF (lflag .GT. 0) GO TO 199
	   
C.............................................................
C SOLVENT
C................................................................
	ELSEIF (key .EQ. 'SOLV') THEN 
	   n = 2   
	   CALL GTTREA(N,SOLVEN,LFLAG,NTOK,ITYP,FVALUE)
	   IF (lflag .GT. 0) GO TO 199

C.................................................................
C No. of residues in a monomer
C NRES
	ELSEIF (key .EQ. 'NRES') THEN
	   n = 2   
           CALL GTTINT(N,NRES,LFLAG,NTOK,ITYP,FVALUE)
	   IF (lflag .GT. 0) GO TO 199

C.................................................................
C No. of nucleotides of RNA present
C NRNA
	ELSEIF (key .EQ. 'NRNA') THEN
	   n = 2   
           CALL GTTINT(N,NRNA,LFLAG,NTOK,ITYP,FVALUE)
	   IF (lflag .GT. 0) GO TO 199	   

C.................................................................
C No. of nucleotides of DNA present
C NRNA
	ELSEIF (key .EQ. 'NDNA') THEN
	   n = 2   
           CALL GTTINT(N,NDNA,LFLAG,NTOK,ITYP,FVALUE)
	   IF (lflag .GT. 0) GO TO 199	   

C.................................................................
C No. of monomers in unit cell
C NMON
	ELSEIF (key .EQ. 'NMON') THEN
	   n = 2 
	   CALL GTTINT(N,NMON,LFLAG,NTOK,ITYP,FVALUE)
	   IF (LFLAG .GT. 0) GO TO 199


C................................................................
C Extra atoms elements in protein - no atoms per monomer
C eg PATM C 5 P 16 N 17 Se 1 
C PATM
C............................................................... 	      
	ELSEIF (key .EQ. 'PATM') THEN	   
	   NELEM = ( NTOK - 1 ) / 2 
	   DO 120, I=1,NELEM
	      N = 2*I + 1	
	      NAME=CVALUE(2*I)
	      Z = 0

	      CALL CCPUPC(NAME)
	      CALL NAME2Z(NAME,Z)
              CALL GTTREA(N,PATM,LFLAG,NTOK,ITYP,FVALUE)
	      IF (LFLAG .GT. 0) GO TO 199
	      PROT(Z)= PROT(Z) + PATM	
120	   CONTINUE

C...............................................................
C SATM
C...............................................................
	   ELSEIF (KEY .EQ. 'SATM') THEN
           NELEM = ( NTOK - 1 ) / 2
           DO 130, I=1,NELEM
              N = 2*I + 1
              NAME=CVALUE(2*I)
	      Z = 0
              CALL CCPUPC(NAME)
              CALL NAME2Z(NAME,Z)
              CALL GTTREA(N,PATM,LFLAG,NTOK,ITYP,FVALUE)
              IF (LFLAG .GT. 0) GO TO 199
              SOLVMM(Z)= SOLVMM(Z) + PATM
130	   continue


C............................................................
C ENERGY keV 
C.........................................................
	ELSEIF (key .EQ. 'ENER') THEN
	 n = 2
	 CALL GTTREA(n,EN,lflag,ntok,ityp,fvalue)
         IF (lflag .GT. 0) GO TO 199

C............................................................
C HENDLI Energy dose limit defaults to 2E7 Gray
C.........................................................
	ELSEIF (key .EQ. 'HEND') THEN
	 n = 2
	 CALL GTTREA(n,HENDLI,lflag,ntok,ityp,fvalue)
         IF (lflag .GT. 0) GO TO 199


C...........................................................
C RANGE to plot for graph - low limit, upper limit, step size
C energies in keV
C..........................................................
	ELSEIF(key .EQ. 'RANG') THEN
	CALL GTTREA(2,LRANGE,lflag,ntok,ityp,fvalue)
	IF (lflag .GT. 0) GO TO 199
	CALL GTTREA(3,URANGE,lflag,ntok,ityp,fvalue)
	IF (lflag .GT. 0) GO TO 199
	IF (NTOK .EQ. 4) THEN 
	   CALL GTTREA(4,GSTEP,lflag,ntok,ityp,fvalue)
	   IF (lflag .GT. 0) GO TO 199
	ENDIF

C...........................................................
C PLOTFILE
C............................................................
	ELSEIF (key .EQ. 'PLOT') THEN
	   IF (NTOK .GT. 1) THEN
	      k = ibeg(2)+10
	      plotf = line(ibeg(2):iend(2))
	   ENDIF

C............................................................
C WAVELENGTH A
C............................................................
	ELSEIF (key .EQ. 'WAVE') THEN
	n = 2
	CALL GTTREA(n,WAVE,lflag,ntok,ityp,fvalue)
	IF (lflag .GT. 0) GO TO 199

C.............................................................
C PHOSEC photons s-1
C...........................................................
	ELSEIF (key .EQ. 'PHOS') THEN
	n = 2
	CALL GTTREA(n,PHOSEC,lflag,ntok,ityp,fvalue)
	IF (lflag .GT. 0) GO TO 199	

C.............................................................
C BEAM xbeam ybeam size in mm or radius in mm
C.............................................................
	ELSEIF (key .EQ. 'BEAM') THEN
		IF (ntok .EQ. 3) THEN
	           CALL GTTREA(2,XBEA,lflag,ntok,ityp,fvalue)
		   IF (lflag .GT. 0) GO TO 199
 		   CALL GTTREA(3,YBEA,lflag,ntok,ityp,fvalue)	
		   IF (lflag .GT. 0) GO TO 199
		ENDIF
C	IMPLEMENT CIRCULAR BEAM LATER

C.............................................................
C GAUSS xbeam-FWHM ybeam-FWHM in mm
C.............................................................
	ELSEIF (key .EQ. 'GAUS') THEN
		IF (ntok .EQ. 3) THEN
	           CALL GTTREA(2,XFWHM,lflag,ntok,ityp,fvalue)
		   IF (lflag .GT. 0) GO TO 199
 		   CALL GTTREA(3,YFWHM,lflag,ntok,ityp,fvalue)	
		   IF (lflag .GT. 0) GO TO 199
		ENDIF



C.............................................................
C CRYSTAL xtal dimensions x y z in mm as seen by the beam
C.............................................................
	ELSEIF (key .EQ. 'CRYS') THEN
		IF (ntok .EQ. 4) THEN
                    CALL GTTREA(2,XCRY,lflag,ntok,ityp,fvalue)
	 	    IF (lflag .GT. 0) GO TO 199	    
	            CALL GTTREA(3,YCRY,lflag,ntok,ityp,fvalue)
                    IF (lflag .GT. 0) GO TO 199
                    CALL GTTREA(4,ZCRY,lflag,ntok,ityp,fvalue)
                    IF (lflag .GT. 0) GO TO 199
		ELSE
		   GO TO 199
		ENDIF

C..............................................................
C IMAGES no of images taken
C.............................................................
	ELSEIF (key .EQ. 'IMAG') THEN
	  CALL GTTINT(2,IMAGES,lflag,ntok,ityp,fvalue)
	  IF (lflag .GT. 0)  GO TO 199

C..............................................................
C EXPOSURE exposure time per image
C.............................................................
	ELSEIF (key .EQ. 'EXPO') THEN
	   CALL GTTREA(2,EXPO,lflag,ntok,ityp,fvalue)
	    IF (lflag .GT. 0)  GO TO 199

C...............................................................
C CONV convective heat transfer coefficient W m**-2 K**-1
C..............................................................
	ELSEIF (key .EQ. 'CONV') THEN
	   CALL GTTREA(2,CONVHE,lflag,ntok,ityp,fvalue)
	     IF (lflag .GT. 0) GO TO 199

C..............................................................
C HEAT  Crystal Heat capacity J / kg / K
C..............................................................
	ELSEIF (key .EQ. 'HEAT') THEN
	   CALL GTTREA(2,CP,lflag,ntok,ityp,fvalue)
	     IF (lflag .GT. 0) GO TO 199

C................................................................
C SPLINOR  chooch splinor file element and filename
C...............................................................
	ELSEIF (key .EQ. 'SPLI') THEN	   
	     ESPLI=CVALUE(2)
	     SPLIF=CVALUE(3)

	      IF (NTOK .GT. 1) THEN
		 SPLIF = line(ibeg(3):iend(3))
	      ENDIF

	      CALL CCPUPC(ESPLI)
	      WRITE (*,*) ESPLI

	      IF (LFLAG .GT. 0) GO TO 199
140	   CONTINUE

C...............................................................
C WEBO Weboutput toggle
C..............................................................
	ELSEIF (key .EQ. 'WEBO') THEN
	   CALL GTTREA(2,WEBO,lflag,ntok,ityp,fvalue)
	     IF (lflag .GT. 0) GO TO 199

C...............................................................
	ELSE
	   WRITE (6, '(/'' Unknown keyword: '',A4/)') key
	ENDIF
	GO TO 100
C...............................................................
 199	IERR = IERR + 1
	GO TO 100
C---------------------------------------------------------------
C END OF INPUT PARSING
C--------------------------------------------------------------	


C===============================================================
C Check Commands
C===============================================================
 200	IF (IERR .GT. 0) THEN
	   CALL CCPERR(1,'** Input Error **')
	ENDIF

C........................................................
C Sort out energy and wavelength	
C.......................................................
	IF (EN .LT. 0.0) THEN
		EN = KEVA / WAVE
	ENDIF	

C.........................................................
C Correct the beam flux for non-uniform beam
C........................................................
	CORREF = BEAMF(XBEA,YBEA,XFWHM,YFWHM)

C........................................................
C Read SPLINOR file if necessary
C........................................................
	IF (ESPLI .NE. 'XX') THEN
	CALL RSPLIN 
	ENDIF




C==============================================================
C Contents of Crystal
C==============================================================

C...............................................................
C Calculate unit cell volume
	IF (CVOL .LE. 0.0) THEN
	   CVOL = CELVOL(CELL)
	ENDIF

C--------------------------------------------------
C Estimate Solvent content if not given explicitly
C--------------------------------------------------
C - density protein = 1.35 g /ml
C - density RNA     = 2.0 g /ml
C - density DNA     = 2.0 g /ml
C 1 ml = 10^24 A**3
C mass = 1.66E-24 * 110* NRES * NMOL 
	IF (SOLVEN .LT. 0) THEN
C PROTEIN
	   X = AMU * AAMASS * REAL(NRES) * REAL(NMON) /(CVOL * PDENSI * 1E-24 )
	   SOLVEN = 1 - X
C RNA
	   X =  AMU * NPMASS  * REAL(NRNA) *  REAL(NMON) /(CVOL * RDENSI * 1E-24 )
	   SOLVEN = SOLVEN - X
C DNA
	   X =  AMU * DNMASS  * REAL(NDNA) *  REAL(NMON) /(CVOL * DDENSI * 1E-24 )
	   SOLVEN = SOLVEN - X
	   IF (SOLVEN .LE. 0.0) CALL CCPERR(1, 
     + '** Solvent Content estimated to be negative!**')
	ENDIF


C.....................................................
C Turn solvent atom concentrations (mM) into no atoms per cell
C.....................................................
C Assume solvent is filled with 55M H2O and that 
C each non H20 atom displaces a water
C no of atoms in the cell	
C 1 litre = 1E27 A**3
	do 210 I=1,elemen
	SOLV(I) = SOLVMM(I) * 1E-3 * NAVO * CVOL * 1E-27 * 
     +           SOLVEN
	NDRY = NDRY + SOLV(I)
210 	CONTINUE	
	NWAT = 55000 * 1E-3  * NAVO * CVOL * 1E-27 *
     +         SOLVEN - NDRY 
C	Add oxygen and hydrogen from water to solvent
C       Water is H2O
	SOLV(8) = SOLV(8) + NWAT
	SOLV(1) = SOLV(1) + 2* NWAT

C....................................................
C take into account NMON for PATM atoms
C....................................................
	DO 211 I=1,ELEMEN
	PROT(I) = PROT(I) * REAL(NMON)
211	CONTINUE


C...................................................
C Add atoms from NRES and NMON to PROT(N)
C....................................................
C Use approx amino acid = 5C + 1.35 N + 1.5 O + 8H
C Carbon
	PROT(6) = PROT(6) +  5    * REAL(NRES) * REAL(NMON)
C Nitrogen
	PROT(7) = PROT(7) +  1.35 * REAL(NRES) * REAL(NMON)
C Oxygen	
	PROT(8) = PROT(8) +  1.5  * REAL(NRES) * REAL(NMON)
C Hydrogen
	PROT(1) = PROT(1) +  8    * REAL(NRES) * REAL(NMON) 

C.......................................................
C Add atoms from NDNA and NRNA and NMON to PROT(N)
C.......................................................
C RNA
C ANP = 13 H + 10 C + 5 N + 6 O + 1 P  
C CNP = 11 H + 9  C + 3 N + 7 O + 1 P  
C GNP = 11 H + 10 C + 5 N + 7 O + 1 P  
C UNP = 10 H + 9 C  + 2 N + 8 O + 1 P
C mean NP = 11.25 H + 9.5 C + 3.75 N + 7 O +1 P 
C Carbon
	PROT(6) = PROT(6) +  9.5 * REAL(NRNA) * REAL(NMON)
C Nitrogen
	PROT(7) = PROT(7) +  3.75 * REAL(NRNA) * REAL(NMON)
C Oxygen	
	PROT(8) = PROT(8) +  7  * REAL(NRNA) * REAL(NMON)
C Hydrogen
	PROT(1) = PROT(1) +  11.25 * REAL(NRNA) * REAL(NMON) 
C Phosphorus 
	PROT(15) = PROT(15) +  1 * REAL(NRNA) * REAL(NMON) 
C DNA
C dNP
C dANP = 13 H + 10 C + 5 N + 5 O + 1 P  
C dCNP = 11 H + 9  C + 3 N + 6 O + 1 P  
C dGNP = 11 H + 10 C + 5 N + 6 O + 1 P  
C dTNP = 12 H + 10 C + 5 N + 7 O + 1 P  
C mean dNP = 11.75 H + 9.75 C + 4  N + 6 O +1 P 
C Carbon
	PROT(6) = PROT(6) +  9.75    * REAL(NDNA) * REAL(NMON)
C Nitrogen
	PROT(7) = PROT(7) +  4 * REAL(NDNA) * REAL(NMON)
C Oxygen	
	PROT(8) = PROT(8) +  6  * REAL(NDNA) * REAL(NMON)
C Hydrogen
	PROT(1) = PROT(1) +  11.75    * REAL(NDNA) * REAL(NMON) 
C Phosphorus 
	PROT(15)  = PROT(15) +  1  * REAL(NDNA) * REAL(NMON) 


C--------------------------------------------------------
C Sum protein and solvent contributions 
C---------------------------------------------------------
	DO 215 I=1,ELEMEN
	TATM(I) = PROT(I) + SOLV(I)
215	CONTINUE


C---------------------------------------------------
C Call the mucal routine to get X-ray cross-sections
C Put this in a subroutine, so i can iterate over an array
C of energies
C----------------------------------------------------

C- Call mucal
	MU = 0
	XSCR_PE=0
	DO 220, I=1,ELEMEN
	IF (TATM(I) .GT. 0) THEN
	CALL Z2NAME(NAME,I)
	CALL MUCAL(EN,NAME,I,UNIT,XSEC,EDGES,FLY, ERF,ER)
C            ============================================	
	MASS = MASS + TATM(i)*XSEC(7)*AMU
	
C Empirical cross-section
	IF (I .EQ. ZSPLI) THEN
	   CALL EXPMU(EN,XSEC(1))
C          ====================
	ENDIF
	XSEC_PE(I) = TATM(I)*XSEC(1)*0.1/CVOL
	XSEC_CO(I) = TATM(I)*XSEC(2)*0.1/CVOL
        XSEC_IN(I) = TATM(I)*XSEC(3)*0.1/CVOL
        XSEC_TO(I) = TATM(I)*XSEC(4)*0.1/CVOL
	
	XSCR_PE = XSCR_PE + XSEC_PE(I) 
        XSCR_CO = XSCR_CO + XSEC_CO(I) 
        XSCR_IN = XSCR_IN + XSEC_IN(I) 
        XSCR_TO = XSCR_TO + XSEC_TO(I)	
	ENDIF

220	CONTINUE

C	Density given as kg / m**3
	DENSI = MASS * 1E-3 * 1E30 / CVOL
C 	Density given as g / cm-3 = DENSI / 1000
	DENSIG = DENSI / 1000

C-------------------------------------------------------------------

        MU = XSCR_PE
C       Correction due to non-uniform beam
	CORREF=BEAMF(XBEA,YBEA,XFWHM,YFWHM)

	CALL DOSE(PHOSEC,MU,XCRY,YCRY,ZCRY,XBEA,YBEA,CORREF,DENSI,
     +		 IMAGES,EXPO,EN,PABS,GRYABS,HENDT,FRAC,FRACC,HENDLI)

C	WRITE (*,*) PHOSEC,MU,XCRY,YCRY,ZCRY,XBEA,YBEA,DENSI,
C     +		 	IMAGES,EXPO,EN,PABS,GRYABS,HENDT,FRAC,FRACC
C	WRITE (*,*) 'PHOSEC,MU,XCRY,YCRY,ZCRY,XBEA,YBEA,DENSI,IMAGES'
C	WRITE (*,*) 'EXPO,ENERGY,PABS,GRYABS,HENDT,FRAC,FRACC'
C	WRITE (*,*) 'ENERGY', EN

	CALL  TEMP(XCRY,YCRY,ZCRY,DENSI,PABS,CONVHE,CP,TSYS,DELTAT)

C Photons per unit cell
	PCELLH = HENDLI * DENSIG * CVOL *1E-27 / (EN * 1.602E-16)
	PCELLS = GRYABS * DENSIG * CVOL *1E-27 / (EN * 1.602E-16)
C Percentage absorbed
	PCABS = 100*(1-EXP(-XSCR_PE*ZCRY)) 
C       monomer concentration
	CONC = 1000 * NMON / (CVOL *1E-27 * NAVO)
	IMDOSE = GRYABS/IMAGES

C       OUTPUT 
	WRITE (*,2028) EN
	WRITE (*,2029) HENDLI
	WRITE (*,2027) CORREF
	WRITE (*,2030) CONC
	WRITE (*,2031) DENSIG
	WRITE (*,2032) SOLVEN * 100
	WRITE (*,2033) INT(NWAT)
	WRITE (*,2034) GRYABS
	WRITE (*,2035) IMDOSE
	WRITE (*,2047) XSCR_TO
	WRITE (*,2048) XSCR_PE
	WRITE (*,2049) PCABS
	WRITE (*,2036) HENDT
	WRITE (*,2037) INT(FRAC*100)
	WRITE (*,2038) INT(FRACC*100)
	WRITE (*,2039) CVOL
	WRITE (*,2041) PCELLS
	WRITE (*,2042) PCELLH
	WRITE (*,2045) TSYS
	WRITE (*,2046) DELTAT
	WRITE (*,*) ' '

2027	FORMAT ('Beam non-uniformity intensity correction ',F8.3)
2028	FORMAT ('Energy (keV)', F8.4)
2029	FORMAT ('Radiation Dose Limit (Henderson)', E8.1)
2030	FORMAT ('Protein Monomer Concentration (mM)    ', F8.1)
2031	FORMAT ('Crystal Density (g/cm^3)      ', F8.2)
2032	FORMAT ('Solvent Content (%)           ', F8.1) 
2033	FORMAT ('No. Water Molecules /cell     ', I8)
2034	FORMAT ('Dose in Grays                 ', E8.2)
2035	FORMAT ('Dose per image                ', E8.2)
2036	FORMAT ('Time to reach Henderson limit (s) ', F8.2)
2037	FORMAT ('Fraction of Beam Seen by the Crystal ', I3,'%')
2038	FORMAT ('Fraction of Crystal Seen by the Beam ', I3,'%')
2039	FORMAT ('Unit Cell Volume (A^3) ', F10.0)
2041	FORMAT ('Absorbed photons per unit cell per dataset ', F5.2)
2042	FORMAT ('Absorbed photons per unit cell for Henderson Limit ', F5.2)
2045	FORMAT ('t_sys System temperature time constant (s)', F8.2)
2046	FORMAT ('Delta T  (K)            ', F8.2)
2047	FORMAT ('Attenuation Coefficient (mm^01)',F8.3)
2048	FORMAT ('Absorption Coefficient (mm^01)',F8.3)
2049	FORMAT ('Absorption by the crystal (%)', F8.2)

C C       Table of atoms per unit cell
C 	WRITE (*,*) 'Attenuation Coefficient by Element and interaction'
C 	WRITE (*,*) ' Z   NAME  PROT   SOLV  TOTAL   XSCR_TO  XSCR_PE  XSCR_CO  XSCR_IN'
C 	DO 3000, I=1, ELEMEN
C 	IF ( TATM(I) .GT. 0) THEN 
C 	     CALL Z2NAME(NAME,i)
C 	     WRITE (*,2040) I,NAME, INT(PROT(I)),INT(SOLV(I)),INT(TATM(I)),XSEC_TO(I), XSEC_PE(I),XSEC_CO(I),XSEC_IN(I)
C 2040	     FORMAT (I3,'  ',A4,3I7,' ',4F9.6)

C 	     ENDIF
	    
C 3000	CONTINUE 
C 	WRITE (*,2050)  XSCR_TO, XSCR_PE, XSCR_CO, XSCR_IN
C 2050	FORMAT ('                         TOTAL ',4F9.6)



C--------------------------------------------
	IF (LRANGE .GT. 0) THEN
C write GRAPH file - free format columns of numbers.
C Iterate over energies

C	WRITE (*,2350) EN, HENDT, IMDOSE,GRYABS,PCELLH,PCELLS,
C     +     DIDOSE,DELTAT,XSCR_TO,XSCR_PE,PCABS
	WRITE (*,*) 'A: Energy (keV)'
	WRITE (*,*) 'B: Total exposure time to reach Henderson dose limit (s)'
	WRITE (*,*) 'C: Dose per image (Gy)'
	WRITE (*,*) 'D: Dose per Dataset (Gy)'
	WRITE (*,*) 'E: Photons absorbed per unit cell at the dose limit'
	WRITE (*,*) 'F: Photons absorbed per unit cell per dataset'
	WRITE (*,*) 'G: The diffracted intensity per absorbed dose (A.U.)'
	WRITE (*,*) 'H: The temperature rise (K)'
	WRITE (*,*) 'I: Attenuation coefficient (mm^-1) '
	WRITE (*,*) 'J: Absorption coefficient (mm^-1) '
	WRITE (*,*) 'K: Percentage beam absorption by the crystal'
	
	
	WRITE (*,*) '    A     B         C          D      E       F        G
     +    H       I        J      K'

C Recalculate Mu for each element	 
	   DO 2300 EN=LRANGE,URANGE,GSTEP
	   
	      MU = 0
	      XSCR_PE=0
	      XSCR_CO=0
	      XSCR_IN=0
	      XSCR_TO=0
	      DO 2150, I=1,ELEMEN
		 IF (TATM(I) .GT. 0) THEN
		    CALL Z2NAME(NAME,I)
		    CALL MUCAL(EN,NAME,I,UNIT,XSEC,EDGES,FLY, ERF,ER)
C                        ============================================	
	
C Empirical cross-section
		    IF (I .EQ. ZSPLI) THEN
		       CALL EXPMU(EN,XSEC(1))
C                      ====================
		    ENDIF
		    XSEC_PE(I) = TATM(I)*XSEC(1)*0.1/CVOL
		    XSEC_CO(I) = TATM(I)*XSEC(2)*0.1/CVOL
		    XSEC_IN(I) = TATM(I)*XSEC(3)*0.1/CVOL
		    XSEC_TO(I) = TATM(I)*XSEC(4)*0.1/CVOL
		    
		    XSCR_PE = XSCR_PE + XSEC_PE(I) 
		    XSCR_CO = XSCR_CO + XSEC_CO(I) 
		    XSCR_IN = XSCR_IN + XSEC_IN(I) 
		    XSCR_TO = XSCR_TO + XSEC_TO(I)	
		 ENDIF

 2150	      CONTINUE

	      MU = XSCR_PE
C       Correction due to non-uniform beam
	      CORREF=BEAMF(XBEA,YBEA,XFWHM,YFWHM)

	CALL DOSE(PHOSEC,MU,XCRY,YCRY,ZCRY,XBEA,YBEA,CORREF,DENSI,
     +		 IMAGES,EXPO,EN,PABS,GRYABS,HENDT,FRAC,FRACC,HENDLI)

	CALL  TEMP(XCRY,YCRY,ZCRY,DENSI,PABS,CONVHE,CP,TSYS,DELTAT)

C Photons per unit cell
	PCELLH = HENDLI * DENSIG * CVOL *1E-27 / (EN * 1.602E-16)
	PCELLS = GRYABS * DENSIG * CVOL *1E-27 / (EN * 1.602E-16)
C Diffracted Intensity per absorbed dose unit. (Arbitraty units)
	DIDOSE = 100*(1/EN)**3 * EXP(-XSCR_TO*ZCRY)/(1-EXP(-XSCR_PE*ZCRY))
C Percent absorption
	PCABS = 100*(1- EXP(-MU*ZCRY)) 
C       monomer concentration
	CONC = 1000 * NMON / (CVOL *1E-27 * NAVO)
	IMDOSE =  GRYABS / IMAGES

	WRITE (*,2350) EN, HENDT, IMDOSE,GRYABS,PCELLH,PCELLS,
     +     DIDOSE,DELTAT,XSCR_TO,XSCR_PE,PCABS

2350	 FORMAT(F8.3,' ',F8.2,' ',E8.2,' ',E8.2,' ',F6.2,' ',F6.2,' ',E8.3,' ',
     +       F5.1,' ',F8.3,' ',F8.3,' ',F8.2)
	
 2300	   CONTINUE                                              
       

	ENDIF



	CALL CCPERR(0, 'Normal termination')
	END

C - END of the main program



C-----------------------------------------------------------------------
	REAL FUNCTION CELVOL(CELL)
C-----------------------------------------------------------------------
	IMPLICIT NONE 
	REAL VOLUME, CELL(6), ALPHA,BETA,GAMMA
	REAL ULT,TMP
	REAL PI
	PARAMETER (PI=3.1415926536)

	ALPHA = CELL(4)*PI/180.0
	BETA = CELL(5)*PI/180.0
	GAMMA = CELL(6)*PI/180.0
	   ULT = 1.0 + 2.0*COS(ALPHA)*COS(BETA)*COS(GAMMA) -
     + 	         COS(ALPHA)**2 - COS(BETA)**2 - COS(GAMMA)**2 
	   IF (ULT .LE. 0.0) CALL CCPERR(1, 
     +'** The cell volume cannot be calculated, check the CELL card **')
	   CELVOL = CELL(1) * CELL(2) * CELL(3) * SQRT(ULT)
	   RETURN
	   END


C----------------------------------------------------------------
	SUBROUTINE TEMP(XCRY,YCRY,ZCRY,DENSI,PABS,CONVHE,CP,TSYS,DELTAT)
C----------------------------------------------------------------
C Convective heat transfer from Kuzay ActaD 57 (2001) 69--81
C accurate? - "Snowflake in hell" model.
	IMPLICIT NONE

	EXTERNAL CCPERR
	
C convective heat transfer coeff. W m-2 K-1
	REAL CONVHE
C heat capacity J / kg
	REAL CP

	REAL ASURF, TSYS, DELTAT, PABS
	REAL DENSI, CVOL

	REAL XCRY, YCRY, ZCRY
	REAL XCRYT,YCRYT,ZCRYT
	XCRYT = XCRY / 1000
	YCRYT = YCRY / 1000
	ZCRYT = ZCRY / 1000
	CVOL = 	XCRYT*YCRYT*ZCRYT
	ASURF= 2*( XCRYT*YCRYT + XCRYT*ZCRYT + YCRYT*ZCRYT)
	TSYS = (CP*DENSI*CVOL)/(CONVHE*ASURF)
	DELTAT = PABS/(CONVHE*ASURF)
	
	END
C------------------------------------------------------------------

C Read CHOOCH SPLINOR FILE
C------------------------------------------------------------------
	SUBROUTINE RSPLIN
C------------------------------------------------------------------
	
	IMPLICIT NONE 
	COMMON /SPLINE/SPLIF,ESPLI,ZSPLI,SPLI,NSPLI,MUSPL
	REAL SPLI(2,1000),MUSPL(2)
	CHARACTER SPLIF*132
	CHARACTER ESPLI*2
	INTEGER ZSPLI
	INTEGER NSPLI

C	COMMON /MUC/XSEC,EDGES,FLUOR,ERF
	REAL XSEC(10),EDGES(9),FLUOR(4)
	LOGICAL ERF
	INTEGER ER

	EXTERNAL  CCPERR,CCPFYP,CCPRCS,MEMOPARSE, NAME2Z, Z2NAME,
     +            CCPDPN,QCLOSE,LENSTR,MUCAL

	INTEGER LFLAG
	REAL JUNK1,JUNK2,JUNK3
	INTEGER NPOINT,I
	PARAMETER (NPOINT=1000)
C       Splinor file has max 1000 points

C Read SPLINOR file
	CALL CCPDPN(3,SPLIF,'READONLY','F',0,LFLAG)
	IF (LFLAG .LT. 0) CALL CCPERR(1, 'ERROR OPENING SPLINOR FILE') 
	NSPLI=1
4999	READ(3,5010,END=5000) SPLI(1,NSPLI),SPLI(2,NSPLI),JUNK1,JUNK2,JUNK3
	NSPLI=NSPLI+1

	GO TO 4999
5000	CONTINUE
	NSPLI=NSPLI-1

C Convert to keV
	DO 5050 I=1,NSPLI,1
	   SPLI(1,I) = SPLI(1,I)/1000.0
5050	CONTINUE

C Normalise f'' to PE crosssection
	CALL NAME2Z(ESPLI,ZSPLI)
	JUNK1=SPLI(1,1)
	CALL MUCAL(SPLI(1,1),ESPLI,ZSPLI,'B',XSEC,EDGES,FLUOR,ERF,ER)
	MUSPL(1) = XSEC(1)
	CALL MUCAL(SPLI(1,NSPLI),ESPLI,ZSPLI,'B',XSEC,EDGES,FLUOR,ERF,ER)
	MUSPL(2)=XSEC(1)
C	WRITE (*,*) SPLI(1,1), MUSPL(1), SPLI(1,NSPLI), MUSPL(2)
C	DO 5020 I=1,NSPLI  ,1
C 	   WRITE (*,*) I,SPLI(1,I), SPLI(2,I)
C 5020	CONTINUE
C	WRITE (*,*) NSPLI

5010	FORMAT(5F13.3)
	
	RETURN

	END

C------------------------------------------------------------------
	SUBROUTINE EXPMU(EN,MU)
C------------------------------------------------------------------
	IMPLICIT NONE 
	COMMON /SPLINE/SPLIF,ESPLI,ZSPLI,SPLI,NSPLI,MUSPL
	REAL SPLI(2,1000), MUSPL(2)
	CHARACTER SPLIF*132
	CHARACTER ESPLI*2
	INTEGER ZSPLI,NSPLI

	COMMON /MUC/XSEC,EDGES,FLUOR
	REAL XSEC(10),EDGES(9),FLUOR(4)
	LOGICAL ERF
	INTEGER ER

	REAL EN, MU
	INTEGER I
	REAL FRAC,NORM
	REAL MUNEW

C Is energy in range of scan?
	IF ( (EN .LE. SPLI(1,1)) .OR. (EN .GE. SPLI(1,NSPLI))) THEN
	   RETURN
	ENDIF

	DO 6000 I=1,NSPLI,1
	   IF (  (EN .GT. SPLI(1,I)) .AND. (EN .LE. SPLI(1,I+1)) ) THEN
	      FRAC = (EN - SPLI(1,I)) / (SPLI(1,I+1) - SPLI(1,I))
	      NORM=FRAC*(SPLI(2,I+1)-SPLI(2,I))+SPLI(2,I)
	      MUNEW=NORM*( MUSPL(2)-MUSPL(1)) + MUSPL(1)
	      MU=MUNEW
	      RETURN
	   ENDIF
 6000	CONTINUE
	
	WRITE (*,*) 'ERROR in EXPMU routine'
	RETURN
	END

C------------------------------------------------------------------
	SUBROUTINE DOSE(PHOSEC,MU,XCRY,YCRY,ZCRY,XBEA,YBEA,CORREF,DENSI,
     *     IMAGES,EXPO,EN,PABS,GRYABS,HENDTI,FRAC,FRACC,HENDLI)
C------------------------------------------------------------------

	IMPLICIT NONE
	
	EXTERNAL CCPERR

	REAL PHOSEC,MU,XCRY,YCRY,ZCRY,XBEA,YBEA,DENS,EXPO,EN,GRYABS
	REAL PABS,CRYABS,CORREF
	REAL DENSI, HENDTI
	INTEGER IMAGES
	REAL TEST

	REAL PCELLH,PCELLS

	REAL FRAC, FRACC, PHOSET, PHOABS, JOUABS, VOLUME
	REAL JPKEV, HENDLI
	PARAMETER (JPKEV = 1000 * 1.602E-19)

C Calculate photons per dataset
	PHOSET = EXPO * IMAGES * PHOSEC	* CORREF

C Calculate relative crystal absorbance
	CRYABS = ( 1.0 - EXP(-1.0*MU*ZCRY) )

C Calculate fraction of beam incident on crystal.
	IF ((XBEA .LE. XCRY) .AND. (YBEA .LE. YCRY))
     +		FRAC = 1.0
	IF ((XBEA .GT. XCRY) .AND. (YBEA .LE. YCRY))
     +		FRAC = XCRY/XBEA
	IF ((XBEA .LE. XCRY) .AND. (YBEA .GT. YCRY))
     +		FRAC = YCRY/YBEA
	IF ((XBEA .GT. XCRY) .AND. (YBEA .GT. YCRY))
     +		FRAC = (XCRY/XBEA)*(YCRY/YBEA)

	IF (FRAC .GT. 1) THEN
		CALL CCPERR(1,'** beam/crystal parameters are invalid **')
	ENDIF

C Calculate photons absorbed
	PHOABS = FRAC * PHOSET * CRYABS

C Calculate joules absorbed
	JOUABS = PHOABS * EN * JPKEV

C Calculate power absorbed
	PABS = JOUABS / (IMAGES * EXPO )

C Calculate aborbed dose in Gray
C Need volume of xtal illuminated VOLUME in m**3
	IF ((XBEA .LE. XCRY) .AND. (YBEA .LE. YCRY))
     +		VOLUME = ZCRY * XBEA * YBEA * 1E-9
	IF ((XBEA .GT. XCRY) .AND. (YBEA .LE. YCRY))
     +         VOLUME = ZCRY * XCRY * YBEA * 1E-9
	IF ((XBEA .LE. XCRY) .AND. (YBEA .GT. YCRY))
     +         VOLUME = ZCRY * XBEA * YCRY * 1E-9
	IF ((XBEA .GT. XCRY) .AND. (YBEA .GT. YCRY))
     +         VOLUME = ZCRY * XCRY * YCRY * 1E-9
	FRACC =  1E9 * VOLUME / (XCRY*YCRY*ZCRY)

C Calculate absorbed dose in Gray and time taken to reach henderson
C limit
	GRYABS = JOUABS / ( DENSI * VOLUME )
	HENDTI = HENDLI * EXPO * REAL(IMAGES) / GRYABS


	RETURN

	END

C------------------------------------------------------------------
        REAL FUNCTION BEAMF(XBEA,YBEA,XFWHM,YFWHM)
C------------------------------------------------------------------
	IMPLICIT NONE

 	EXTERNAL CCPERR
C FWHM of a Gaussian is 2 sqrt ( 2 ln 2) = 2.35482 sigma	
	REAL SIGFW
	REAL SIGX,SIGY,C,XFWHM,YFWHM,X,Y,XBEA,YBEA
	INTEGER GRID,A,B
	PARAMETER (GRID=21)
	PARAMETER (SIGFW= 2.35482)
	REAL PROF(21,21)
	INTEGER CENTRE
	
	CENTRE=GRID/2 
	SIGX = XFWHM/SIGFW
	SIGY = YFWHM/SIGFW
C
	IF (XFWHM .LT. 0) THEN 
	   BEAMF=1.0
	   RETURN
	ENDIF
C
	C=0
	DO 4000 A=1,GRID,1
	   DO 4010 B=1,GRID,1
 	      X = (XBEA/GRID)* (A-CENTRE)
 	      Y = (YBEA/GRID)* (B-CENTRE)
 	      C = C + EXP( -(X*X)/(2*SIGX**2) - (Y*Y)/(2*SIGY**2) )
 	      PROF(A,B) = EXP( - (X*X)/(2*SIGX**2) - (Y*Y)/(2*SIGY**2) )
4010	   CONTINUE
4000	CONTINUE

 	DO 4020 A=1,GRID,1
 	   DO 4030 B=1,GRID,1
 	      PROF(A,B) = PROF(A,B)/ C
4030	   CONTINUE
4020	CONTINUE

	BEAMF = GRID*GRID/C
	
	RETURN
	END

