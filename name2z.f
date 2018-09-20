c-----------------------------------------------------------------------                  
      subroutine name2Z(name,Z)
      implicit none
      character name*2,elements(105)*2
      integer Z,Znumbers(105),j
      data elements/'H ','HE','LI','BE','B ','C ','N ','O ','F ','NE',
     *              'NA','MG','AL','SI','P ','S ','CL','AR',
     *              'K ','CA','SC','TI','V ','CR','MN','FE','CO','NI',
     *              'CU','ZN','GA','GE','AS','SE','BR','KR',
     *              'RB','SR','Y ','ZR','NB','MO','TC','RU','RH','PB',
     *              'AG','CD','IN','SN','SB','TE','I ','XE',
     *              'CS','BA','LA','CE','PR','ND','PM','SM','EU','GD',
     *              'TB','DY','HO','ER','TM','YB','LU','HF','TA','W ',
     *              'RE','OS','IR','PT','AU','HG','TL','PB','BI','PO',
     *              'AT','RN','FR','RA','AC','TH','PA','U ','NP','PU',
     *              'AM','CM','BK','CF','ES','FM','MD','NO','LR','RF',
     *              'HA'/
      data Znumbers/   1,   2,   3,   4,   5,   6,   7,   8,   9,  10,
     *                11,  12,  13,  14,  15,  16,  17,  18,
     *                19,  20,  21,  22,  23,  24,  25,  26,  27,  28,
     *                29,  30,  31,  32,  33,  34,  35,  36,
     *                37,  38,  39,  40,  41,  42,  43,  44,  45,  46,
     *                47,  48,  49,  50,  51,  52,  53,  54,
     *                55,  56,  57,  58,  59,  60,  61,  62,  63,  64,
     *                65,  66,  67,  68,  69,  70,  71,  72,  73,  74,
     *                75,  76,  77,  78,  79,  80,  81,  82,  83,  84,
     *                85,  86,  87,  88,  89,  90,  91,  92,  93,  94,
     *                95,  96,  97,  98,  99, 100, 101, 102, 103, 104,
     *               105/
      call upcase(name)
      do j=1,105
        if (name(1:2).eq.elements(j)(1:2)) then
          z=Znumbers(j)
          go to 100
        endif
      enddo
c-- Error if the prog gets here                                                 
      write(*,1000) name
 1000  format('ERROR: Element ',a2,' is not recognised.')
      stop
 100   continue
      return

      end
