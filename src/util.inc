#INLINE F90_UTIL

SUBROUTINE InitVal (concs)

USE model_Global
use model_Monitor
INTEGER :: i,q,found
REAL(kind=dp) :: x
CHARACTER(LEN=10000) :: line
CHARACTER(LEN=15) :: test
REAL(kind=dp) :: concs(NSPEC)
REAL(kind=dp) :: oldvar(NSPEC)
INTEGER :: COUNT_NOX_CONSTRAINTS, i_EMISS
CHARACTER(len=6), dimension(8) :: NOXARRAY
NOXARRAY(:)= (/'NO2   ','NO    ','NO3   ','HONO  ','HNO2  ',&
    &'HO2NO2','HNO4  ','PNA   '/)
    OLDVAR(1:NVAR)=C(1:NVAR)
    VAR(1:NVAR)=0.
    C(1:NVAR)=0.
    i_EMISS=0
    DEPOS = 0.
    SAREA=0.
    ALBEDO=0.
    RP1=0.

    do I=1,NSPEC+1
    FOUND=0
    if (SPEC_NAME(I) .NE. '') then
    SPEC_NAME(I)=ADJUSTL(SPEC_NAME(I))
    do Q=1,NVAR

        if (TRIM(SPEC_NAME(I)) .EQ. TRIM(SPC_NAMES(Q))) then
            FOUND=1
            VAR(Q)=CONCS(I)


            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            print*,  achar(27)//'[95m' , spec_name(i),&
               & concs(i),achar(27)//'[97m'
            
            



            if (const_method(i) > 0) then
                 CONSTRAIN(Q)=CONCS(I)
                 print*, 'soft-constrained ', spec_name(i), 'at', constrain(Q)
            end if

            !populate spin up array with any concentrations in ics
            spin_const(Q)=concs(i)


            if (trim(SPEC_NAME(I)) .EQ. 'CH4') then
                SPEC_CH4=.TRUE.


                else if (trim(SPEC_NAME(I)) .EQ. 'H2') then
                SPEC_H2=.TRUE.

            end if

        endif

    enddo




    if (FOUND .ne. 1) then

    SELECT CASE (TRIM(SPEC_NAME(I)))

    CASE('H2O')
    H2O=CONCS(I)

    CASE('PRESS')
    PRESS=CONCS(I)

    CASE('NOx')
        if (concs(i)>0) then
        CONSTRAIN_NOX=.TRUE.
        WRITE (OUTPUT_UNIT,*) 'Constraining total NOx concentation'
        end if

    CASE('LAT')
    LAT=CONCS(I)

    CASE('LON')
    LON=CONCS(I)

    CASE('TEMP')
    TEMP=CONCS(I)

    CASE('JDAY')
    JDAY=CONCS(I)

    CASE('TIME(GMTs)')
    continue

    CASE('O3COL')
    O3COL=CONCS(I)

    CASE('ALBEDO')
    ALBEDO=CONCS(I)

    CASE('SAREA')
    SAREA=CONCS(I)

    CASE('RP1')
    RP1=CONCS(I)

    CASE('JNO2')
        if (CONST_METHOD(I) .GE. 1) JNO2=CONCS(I)

        JREPEAT=0

        if (CONST_METHOD(I) .EQ. 2) JREPEAT=1


    CASE('JO1D')

        if (CONST_METHOD(I) .GE. 1) JO1D=CONCS(I)

        JREPEAT=0

        if (CONST_METHOD(I) .EQ. 2) JREPEAT=1


    CASE('DEPOS')

        if (CONCS(I) > 0) then
            DEPOS = 1.
            print *, 'DEPOSITON ENABLED'
        end if

    CASE('SPINUP')
        spinup = CONCS(I)

    CASE DEFAULT

            !! for any species in the fixed array
            if (SPEC_NAME(I)(1:1) .EQ. 'F')  then
                SPEC_NAME(i)  = SPEC_NAME(i)(2:)


                do Q=NVAR+1,NVAR+NFIX

                    if (TRIM(SPEC_NAME(I)) .EQ. TRIM(SPC_NAMES(Q))) then
                        FOUND=1
                        FIX(Q-NVAR)=CONCS(I)
                        print*, 'hard-constrained ', spec_name(i),&
                            &' at ', concs(I)

                        if (trim(SPEC_NAME(I)) .EQ. 'EMISS') then
                            i_EMISS = (Q-NVAR)
                        else if (trim(SPEC_NAME(I)) .EQ. 'CH4') then
                            SPEC_CH4=.TRUE.
                        else if (trim(SPEC_NAME(I)) .EQ. 'H2') then
                            SPEC_H2=.TRUE.         
                        end if            
            
                        exit 
                   else 
                     WRITE(OUTPUT_UNIT,*) SPEC_NAME(I),' NOT FOUND - stopping'
                     WRITE(6,*) SPEC_NAME(I),' NOT FOUND (FIXED) - stopping'
                     stop
                    endif
                    
                  
                    
                enddo

            end if

            if (FOUND .ne. 1 .and. SPEC_NAME(I)(1:1) .NE. 'X') then
                WRITE(OUTPUT_UNIT,*) 'spec_name:',SPEC_NAME(I),&
                        &' NOT FOUND - stopping'
                WRITE(6,*) SPEC_NAME(I),' NOT FOUND - stopping '
                stop
            end if
            WRITE (OUTPUT_UNIT,*) 'Starts with X => ignored'


    END SELECT
    endif
    endif


    enddo


    !convert to mol cm-3s-1
    CFACTOR=PRESS*1e2*1e-6/(8.314*TEMP)*6.022E23
    H2O=H2O*CFACTOR
    VAR(:)=VAR(:)*CFACTOR
    CONSTRAIN(:)=CONSTRAIN(:)*CFACTOR


    spin_const(:)=spin_const(:) * CFACTOR


    FIX(:)=FIX(:)*CFACTOR
    if (i_EMISS>0) FIX(i_EMISS)= FIX(i_EMISS)/CFACTOR
    !correct emissions ^^


    if (CONSTRAIN_NOX) then
        COUNT_NOX_CONSTRAINTS=0


        do I=1,NVAR
        do Q= 1,8

        if (TRIM(SPC_NAMES(I)) .EQ. TRIM(NOXARRAY(Q))) then
        NOX(I)=1
        CONSTRAIN(I)=-1.*CONSTRAIN(I)
        if (CONSTRAIN(I) .NE. 0 .and.&
            TRIM(SPC_NAMES(I)).ne. 'PNA')&
            COUNT_NOX_CONSTRAINTS=COUNT_NOX_CONSTRAINTS+1
            WRITE (OUTPUT_UNIT,*) SPC_NAMES(I),' IN NOX FAMILY'
            exit
            end if
            end do

            if (TRIM(SPC_NAMES(I)) .EQ. 'N2O5') then
            NOX(I)=2
            CONSTRAIN(I)=-1.*CONSTRAIN(I)
            if (CONSTRAIN(I) .NE. 0)  COUNT_NOX_CONSTRAINTS=&
            COUNT_NOX_CONSTRAINTS+1
            endif

            end do

            if (COUNT_NOX_CONSTRAINTS .GT. 1) then
            WRITE (OUTPUT_UNIT,*) 'You can not contrain multiple NOX species'
            STOP
            endif

    end if


        ! FIND CH4 and H2 species
        if (SPEC_CH4 .EQv. .FALSE.) then
        WRITE (OUTPUT_UNIT,*) 'No CH4 specified assuming 1770 ppbv'

     print*, 'fix no ch4 default, add in model subset?,&
     &add generic reaction of inorganics?'
        !if (ind_CH4 .gt. NVAR) then
        !    FIX((ind_CH4-NVAR)) = 1770e-9*CFACTOR
        !else
        !    VAR(ind_CH4)=1770e-9*CFACTOR
        !    CONSTRAIN(ind_CH4)=VAR(ind_CH4)
        !end if
        !end if

        !if (SPEC_H2 .EQv. .FALSE.) then
        WRITE (OUTPUT_UNIT,*) 'No H2 specified assuming 550 ppbv'
          if (ind_H2 .gt. NVAR) then 
            !FIX((ind_H2-NVAR)) = 550e-9*CFACTOR
        else 
            VAR(ind_H2)= 550e-9*CFACTOR
            CONSTRAIN(ind_H2)=VAR(ind_H2)
        end if
        end if

        !set R value to 0 
        FIX(indf_R)= 1

        END SUBROUTINE Initval
        ! End of InitVal function
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


        ! ****************************************************************
        ! ****************************************************************
        !
        ! SaveData - Write LOOKAT species in the data file
        !   Parameters :  index_loc = daycounter
        !
        ! ****************************************************************

        subroutine SaveOut(index_loc)
        USE model_Global
        USE model_Monitor
        USE model_Integrator
        integer :: index_loc

        !print *, 'adjust to reflect multiple times'
        !spec
        output_s(index_loc,1)=Jday*86400  + DAYCOUNTER*dt
        output_s(index_loc,2)=Lat
        output_s(index_loc,3)=Lon
        output_s(index_loc,4)=Press
        output_s(index_loc,5)=Temp
        output_s(index_loc,6)=H2O
        output_s(index_loc,7)=Cfactor
        output_s(index_loc,8)=JfactNO2
        output_s(index_loc,9)=JfactO1D
        output_s(index_loc,10)=RO2
        output_s(index_loc,11:) = C
        !rate
        output_r(index_loc,1)=Jday*86400  + DAYCOUNTER*dt
        output_r(index_loc,2)=Lat
        output_r(index_loc,3)=Lon
        output_r(index_loc,4)=Press
        output_r(index_loc,5)=Temp
        output_r(index_loc,6)=H2O
        output_r(index_loc,7:) = Rconst

        end subroutine SaveOut


        ! ****************************************************************


        SUBROUTINE ICNVRT(WAY,NUM,STRING,LENGTH,IERR)
        !
        !       FUNCTION:
        !F
        !F        This subroutine does an integer-to-character conversion
        !F        or a characater-to-integer conversion depending on the
        !F        integer WAY:
        !F                If WAY = 0 then an integer-to-character conversion
        !F                is done. If WAY .NE. 0 then a character-to-integer
        !F                conversion is done.
        !F
        !       USAGE:
        !U
        !U        !ALL I!NVRT(WAY,NUM,STRING)[
        !U             where WAY, NUM, STRING, and LENGTH are defined below.
        !U
        !U        Example: !ALL I!NVRT(0,1000,STRING,LENGTH)
        !U                 on return STRING = '1000' and
        !U                 LENGTH = 4.
        !U
        !       INPUTS:
        !I
        !I        WAY - INTEGER; Determines which way the conversion goes:
        !I              if WAY = 0 then an integer-to-character conversion
        !I                         is performed;
        !I              if WAY.NE.0 then a character-to-integer conversion
        !I                         is performed.
        !I
        !I         NUM - INTEGER; an input only if WAY = 0. NUM is the integer
        !I               number to be converted to a character expression.
        !I
        !I         STRING - !HARACTER; an input only if WAY .NE. 0. STRING
        !I                is the character expression to be converted to an
        !CI                integer value. It contain no decimal points or
        !I                non-numeric characters other than possibly a
        !I                sign. If STRING contains  a '+' sign, it will be
        !I                stripped of it on return.
        !I
        !       OUTPUTS:
        !O
        !O         NUM - INTEGER; contains the INTEGER representation of
        !O                STRING.
        !O
        !O         STRING - !HARA!TER; contains the CHARACTER representation of
        !O                  NUM.
        !CO
        !CO         LENGTH - INTEGER; The length of STRING to the first blank.
        !CO                  The significant part of STRING can be accessed with
        !CO                  the declaration STRING(1:LENGTH).
        !CO
        !CO         IERR - INTEGER variable giving return condition:
        !CO                IERR = 0 for normal return;
        !CO                IERR = 1 if NUM cannot be converted to STRING because
        !CO                       STRING is too short or STRING cannot be
        !CO                       converted to NUM because STRING is too long.
        !CO                IERR = 2 if STRING contained a non-numeric character
        !CO                       other than a leading sign or something went
        !CO                       wrong with an integer-to-character conversion.
        !CO
        !C       ALGORITHM:
        !CA
        !CA         Nothing noteworthy, except that this subroutine will work
        !CA          for strange character sets where the character '1' doesn't
        !CA          follow '0', etc.
        !CA
        !C       MACHINE DEPENDENCIES:
        !CM
        !CM          The parameter MAXINT (below) should be set to the
        !CM          number of digits that an INTEGER data type can have
        !CM          not including leading signs. For VAX FORTRAN V4.4-177
        !CM          MAXINT = 10.
        !CM
        !CM          NOTE: Under VAX FORTRAN V4.4-177, the
        !CM          error condition IERR = 1 will never occur for an
        !CM          integer-to-character conversion if STRING
        !CM          is allocated at least 11 bytes (CHARACTER*11).
        !CM
        !C       HISTORY:
        !CH
        !CH      written by:             bobby bodenheimer
        !CH      date:                   september 1986
        !CH      current version:        1.0
        !CH      modifications:          none
        !CH
        !C       ROUTINES CALLED:
        !CC
        !CC          NONE.
        !CC
        !C----------------------------------------------------------------------
        !C       written for:    The CASCADE Project
        !C                       Oak Ridge National Laboratory
        !C                       U.S. Department of Energy
        !C                       contract number DE-AC05-840R21400
        !C                       subcontract number 37B-7685 S13
        !C                       organization:  The University of Tennessee
        !C----------------------------------------------------------------------
        !C       THIS SOFTWARE IS IN THE PUBLIC DOMAIN
        !C       NO RESTRICTIONS ON ITS USE ARE IMPLIED
        !C----------------------------------------------------------------------
        !C
        !C Global Variables.
        !C
        INTEGER       WAY
        INTEGER       LENGTH
        INTEGER       NUM
        INTEGER       IERR
        CHARACTER*(*) STRING
        !C
        !C Local Variables
        !C
        INTEGER       I
        INTEGER       MAXINT
        INTEGER       MNUM
        INTEGER       M
        LOGICAL       NEG
        character(1) dummy
        integer    dummyvar
        !C
        PARAMETER(MAXINT=10)
        !C
        NEG = .FALSE.
        IERR = 0
        !C
        !C  Integer-to-character conversion.
        !C
        IF (WAY.EQ.0) THEN
        STRING = ' '
        IF (NUM.LT.0) THEN
        NEG = .TRUE.
        MNUM = -NUM
        LENGTH = INT(LOG10(REAL(MNUM))) + 1
        ELSE IF (NUM.EQ.0) THEN
        MNUM = NUM
        LENGTH = 1
        ELSE
        MNUM = NUM
        LENGTH = INT(LOG10(REAL(MNUM))) + 1
        END IF
        IF (LENGTH.GT.LEN(STRING)) THEN
        IERR = 1
        RETURN
        END IF


        DO 10, I=LENGTH,1,-1
        M=INT(REAL(MNUM)/10**(I-1))

        !jkhjkhjkhjkhljhlkjhj
        if (m.lt.10) then
        Write( dummy, '(i1)' ) M
        STRING(LENGTH-I+1:LENGTH-I+1) = dummy
        ELSE
        IERR = 2
        RETURN
        END IF


        MNUM = MNUM - M*10**(I-1)
        10       CONTINUE
        IF (NEG) THEN
        STRING = '-'//STRING
        LENGTH = LENGTH + 1
        END IF
        !C
        !C  Character-to-integer conversion.
        !!C
        ELSE
        IF (STRING(1:1).EQ.'-') THEN
        NEG = .TRUE.
        STRING = STRING(2:LEN(STRING))
        END IF
        IF (STRING(1:1).EQ.'+') STRING = STRING(2:LEN(STRING))
        NUM = 0
        LENGTH = INDEX(STRING,' ') - 1
        IF (LENGTH.GT.MAXINT) THEN
        IERR = 1
        RETURN
        END IF


        DO 20, I=LENGTH,1,-1
        read( STRING(LENGTH-I+1:LENGTH-I+1), '(i10)' ) dummyvar
        if (dummyvar.lt.10) then
        M= (dummyvar)
        ELSE
        IERR = 2
        RETURN
        END IF

        NUM = NUM + INT(10**(I-1))*M
        20       CONTINUE

        IF (NEG) THEN
        NUM = -NUM
        STRING = '-'//STRING
        LENGTH = LENGTH + 1
        END IF
        END IF
        !C
        !C  Last lines of ICNVRT
        !C
        RETURN
        END SUBROUTINE



        #ENDINLINE
