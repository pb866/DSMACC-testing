PROGRAM driver
USE model_global
USE model_Parameters
USE model_Rates,       ONLY: Update_SUN, Update_RCONST
USE model_integrator,  ONLY: integrate
USE model_monitor,     ONLY: spc_names,Eqn_names
USE model_Util
USE constants


IMPLICIT NONE
REAL(dp) :: ENDSTATE(NVAR), total, RATIO, TNOX, TNOX_OLD
REAL(dp) :: STARTSTATE(NVAR), TIMESCALE, RH, RSTATE(20)
REAL(dp) :: DIURNAL_OLD(NVAR,3000), DIURNAL_NEW(NVAR,3000)
REAL(dp) :: DIURNAL_RATES(NREACT, 3000)
REAL(dp) :: FULL_CONCS(NSPEC,999999), concs(NSPEC)
! Photolysis calculation variables
character(len=50), allocatable ::  s_names(:), r_names(:)


! DO NOT NEED ABOVE

REAL(dp) :: NOXRATIO,Alta,Fracdiff,SpeedRatio,oldfracdiff,FRACCOUNT, newtime
character(50) :: counter, cw,filename
character (3) :: ln
INTEGER  :: ERROR, IJ, PE ,runtimestep,ICNTRL_U(20)
Integer  :: CONSTNOXSPEC, JK, full_counter, line, nc_set, nc_counter
character(200) :: dummychar
integer :: run_counter = 0

!!!~~~~ BEGIN TEMPERTURE/CONCENTRATION CONSTRAINTS Edwards et al. (2014) ~~~~!!!
! find the correct half hour of the model integration
INTEGER :: HALFHOUR
REAL    :: OBS_TEMP(48)=(/&
     -4.90988,-5.998,-6.76708,-8.16098,-8.60077,-8.97174,-9.09161,-9.34328,-9.2665,-9.42021,&
     -9.32238,-9.36907,-9.49675,-9.27415,-9.07073,-9.18774,-10.3613,-10.3786,-10.3822,-10.5816,&
     -10.5055,-10.5738,-10.9845,-11.2199,-11.5564,-11.5178,-11.4663,-11.4704,-11.4965,-11.5318,&
     -11.6473,-11.7589,-11.2993,-10.784,-10.0277,-9.14637,-7.99273,-7.53081,-6.87639,-6.05381,&
     -5.63831,-5.54391,-4.88845,-4.97976,-4.86069,-4.50494,-4.76706,-4.63537/)



! 2012 HONO
REAL     :: HONO(48)=(/&
     28.88,27.85,24.45,24.59,29.24,28.68,30.13,31.81,32.57,&
     41.21,43.53,44.03,47.46,47.28,49.,54.26,59.4,58.79,58.81,&
     59.1, 61.15,60.21,62.51,63.9, 62.79, 62.91, 67.86,68.44, 69.7,&
      74.19, 76.06 ,92.93, 92.57, 86.27, 68.17, 45.34, 37.58, 30.85, 29.,&
      26.8, 30.18, 38.34, 33.73, 33.6, 29.07, 26.01, 29.28,32.42/)

REAL     :: CLNO2(48)=(/&
     21.0994,26.3102,45.4727,78.7776,153.944,191.638,233.375,272.866,275.795,267.478,&
     272.051,282.215,259.541,230.965,235.998,259.276,301.291,297.618,300.514,340.369,&
     301.005,279.349,261.297,279.124,296.919,291.245,287.455,276.348,279.633,292.379,&
     244.924,192.209,160.123,119.415,93.7459,59.6054,39.4277,26.4976,19.0286,15.7525,&
     13.7339,13.5311,14.6732,15.9641,15.8982,15.7743,15.1901,16.8845/)

! Have called C10 aromatics 5-ethyl-m-xylene (DIME35EB) as this is the only C10 in MCMv3.2
REAL      ::  DIME35EB(48)=(/&
     69.452,70.3056,67.4846,74.1145,74.8584,73.826,82.6164,88.9278,76.7496,83.3863,&
     75.3736,67.7071,67.4484,64.812,64.2898,67.9836,70.5326,65.4262,69.1544,76.5433,&
     85.9401,75.6455,79.2941,75.1509,82.3132,80.2801,75.915,87.1978,81.9256,81.6976,&
     104.801,104.91,100.605,108.964,101.727,98.7474,102.016,103.402,94.3561,84.8273,&
     74.6906,80.2924,81.3013,77.1724,76.3041,74.0445,64.804,66.351/)

! Have grouped C11 and C12 aromatics as 3,5-diethyl toluene (DIET35TOL) as this is the only C10+ in MCMv3.2
REAL      ::  DIET35TOL(48)=(/&
     36.8734,37.0257,34.21,35.1549,37.7419,34.8994,37.9065,40.8461,36.5408,37.1806,&
     34.3297,32.7743,32.0106,31.163,29.9589,32.8109,32.9365,32.3456,34.3427,36.1953,&
     38.1519,35.3669,34.9246,34.7808,35.198,35.4073,34.135,38.581,36.941,35.414,&
     44.537,44.8959,45.9063,47.5709,46.9998,46.6037,48.8989,48.9082,46.6074,44.5145,&
     41.5593,44.9039,46.0723,43.1969,40.4023,39.1257,34.5272,35.63/)

!!!~~~~ END CONCENTRATION CONSTRAINTS Edwards et al. (2014) ~~~~!!!


STEPMIN = 0.0_dp
STEPMAX = 0.0_dp
RTOL(:) = 1.0E-5_dp !-5
ATOL(:) = 1.0_dp    !-16?

!desired |true-computed| < RTOL*|TRUE| + ATOL
!want ATOL/calc_value < RTOL

!rtol - #sig fig

mechanism=''
LAST_POINT=.False.
CONSTRAIN_NOX=.False.
CONSTRAIN_RUN=.FALSE.
SAVE_LEGACY=.false.
time=tstart
!tuv defined globally


!spinup time 0 for false, number timesteps otherwise


!dt is the output timestep and the timestep between times
!rate constants and notably photolysis rates are calcualted " 600 = ten minutes
dt = 600.
nc_set=1!36 ! the grouping factor that decides how often to write to file (modulo daycounter/nc_set)
!use nc_set = 0 for a single memory dump at the end of simulation.


spinup=0 !allocate this within ics
allocate(spin_const(NVAR))

call getarg(1,counter)!name
call getarg(2,ln)!location in Init Cons
read(ln, *) line

call getarg(3,ln)
read(ln, *) obs



!set OUTPUT_UNIT to 6 in globals file.
open(UNIT=output_unit,FILE='Outputs/'//trim(counter)//'.sdout')


!all initialisation calculations:
INCLUDE './src/initialisations.inc'


!i'cs copied from python initiation program

!so T=0 of the output file gives the initial condition
!i'cs copied from python initiation program

WRITE (SPEC_UNIT) newtime,LAT, LON, PRESS, TEMP,H2O, CFACTOR, RO2, C(:NSPEC)
WRITE (RATE_UNIT) newtime,LAT, LON, PRESS, TEMP, M, RCONST(:NREACT)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
TEND = TEND + spinup*dt !additional spinup time added if included


print*,  achar(27)//'[2J'
print*, 'Run Progress'
time_loop: DO WHILE (time < TEND)! This is the main loop for integrations
run_counter = run_counter+1

!!!~~~~ BEGIN DEPOSITION CONSTRAINTS Edwards et al. (2014) ~~~~!!!
 DEPOSITION=(1e-6+1.8e-7)*1.5
IF (MOD(TIME/(60*60),24.) .GT. 10. .AND. &
     MOD(TIME/(60*60),24.) .LT. 19.) THEN
   DEPOSITION=(1e-6+1.8e-7)*10.*1.5
ENDIF

O3DEPOSITION=(1.e-6+1.8e-7-3.e-7)*1.5
IF (MOD(TIME/(60*60),24.) .GT. 10. .AND. &
     MOD(TIME/(60*60),24.) .LT. 19.) THEN
   O3DEPOSITION=(1.e-6+1.8e-7-3.e-7)*10.*1.5
ENDIF

HNO3DEPOSITION=(1e-6+1.8e-7)*1.5
IF (MOD(TIME/(60*60),24.) .GT. 10. .AND. &
     MOD(TIME/(60*60),24.) .LT. 19.5) THEN
           HNO3DEPOSITION=(1e-6+1.8e-7)*10.*1.5
ENDIF

PNADEPOS = DEPOSITION
N2O5DEPOSITION=0.

N2O5HYD=(1.2e-3)
IF (MOD(TIME/(60*60),24.) .GT. 1. .AND. &
     MOD(TIME/(60*60),24.) .LT. 15.5) THEN
   N2O5HYD=(9.e-3)
ENDIF

KNOSCAV=1.0e-15
KNO2SCAV=1.0e-16
NO2SCAVDEP=0.
NOSCAVDEP=0.
!!!~~~~ END DEPOSITION CONSTRAINTS Edwards et al. (2014) ~~~~!!!

CALL Update_RCONST()! Update the rate constants

ICNTRL_U(:)=0
CALL INTEGRATE( TIN = time, TOUT = time+DT, RSTATUS_U = RSTATE, &! Integrate the model +1 timestep
ICNTRL_U = ICNTRL_U,IERR_U=ERROR)

! Traps for NaN
DO I=1,NVAR
IF (ISNAN(C(I)) .or. (ERROR .NE. 1)) then
    write (OUTPUT_UNIT,*) 'Integration error / NaN, Skipping Point'
    !print *, ''//achar(27)//'[95m', c,''//achar(27)//'[97m'
    !if (i .eqv. ind_DUMMY) then
    !    c(i)= 0.
    !    cycle
    !end if

    C(1:NVAR)=0.
    GOTO 1000
ENDIF
End Do

! Update the time to reflect the integration has taken place
time = RSTATE(1)
Daycounter=Daycounter+1


IF (CONSTRAIN_NOX) THEN
    write (OUTPUT_UNIT,*) 'Constraining NOx'
    TNOX=0.! Calcualte the total NOx in the box
    TNOX=TNOX+sum(C(1:NVAR)*NOX(1:NVAR))
    ! Update all NOx variables so that the total NOx in the box is the same as it was
    DO I=1,NVAR
      IF (NOX(I) .NE. 0)  C(I)=C(I)*TNOX_OLD/TNOX
    ENDDO
ENDIF


!!!~~~~ BEGIN SAVE CONSTRAINTS TO DSMACC/KPP VARIABLES Edwards et al. (2014) ~~~~!!!
! mje for Pete update
! tstart is local time, but Pete's numbers are GMT
   HALFHOUR=INT(MOD((TIME-LON/360.*24*60*60)*2./(60.*60.),48.))+1
   TEMP=OBS_TEMP(HALFHOUR)+273.15
   C(ind_HONO)=HONO(HALFHOUR)*1e-12*CFACTOR
   C(ind_CLNO2)=CLNO2(HALFHOUR)*1e-12*CFACTOR*0.6
   C(ind_DIME35EB)=DIME35EB(HALFHOUR)*1e-12*CFACTOR
   C(ind_DIET35TOL)=DIET35TOL(HALFHOUR)*1e-12*CFACTOR
!!!~~~~ END SAVE CONSTRAINTS TO DSMACC/KPP VARIABLES Edwards et al. (2014) ~~~~!!!



!!If constrain species concentrations if necessary

DO I=1,NVAR
    IF (CONSTRAIN(I) .GT. 0) C(I)=CONSTRAIN(I)
END DO






newtime = Jday*86400 + DAYCOUNTER*dt

WRITE (SPEC_UNIT) newtime,LAT, LON, PRESS, TEMP,H2O, CFACTOR, RO2, C(:NSPEC)
WRITE (RATE_UNIT) newtime,LAT, LON, PRESS, TEMP, M, RCONST(:NREACT)


if (run_counter > nc_set) then !increased at start

    DO i = 1,daycounter
    !write (SPEC_UNIT, 999) output_s(i,:) !(output_s(i,ij), ij=1,NSPEC)
    !write (RATE_UNIT, 999) output_r(i,:) !(output_r(i,ij), ij=1,NREACT)
    end do

    if (mod(run_counter/nc_set,20)==0) then
        print*, achar(27)//'[f', achar(27)//'['//trim(ln)//'B', achar(27)//'[94m |',repeat('#',floor(time/TEND*20)), repeat(' ',int(20-floor(time/TEND*20))),'| '//achar(27)//'[97m'//trim(counter)
    end if

else
    run_counter = 1
    output_s(:,:)=0.
    output_r(:,:)=0.
    !SaveOut(run_counter)


end if


! If we are doing a constrained run we need to store the diurnal profile of all the species
IF (CONSTRAIN_RUN .EQv. .TRUE.) THEN
    DIURNAL_NEW(1:NVAR,DAYCOUNTER)=C(1:NVAR)
    DIURNAL_RATES(1:NREACT,DAYCOUNTER)=RCONST(1:NREACT)
    FULL_CONCS(1:NVAR,full_counter+daycounter)=C(1:NVAR)
    full_counter=full_counter+daycounter+1.
    ! Are we at the end of a day?
    ! If so we need to 1) fiddle with the NOX to ensure it has the right concentrations see if we have reached a steady state
    IF (DAYCOUNTER*DT .GE. 24.*60.*60.) THEN
        ! Sort out the NOx. Need to increase the NOx concentration so that the constrained species is right
        ! What is  the constrained NOx species? Put result into CONSTNOXSPEC
        ! Calculate the ratio between the value we the constrained NOx species and what we have
        ! Remember the constrained NOx species is given by the negative constrained value
        DO I=1,NVAR
        IF (CONSTRAIN(I) .LT. 0)     CONSTNOXSPEC=I
        IF (NOX(I) .NE. 0)        C(I)=C(I)*NOXRATIO
        ENDDO
        NOXRATIO=-CONSTRAIN(CONSTNOXSPEC)/C(CONSTNOXSPEC)
        ! Multiply all the NOx species by the ratio so
        ! Update the total amount of NOx in box
        TNOX_OLD=TNOX_OLD*NOXRATIO
        ! Lets see how much the diurnal ratios have changed since the last itteration
        ! Frac diff is our metric for how much it has changed
        FRACDIFF=0.
        FRACCOUNT=0.
        ! Add up for all species and for each time point in the day
        DO I=1,NVAR
        DO JK=1,DAYCOUNTER
        !If there is a concentration calculated
        IF (DIURNAL_NEW(I,JK) .GT. 1.e2 .AND. &
        TRIM(SPC_NAMES(I)) .NE. 'DUMMY') THEN
            !Calculate the absolute value of the fractional difference and add it on
            ! Increment the counter to calculate the average
            FRACDIFF=FRACDIFF+&
            ABS(DIURNAL_OLD(I,JK)-DIURNAL_NEW(I,JK))/&
            DIURNAL_NEW(I,JK)
            FRACCOUNT=FRACCOUNT+1
        ENDIF
        ENDDO
        ENDDO

    FRACDIFF=FRACDIFF/FRACCOUNT !average fractional difference

    write (OUTPUT_UNIT,*) 'Fraction difference in the diurnal profile:', FRACDIFF! Output diagnostic


    ! Store the new diurnal profile as the old one so we can compare with the next day
    DIURNAL_OLD(1:NVAR,1:Daycounter)=DIURNAL_NEW(1:NVAR,1:DAYCOUNTER)

    print *, 'fractional difference aim 0 :', (fracdiff - 1e-3)
    IF (FRACDIFF .LE. 1e-3)  GOTO 1000    ! if system has converged end simulation


    DAYCOUNTER=0! reset the day counter to 0
    OLDFRACDIFF=FRACDIFF
    ENDIF

ENDIF


ENDDO time_loop
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

1000  print *, ' '//achar(27)//'[91m exit condition 1000 '//achar(27)//'[97m'
!if (CONSTRAIN_RUN .EQ. .true.) .AND. (OUTPUT_LAST .EQ. .false.)  print*, 'why not SaveOut(run_counter) work dammit'




    !NETCDF // WRITE TO FILE //
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    print*, 'special case if so, deallocate and allocate array'



    if(OUTPUT_LAST .eqv. .true.) then
        do I=1,DAYCOUNTER
        !not yet tested.
        output_s(I,1)= Jday*86400  + I*dt
        output_s(I,2)=Lat
        output_s(I,3)=Lon
        output_s(I,4)=Press
        output_s(I,5)=Temp
        output_s(I,6)=H2O
        output_s(I,7)=Cfactor

        output_s(I,8)=RO2
        output_s(I,9:)=(DIURNAL_NEW(1:NVAR,I))

        output_r(I,1)= Jday*86400  + I*dt
        output_r(I,2)=Lat
        output_r(I,3)=Lon
        output_r(I,4)=Press
        output_r(I,5)=Temp
        output_r(I,6)=H2O
        output_r(I,7)=Cfactor
        output_r(I,8:)=(DIURNAL_RATES(1:NREACT,I))
        enddo
    endif



    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !adjust for special cases


    CLOSE(10)
    CLOSE(12)



    close(output_unit)
    deallocate(output_s)
    deallocate(output_r)
    deallocate(s_names)
    deallocate(r_names)


    END PROGRAM driver
