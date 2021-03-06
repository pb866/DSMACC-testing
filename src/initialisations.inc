!!!!!Input arguments!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!open the init cons file.
open(UNIT=21,FILE='Init_cons.dat')
READ(21,'(i10)') IntTime
! check if a special run has been selected.
if (INTTIME .LE. 0)  then
    CONSTRAIN_RUN=.TRUE.
    print*, 'reverting to legacy output.'
    SAVE_LEGACY= .TRUE.

    WRITE (OUTPUT_UNIT,*) 'Integration to convergence'

    if (INTTIME .EQ. -1) then
        INTTIME=int(4380*24.*60.*60.)
        OUTPUT_LAST=.FALSE.


    else if (INTTIME .EQ. -2) then
        INTTIME=int(50*24.*60.*60.)
        OUTPUT_LAST=.TRUE.


    else
        WRITE (OUTPUT_UNIT,*) 'Negative Integration Time \n &
        &But not a special case'
        STOP
    endif


endif

runtimestep=int(Inttime / dt) !used for netcdf


!read species and concentrations
!number of read == 10000


READ(21,'(*(a15,x))') spec_name
READ(21,'(*(i15,x))') const_method

print*, 'if -checks all this is an incorrect conversion'

if (line.gt.1) then !skip lines
    do 133 i = 1,(line-1)
    133 read(21,*)
end if

!get required concentrations
READ (21,'(*(e15.4,x))') concs
 close(21)


if (save_legacy) nc_set = 144 ! 1day
if (nc_set .eq. 0)  nc_set = runtimestep

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!Variable Initialisation !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
allocate(output_s(1:nc_set, Nspec+10))
allocate(output_r(1:nc_set, Nreact+6))
print*, 'rename ncset'
allocate(s_names(Nspec+10))
allocate(r_names(Nreact+6))
output_s(:,:) = 0.
output_r(:,:) = 0.
call initVal(concs)!initialise values
M   = CFACTOR
O2 = 0.21 * CFACTOR
N2 = 0.78 * CFACTOR
tstart = (mod(jday,1.))*24.*60.*60. ! time start
tstart = tstart+LON/360.*24.*60.*60.!conv 2 localtime
tend = tstart + IntTime  !end time
nc_counter=0


write (OUTPUT_UNIT,*) 'Initializing model ', counter
write (OUTPUT_UNIT,*) 'Starting Jday:',jday
write (OUTPUT_UNIT,*) 'Starting time:',tstart
write (OUTPUT_UNIT,*) 'Ending time:', tend
write (OUTPUT_UNIT,*) 'Time step:', dt
write (OUTPUT_UNIT,*) 'hvhvhvhvhvhvhvhvhvhvhvhvhvhvhvhvhvhvhvhvhvhvhvhvhvhvhv'
write (OUTPUT_UNIT,*) 'Using TUV to calculate photolysis rates as a function of SZA'
alta=(1-(press/1013.25)**0.190263)*288.15/0.00198122*0.304800/1000.
write (OUTPUT_UNIT,*) 'Aerosol surface area', SAREA
write (OUTPUT_UNIT,*) 'Aerosol particle radius 1', RP1
write (OUTPUT_UNIT,*) 'Altitude =', alta
write (OUTPUT_UNIT,*) 'Pressure =', Press
write (OUTPUT_UNIT,*) 'Temperature =', Temp
write (OUTPUT_UNIT,*) 'Latitude =', Lat
write (OUTPUT_UNIT,*) 'Lon =', Lon
write (OUTPUT_UNIT,*) 'Local Time =', Tstart/(60.*60.)
write (OUTPUT_UNIT,*) 'SZA =',ZENANG(int(jday),Tstart/(60.*60.),lat)*180./(4*ATAN(1.))

if (o3col .eq. 0) then
    o3col=260.
    write (OUTPUT_UNIT,*) 'Ozone column not specified using 260 Dobsons'
else
    write (OUTPUT_UNIT,*) 'Ozone column =', o3col
endif

if (albedo .eq. 0) then
    albedo=0.1
    write (OUTPUT_UNIT,*) 'Albedo not specified using 0.1'
else
    write (OUTPUT_UNIT,*) 'Albedo =', albedo
endif

call set_up_photol(O3col,Albedo,alta,temp,bs,cs,ds,szas,svj_tj)


write (OUTPUT_UNIT,*) 'Photolysis rates calculated'
write (OUTPUT_UNIT,*) 'hvhvhvhvhvhvhvhvhvhvhvhvhvhvhvhvhvhvhvhvhvhvhvhvhvhvhv'
time = tstart
OLDFRACDIFF=0.

! If NOx is being constrained calculate the total NOx in the model
IF (CONSTRAIN_NOX) THEN
    TNOX_OLD=0.
    TNOX_OLD=TNOX_OLD+sum(C(1:NVAR)*NOX(1:NVAR))
ENDIF

! Define the initial state of the model
STARTSTATE(1:NVAR)=C(1:NVAR)
DIURNAL_OLD(:,:)=0.
! Calculate clear sky photolysis rates
JFACTNO2=1.
JFACTO1D=1.

! Update the rate constants
CALL Update_RCONST()
! Calcualte correction factors for the model photolysis rates
if (JO1D .NE. 0.) then
    if (J(1) .GT. 0.) JFACTO1D=JO1D/J(1)
    if (JNO2 .NE. 0.) JFACTO1D=JFACTNO2
endif

if (JNO2 .NE. 0.) then
    if (J(4) .GT. 0.) JFACTNO2=JNO2/J(4)
    if (JO1D .NE. 0.) JFACTNO2=JFACTO1D
endif

! Fix JFACTO1D and JFACTNO2 to observed
 JFACTO1D=0.5167336
 JFACTNO2=0.6643193

write (OUTPUT_UNIT,*) 'JO1D Calc=', J(1)
write (OUTPUT_UNIT,*) 'JO1D Measre =', JO1D
write (OUTPUT_UNIT,*) 'Correction JO1D and JNO2 by', JFACTO1D,JFACTNO2

! Set up a counter to count the number time that the model has been run for
Daycounter=0
full_counter=0


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!NETCDF INITIALISATIONS

        filename='./Outputs/'//trim(counter)//'_'//trim(mechanism)//'.spec'
        open(SPEC_UNIT, file=filename,action="write",form='unformatted')
        filename='./Outputs/'//trim(counter)//'_'//trim(mechanism)//'.rate'
        open(RATE_UNIT, file=filename,action="write",form='unformatted')


        open (100,file='./Outputs/spec.names')
        write(100,*)'TIME,LAT,LON,PRESS,TEMP,H2O,M,RO2', (','//trim(SPC_NAMES(I)),I=1,NSPEC)
        close(100)
        open (100,file = './Outputs/rate.names')
        write(100,*)'TIME,LAT,LON,PRESS,TEMP,M', (','//trim(EQN_NAMES(I)),I=1,NREACT)


         999 FORMAT(*(E24.16,"!"))
