!=======================================================================
! * Driver for FFTLog
! *   Copyright by Guido Walter Pettinari
! *   Last modofied on 30 June 2013
!=======================================================================


PROGRAM fftlog_f90

	IMPLICIT NONE

!	If the input file has N_MAX lines then the calculation would take 2GB of RAM...
!	it's better to avoid it!
	INTEGER, PARAMETER :: N_MAX=67108864
	REAL ( KIND = 8 ), PARAMETER :: PI=3.141592653589793238462643383279502884197d0, CONSTANT=1/(2*PI*PI)
!	xx and yy will contain input; after FFTL has been called, they will contain output (if INTERPOLATION = 0)
	REAL ( KIND = 8 ), DIMENSION(:), ALLOCATABLE :: xx, yy, yy_second_derivative, yy_interp, wsave
	REAL ( KIND = 8 ) :: xx_inf_limit = -1d300, xx_sup_limit = 1d300, dln_xx, dlog_xx, log_xmedian,&
                       &log_xxmedian, log_xxmax, log_xxmin, rk, kr, mu, CENTRAL_INDEX, q, x, temp1, temp2
	INTEGER :: N_ARGUMENTS, DIRECTION, KR_OPTION, UNIT, error, N_INFILE=0, N_USED=0, FIRST_USED_INDEX = 0,&
             &LAST_USED_INDEX=0, N_OUTFILE=0, i=0
	LOGICAL :: INTERPOLATION = .FALSE., FFT_OK
	CHARACTER(512) :: input_filename, output_filename
	CHARACTER(128) :: buffer
	
	N_ARGUMENTS = COMMAND_ARGUMENT_COUNT()
	
	IF( N_ARGUMENTS.GT.7 .OR. N_ARGUMENTS.LT.2 ) STOP &
	                & "I will calculate the Fourier Bessel integral of the input file&
					& for various values of frequency. Please give input_filename and output_filename as arguments.&
          & If a positive number is provided as the third argument (N_OUTFILE), the output file will contain&
          & that many lines; these extra points are generated using spline interpolation.&
					& 4th optional argument is the bessel function order (default = 0.5, that is sine);&
					& 5th optional argument is the bias (default = 0.0);&
					& 6th & 7th optional arguments are inf and sup integration limits (default is the whole file domain).&
					& If you do not want interpolation, just insert a value <= 0 as third argument.&
          & The reference for the FFTlog algorithm by Andrew Hamilton can be found here&
          & http://casa.colorado.edu/~ajsh/FFTLog"


	CALL GET_COMMAND_ARGUMENT( 1, input_filename )
	CALL GET_COMMAND_ARGUMENT( 2, output_filename )

	IF ( N_ARGUMENTS .GE. 3 ) THEN
		CALL GET_COMMAND_ARGUMENT( 3, buffer )
		READ ( buffer, * ) N_OUTFILE
	END IF

	!        order of Bessel function ( mu = 0.5 for sine transform )
  mu = 0.5d0
	IF ( N_ARGUMENTS .GE. 4 ) THEN
		CALL GET_COMMAND_ARGUMENT( 4, buffer )
		READ ( buffer, * ) mu
	END IF

	!        bias exponent: q = 0 is unbiased ( good for power spectrum <-> correlation function )
  q = 0.d0
	IF ( N_ARGUMENTS .GE. 5 ) THEN
		CALL GET_COMMAND_ARGUMENT( 5, buffer )
		READ ( buffer, * ) q
	END IF

	IF ( N_ARGUMENTS == 6 ) THEN
		CALL GET_COMMAND_ARGUMENT( 6, buffer )
		READ ( buffer, * ) xx_sup_limit
	END IF
		
	IF ( N_ARGUMENTS == 7 ) THEN
		CALL GET_COMMAND_ARGUMENT( 6, buffer )
		READ ( buffer, * ) xx_inf_limit
		CALL GET_COMMAND_ARGUMENT( 7, buffer )
		READ ( buffer, * ) xx_sup_limit
	END IF
	
!	Select the entries of input_filename to process
	OPEN ( UNIT=50, FILE=input_filename, STATUS='OLD', IOSTAT=error )
	IF ( error /= 0 ) STOP "Input file could not be opened, exiting..."
	DO
		READ (UNIT=50, FMT=*, IOSTAT=error) temp1
		IF ( error < 0 ) EXIT
		N_INFILE = N_INFILE + 1		
		IF ( temp1 < xx_inf_limit .OR. temp1 > xx_sup_limit ) CYCLE
		N_USED = N_USED + 1
		LAST_USED_INDEX = N_INFILE
	END DO
	CLOSE(UNIT = 50)
	FIRST_USED_INDEX = LAST_USED_INDEX - N_USED + 1
	
!	Control block	
	IF ( N_USED == 0 ) STOP "Either the input file is empty or incompatible with the specified integration limits."	
	IF ( N_OUTFILE .LE. 0 ) THEN	
		N_OUTFILE = N_USED
		INTERPOLATION = .FALSE.
		PRINT *, "No interpolation of input data will be performed."
	ELSE
		INTERPOLATION = .TRUE.
	END IF		
    IF ( (N_USED > N_MAX) .OR. (N_OUTFILE > N_MAX) ) &
    &STOP "The number you specified or the number of elements of the input file are greater that N_MAX"

!	It is now safe to allocate memory to the arrays
	ALLOCATE( xx(N_USED), yy(N_USED), yy_second_derivative(N_USED) )
	ALLOCATE( yy_interp(N_OUTFILE), wsave ( 2*N_OUTFILE + 3*(N_OUTFILE/2) + 19 ) )

!	We now read the first two columsn of the input file, remembering our previously set integration limits
	OPEN ( UNIT=50, FILE=input_filename, STATUS='OLD', IOSTAT=error )
	IF ( error /= 0 ) STOP "Input file could not be opened, exiting..."
!	Let's dump the rows of the file that are smaller than xx_inf_limit
	DO i = 1, FIRST_USED_INDEX - 1
		READ ( 50, FMT=*, IOSTAT=error )
	END DO
	DO i = 1, N_USED
		READ ( 50, FMT=*, IOSTAT=error ) temp1, temp2
		IF ( error < 0 ) EXIT
		xx(i) = temp1
		yy(i) = temp2
	END DO
	CLOSE( UNIT = 50 )
	
	PRINT *, "Inferior integration limit = ", xx(1)
	PRINT *, "corresponding in input file to row = ", FIRST_USED_INDEX
	PRINT *, "Superior integration limit = ", xx(N_USED)
	PRINT *, "corresponding in input file to row = ", LAST_USED_INDEX	
	PRINT *, "Order of Bessel function mu = ", mu
	PRINT *, "Number of elements in input file = ", N_INFILE
  PRINT *, "Number of elements used = ", N_USED
	PRINT *, "Number of elements in output file = ", N_OUTFILE

!---CREATION OF INPUT LOGARITHMIC ARRAY; SHOULD HAVE N_OUTFILE ELEMENTS INDEPENDENTLY OF INTERPOLATION OR NOT
	log_xxmin = LOG10( xx(1) )
	log_xxmax = LOG10( xx(N_USED) )

	dlog_xx = (log_xxmax-log_xxmin) / (N_OUTFILE-1)
	!  central index (1/2 integral if N_OUTFILE is even)
	CENTRAL_INDEX = dble( N_OUTFILE+1 ) / 2.d0
	
	!		 logarithmical spacing between points (needed by fhti). Was dlnr
	dln_xx = dlog_xx * LOG(10.d0)
	!        central point of periodic interval at log10 (xxmedian). Was logrc
	log_xxmedian = ( log_xxmin + log_xxmax) / 2.d0
	!        sensible approximate choice of k_c r_c
	kr = 1.d0
	!        tell fhti to change kr to low-ringing value
	KR_OPTION = 2
	!        forward transform
	DIRECTION = 1
	!		 central point in k-space. Was logkc
	log_xmedian = LOG10 (kr) - log_xxmedian
	!        rk = r_c/k_c
  rk = 10.d0 ** ( log_xxmedian - log_xmedian )


	IF( INTERPOLATION ) THEN
		! spline_cubic_set ( n, t, y, ibcbeg, ybcbeg, ibcend, ybcend, ypp )
		CALL SPLINE_CUBIC_SET ( N_USED, xx, yy, 0, 0.0d0, 0, 0.0d0, yy_second_derivative )
		DO i = 1, N_OUTFILE ! the x's are the input domain
			x = 10.d0 ** ( log_xxmedian + ( i - CENTRAL_INDEX )*dlog_xx )
			! spline_cubic_val ( n, t, y, ypp, tval, yval, ypval, yppval )
			CALL SPLINE_CUBIC_VAL ( N_USED, xx, yy, yy_second_derivative, x, yy_interp(i), temp1, temp2 )
			yy_interp(i) = yy_interp(i)*x
		END DO
		DEALLOCATE ( yy_second_derivative )
	ELSE
		DO i = 1, N_USED
			yy(i) = xx(i) * yy(i)
		END DO
	END IF
!--------initialize FFTLog transform - note fhti resets kr
	CALL FHTI ( N_OUTFILE, mu, q, dln_xx, kr, KR_OPTION, wsave, FFT_OK )

	IF( .NOT. FFT_OK ) STOP "FHTI not ok!"
!--------transform
	IF( INTERPOLATION ) THEN 
		CALL FFTL ( N_OUTFILE, yy_interp, rk, DIRECTION, wsave )
	ELSE
		CALL FFTL ( N_OUTFILE, yy, rk, DIRECTION, wsave )
	END IF
	
!--------print/write result
	OPEN( UNIT=40, FILE = output_filename, STATUS = "REPLACE", IOSTAT = error )
	IF ( error /= 0 ) STOP "Output file could not be opened, exiting..."
	IF ( INTERPOLATION ) THEN 
		DO i = 1, N_OUTFILE ! now the x's are the output domain
			x = 10.d0 ** ( log_xmedian + ( i - CENTRAL_INDEX ) * dlog_xx )
			WRITE ( 40, '(3g24.16)' ) x, CONSTANT * yy_interp(i)/x
		END DO
		DEALLOCATE ( yy_interp )
	ELSE
		DO i = 1, N_OUTFILE ! now the x's are the output domain
			x = 10.d0 ** ( log_xmedian + ( i - CENTRAL_INDEX )*dlog_xx )
			WRITE( 40, '(3g24.16)' ) x, CONSTANT * yy(i)/x
		END DO
	END IF
	DEALLOCATE(xx,yy,wsave)
	
END PROGRAM fftlog_f90