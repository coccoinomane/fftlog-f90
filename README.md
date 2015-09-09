FFTLog-f90 is the Fortran 90 version of [FFTLog], a Fortran 77 code by [Andrew Hamilton] to compute the Fourier-Bessel, or Hankel, transform of a periodic sequence of logarithmically spaced points.

## QUICK START
To use FFTLog-f90:

1. Customise the Makefile to use your preferred Fortran compiler.

2. Run `make fftlog-f90`.

3. If you use parallelisation with iFort, make sure to place the file `libiomp5.dylib` either in your $PATH or in the working directory; you can find said file in your ifort folder, usually `/opt/intel/` in a Unix system.

4. To run a quick test of FFTLog-f90, execute the command `./fftlog-f90 p_kappa.dat xi.dat`, which computes the Fourier transform of the function contained in `p_kappa.dat`. This is the power spectrum P(k) of the matter field in our Universe, estimated using the [CLASS] code for the standard LCDM cosmological model. The Fourier transform of P(k) is the correlation function xi(r), which will be stored in the file `xi.dat` and will have the domain 1/k_max < r < 1/k_min. The results at the edges should not be trusted; see the documentation of [FTTLog] for further details.

5. You can then plot the result in [Gnuplot] with `set log; plot [0.01:1e4] "xi.dat" u 1:(abs($2)) w li`. The feature at r~120 Mpc is called the baryon acoustic peak; you can zoom in it with `unset log; plot [80:180] "xi.dat" u 1:2 w li`.

The function supports spline integration of the result with the syntax

    ./fftlog-f90 in.dat out.dat N_POINTS

where `N_POINTS` is the number of points you want in the output file. For example, with

    ./fftlog-f90 p_kappa.dat xi_splined.dat 2000

the peak is much smoother:

    unset log; plot [80:180] "xi_splined.dat" u 1:2 w li

FFTLog-f90 also supports custom integration ranges and arbitrary Bessel order. Refer to the usage message in the file `fftlog_driver.f90` for more details:

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


## CONTRIBUTE
Please feel free to improve FFTLog-f90. You can do so via the Github project page: <https://github.com/coccoinomane/fftlog-f90>. Fork the repository, make your modifications and send a pull request.

          
[FFTLog]: http://casa.colorado.edu/~ajsh/FFTLog
[Andrew Hamilton]: http://casa.colorado.edu/~ajsh
[Gnuplot]: http://www.gnuplot.info/