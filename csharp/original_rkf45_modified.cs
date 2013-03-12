// Adaptive RKF45 solver in C#, based on C code from
// http://people.sc.fsu.edu/~jburkardt/c_src/rkf45/rkf45.html

// sestoft@itu.dk * 2013-01-31

// Part of the Actulus project.

// DO NOT DISTRIBUTE: Contains project-internal information.

/* Original license and author statement for the RKF45 solver:

  Licensing:
    This code is distributed under the GNU LGPL license. 

  Modified:
    05 April 2011

  Author:
    Original FORTRAN77 version by Herman Watts, Lawrence Shampine.
    C++ version by John Burkardt.
*/

/* 

The coefficients used in this solver are from Fehlberg's 1969
publication (NASA TR R-315) although many better coefficient sets have
been found later.  See file rkfsolvers.txt.

*/


using System;
using System.Diagnostics;

class Rkf45 {

  static readonly double DoubleEpsilon = FindDoubleEpsilon();
  const int MAXNFE = 3000;

  private readonly Action<double, double[], double[]> f; 
  private readonly int neqn;
  private readonly double[] f1, f2, f3, f4, f5;

  // These were C "static" variables in R8_RKF45:
  private double abserr_save = -1.0;
  private int flag_save = -1000;
  private double h = -1.0;
  private bool init = false;
  private int kflag = -1000;  // The flag value (3, 4, 6) last returned
  private int kop = -1;
  private int nfe = -1;
  private double relerr_save = -1.0;
  private double remin = 1.0E-12;
  // This was a parameter in R8_RKF45:
  private readonly double[] yp;

  public Rkf45(Action<double, double[], double[]> f, int neqn) {
    this.f = f;
    if ( neqn < 1 )
      throw new Exception("There must be at least one equation >= 1");
    this.neqn = neqn;
    this.f1 = new double[neqn]; 
    this.f2 = new double[neqn]; 
    this.f3 = new double[neqn]; 
    this.f4 = new double[neqn]; 
    this.f5 = new double[neqn];
    this.yp = new double[neqn];
  }

  /******************************************************************************/

  private void r8_fehl (double[] y, double t, double h, double[] yp, double[] s)

  /******************************************************************************/
  /*
  Purpose:

    R8_FEHL takes one Fehlberg fourth-fifth order step.

  Discussion:

    This version of the routine uses DOUBLE real arithemtic.

    This routine integrates a system of NEQN first order ordinary differential
    equations of the form
      dY(i)/dT = F(T,Y(1:NEQN))
    where the initial values Y and the initial derivatives
    YP are specified at the starting point T.

    The routine advances the solution over the fixed step H and returns
    the fifth order (sixth order accurate locally) solution
    approximation at T+H in array S.

    The formulas have been grouped to control loss of significance.
    The routine should be called with an H not smaller than 13 units of
    roundoff in T so that the various independent arguments can be
    distinguished.


  Reference:

    Erwin Fehlberg,
    Low-order Classical Runge-Kutta Formulas with Stepsize Control,
    NASA Technical Report R-315, 1969.

    Lawrence Shampine, Herman Watts, S Davenport,
    Solving Non-stiff Ordinary Differential Equations - The State of the Art,
    SIAM Review,
    Volume 18, pages 376-411, 1976.

  Parameters:

    Input, external F, a user-supplied subroutine to evaluate the
    derivatives Y'(T), of the form:

      void f ( double t, double y[], double yp[] )

    Input, int NEQN, the number of equations to be integrated.

    Input, double Y[NEQN], the current value of the dependent variable.

    Input, double T, the current value of the independent variable.

    Input, double H, the step size to take.

    Input, double YP[NEQN], the current value of the derivative of the
    dependent variable.

    Output, double F1[NEQN], F2[NEQN], F3[NEQN], F4[NEQN], F5[NEQN], derivative
    values needed for the computation.

    Output, double S[NEQN], the estimate of the solution at T+H.
  */
  {
    double ch;
    ch = h / 4.0;
    for (int i = 0; i < neqn; i++ )
      f5[i] = y[i] + ch * yp[i];
    f ( t + ch, f5, f1 );

    ch = 3.0 * h / 32.0;
    for (int i = 0; i < neqn; i++ )
      f5[i] = y[i] + ch * ( yp[i] + 3.0 * f1[i] );
    f ( t + 3.0 * h / 8.0, f5, f2 );

    ch = h / 2197.0;
    for (int i = 0; i < neqn; i++ )
      f5[i] = y[i] + ch * ( 1932.0 * yp[i] + ( 7296.0 * f2[i] - 7200.0 * f1[i] ) );
    f ( t + 12.0 * h / 13.0, f5, f3 );
    
    ch = h / 4104.0;
    for (int i = 0; i < neqn; i++ )
      f5[i] = y[i] + ch * ( ( 8341.0 * yp[i] - 845.0 * f3[i] ) + 
              ( 29440.0 * f2[i] - 32832.0 * f1[i] ) );
    f ( t + h, f5, f4 );
    
    ch = h / 20520.0;
    for (int i = 0; i < neqn; i++ )
      f1[i] = y[i] + ch * ( ( -6080.0 * yp[i] + 
	      ( 9295.0 * f3[i] - 5643.0 * f4[i] ) ) + ( 41040.0 * f1[i] - 28352.0 * f2[i] ) );
    f ( t + h / 2.0, f1, f5 );

    /*
      Ready to compute the approximate solution at T+H.
    */
    ch = h / 7618050.0;
    for (int i = 0; i < neqn; i++ )
      s[i] = y[i] + ch * ( ( 902880.0 * yp[i] + 
	     ( 3855735.0 * f3[i] - 1371249.0 * f4[i] ) ) + ( 3953664.0 * f2[i] + 277020.0 * f5[i] ) );
  }
  

  /******************************************************************************/

  // Main public function.

  public int r8_rkf45(double[] y, ref double t, double tout, 
		      ref double relerr, double abserr, int flag )

  /******************************************************************************/
  /*
    Purpose:

      R8_RKF45 carries out the Runge-Kutta-Fehlberg method.

    Discussion:

      This version of the routine uses DOUBLE real arithmetic.

      This routine is primarily designed to solve non-stiff and mildly stiff
      differential equations when derivative evaluations are inexpensive.
      It should generally not be used when the user is demanding
      high accuracy.

      This routine integrates a system of NEQN first-order ordinary differential
      equations of the form:

	dY(i)/dT = F(T,Y(1),Y(2),...,Y(NEQN))

      where the Y(1:NEQN) are given at T.

      Typically the subroutine is used to integrate from T to TOUT but it
      can be used as a one-step integrator to advance the solution a
      single step in the direction of TOUT.  On return, the parameters in
      the call list are set for continuing the integration.  The user has
      only to call again (and perhaps define a new value for TOUT).

      Before the first call, the user must 

      * supply the subroutine F(T,Y,YP) to evaluate the right hand side;
	and declare F in an EXTERNAL statement;

      * initialize the parameters:
	NEQN, Y(1:NEQN), T, TOUT, RELERR, ABSERR, FLAG.
	In particular, T should initially be the starting point for integration,
	Y should be the value of the initial conditions, and FLAG should 
	normally be +1.

      Normally, the user only sets the value of FLAG before the first call, and
      thereafter, the program manages the value.  On the first call, FLAG should
      normally be +1 (or -1 for single step mode.)  On normal return, FLAG will
      have been reset by the program to the value of 2 (or -2 in single 
      step mode), and the user can continue to call the routine with that 
      value of FLAG.

      (When the input magnitude of FLAG is 1, this indicates to the
      program that it is necessary to do some initialization work.  An
      input magnitude of 2 lets the program know that that
      initialization can be skipped, and that useful information was
      computed earlier.)  CHANGE: The setup implied by flag=1 should be
      done by solver object construction (initialization).

      The routine returns with all the information needed to continue
      the integration.  If the integration reached TOUT, the user need only
      define a new TOUT and call again.  In the one-step integrator
      mode, returning with FLAG = -2, the user must keep in mind that 
      each step taken is in the direction of the current TOUT.  Upon 
      reaching TOUT, indicated by the output value of FLAG switching to 2,
      the user must define a new TOUT and reset FLAG to -2 to continue 
      in the one-step integrator mode.

      In some cases, an error or difficulty occurs during a call.  In that case,
      the output value of FLAG is used to indicate that there is a problem
      that the user must address.  These values include:

      * 3, integration was not completed because the input value of RELERR, the 
	relative error tolerance, was too small.  RELERR has been increased 
	appropriately for continuing.  If the user accepts the output value of
	RELERR, then simply reset FLAG to 2 and continue.  CHANGE: Warning.

      * 4, integration was not completed because more than MAXNFE derivative 
	evaluations were needed.  This is approximately (MAXNFE/6) steps.
	The user may continue by simply calling again.  The function counter 
	will be reset to 0, and another MAXNFE function evaluations are allowed. 
	CHANGE: Warning.

      * 5, integration was not completed because the solution vanished, 
	making a pure relative error test impossible.  The user must use 
	a non-zero ABSERR to continue.  Using the one-step integration mode 
	for one step is a good way to proceed.  CHANGE: Warning.

      * 6, integration was not completed because the requested accuracy 
	could not be achieved, even using the smallest allowable stepsize. 
	The user must increase the error tolerances ABSERR or RELERR before
	continuing.  It is also necessary to reset FLAG to 2 (or -2 when 
	the one-step integration mode is being used).  The occurrence of 
	FLAG = 6 indicates a trouble spot.  The solution is changing 
	rapidly, or a singularity may be present.  It often is inadvisable 
	to continue.  CHANGE: Exception.

      * 7, it is likely that this routine is inefficient for solving
	this problem.  Too much output is restricting the natural stepsize
	choice.  The user should use the one-step integration mode with 
	the stepsize determined by the code.  If the user insists upon 
	continuing the integration, reset FLAG to 2 before calling 
	again.  Otherwise, execution will be terminated.  

      * 8, invalid input parameters, indicates one of the following:
	NEQN <= 0;
	T = TOUT and |FLAG| /= 1;
	RELERR < 0 or ABSERR < 0;
	FLAG == 0  or FLAG < -2 or 8 < FLAG.
	CHANGE: Exception.

    Licensing:

      This code is distributed under the GNU LGPL license. 

    Modified:

      27 March 2004

    Author:

      Original FORTRAN77 version by Herman Watts, Lawrence Shampine.
      C++ version by John Burkardt.

    Reference:

      Erwin Fehlberg,
      Low-order Classical Runge-Kutta Formulas with Stepsize Control,
      NASA Technical Report R-315, 1969.

      Lawrence Shampine, Herman Watts, S Davenport,
      Solving Non-stiff Ordinary Differential Equations - The State of the Art,
      SIAM Review,
      Volume 18, pages 376-411, 1976.

    Parameters:

      Input, external F, a user-supplied subroutine to evaluate the
      derivatives Y'(T), of the form:

	void f ( double t, double y[], double yp[] )

      Input, int NEQN, the number of equations to be integrated.

      Input/output, double Y[NEQN], the current solution vector at T.

      Input/output, double YP[NEQN], the derivative of the current
      solution vector at T.  The user should not set or alter this
      information!  CHANGE: This looks like a candidate for being a
      field, holding in-between-calls information, rather than an
      input/output parameter.  DONE.

      Input/output, ref double T, the current value of the independent
      variable.  CHANGE: This should be a ref parameter.

      Input, double TOUT, the output point at which solution is desired.  
      TOUT = T is allowed on the first call only, in which case the routine
      returns with FLAG = 2 if continuation is possible.

      Input, ref double RELERR, double ABSERR, the relative and absolute error 
      tolerances for the local error test.  At each step the code requires:
	abs ( local error ) <= RELERR * abs ( Y ) + ABSERR
      for each component of the local error and the solution vector Y.
      RELERR cannot be "too small".  If the routine believes RELERR has been
      set too small, it will reset RELERR to an acceptable value and return
      immediately for user action.
      CHANGE: The RELERR should be a ref parameter.

      Input, int FLAG, indicator for status of integration. On the first
      call, set FLAG to +1 for normal use, or to -1 for single step
      mode.  On subsequent continuation steps, FLAG should be +2, or -2
      for single step mode.  CHANGE: Keep this flag to begin with, but
      do not consider the initialization case (+1 and -1).

      Output, int RKF45_D, indicator for status of integration.  A value of 2 
      or -2 indicates normal progress, while any other value indicates a 
      problem that should be addressed.
  */
  {
    double ae;
    double dt;
    double ee;
    double eeoet;
    double eps;
    double esttol;
    double et;
    int flag_return;
    bool hfaild;
    double hmin;
    int mflag;
    bool output;
    double relerr_min;
    double s;
    double scale;
    double tol;
    double ypk;

    flag_return = flag;
    /*
      Check the input parameters.
    */
    eps = DoubleEpsilon;

    if ( relerr < 0.0 )
      throw new Exception("relative error is negative");

    if ( abserr < 0.0 )
      throw new Exception("absolute error is negative");

    if ( flag_return == 0 || flag_return == 1 || flag_return == -1
	 || 8 < flag_return  || flag_return < -2 )
      throw new Exception("illegal flag value");

    mflag = Math.Abs( flag_return );

  /*
    Is this a continuation call?
  */
    if ( mflag != 1 )
    {
      if ( t == tout && kflag != 3 )
      {
	flag_return = 8;
	throw new Exception("Warning: t == tout");
	// return flag_return;
      }
  /*
    FLAG = -2 or +2:
  */
      if ( mflag == 2 )
      {
	if ( kflag == 3 )
	{
	  flag_return = flag_save;
	  mflag = Math.Abs( flag_return );
	}
	else if ( !init )
	{
	  flag_return = flag_save;
	}
	else if ( kflag == 4 )
	{
	  nfe = 0;
	}
	else if ( kflag == 5 && abserr == 0.0 )
	  throw new Exception("R8_RKF45 - Fatal error!\nKFLAG = 5 and ABSERR = 0.0.\n" );
	else if ( kflag == 6 && relerr <= relerr_save && abserr <= abserr_save )
	  throw new Exception("R8_RKF45 - Fatal error!\nKFLAG = 6 and\nRELERR <= RELERR_SAVE and\nABSERR <= ABSERR_SAVE.\n" );
      }
  /*
    FLAG = 3, 4, 5, 6, 7 or 8.
  */
      else
      {
	if ( flag_return == 3 )
	{
	  flag_return = flag_save;
	  if ( kflag == 3 )
	  {
	    mflag = Math.Abs( flag_return );
	  }
	}
	else if ( flag_return == 4 )
	{
	  nfe = 0;
	  flag_return = flag_save;
	  if ( kflag == 3 )
	  {
	    mflag = Math.Abs( flag_return );
	  }
	}
	else if ( flag_return == 5 && 0.0 < abserr )
	{
	  flag_return = flag_save;
	  if ( kflag == 3 )
	  {
	    mflag = Math.Abs( flag_return );
	  }
	}
  /*
    Integration cannot be continued because the user did not respond to
    the instructions pertaining to FLAG = 5, 6, 7 or 8.
  */
	else
	  throw new Exception("R8_RKF45 - Fatal error!R4_RKF45 - Fatal error!\nIntegration cannot be continued.\nThe user did not respond to the output value\nFLAG = 5, 6, 7 or 8.\n" );
      }
    }
  /*
    Save the input value of FLAG.  
    Set the continuation flag KFLAG for subsequent input checking.
  */
    flag_save = flag_return;
    kflag = 0;
  /*
    Save RELERR and ABSERR for checking input on subsequent calls.
  */
    relerr_save = relerr;
    abserr_save = abserr;
  /*
    Restrict the relative error tolerance to be at least 

      2*EPS+REMIN 

    to avoid limiting precision difficulties arising from impossible 
    accuracy requests.
  */
    relerr_min = 2.0 * DoubleEpsilon + remin;
  /*
    Is the relative error tolerance too small?
  */
    if ( relerr < relerr_min )
    {
      relerr = relerr_min;
      kflag = 3;
      flag_return = 3;
      throw new Exception("Warning: RELERR too small");
      // return flag_return;
    }

    dt = tout - t;
  /*
    Initialization:

    Set the initialization completion indicator, INIT;
    set the indicator for too many output points, KOP;
    evaluate the initial derivatives
    set the counter for function evaluations, NFE;
    estimate the starting stepsize.
  */

    // CHANGE: Get rid of the explicit initialization step, by
    // assuming that init is false when the object is created, by
    // testing init in every call, and performing the requisite
    // initialization if needed.
    
    if ( !init )
    {
      init = true;
      kop = 0;
      f ( t, y, yp );
      nfe = 1;

      if ( t == tout )
      {
	flag_return = 2;
	return flag_return;
      }

      h = Math.Abs( dt );
      double toln = 0.0;

      for (int k = 0; k < neqn; k++ )
      {
	tol = relerr * Math.Abs( y[k] ) + abserr;
	if ( 0.0 < tol )
	{
	  toln = tol;
	  ypk = Math.Abs( yp[k] );
	  if ( tol < ypk * Math.Pow( h, 5 ) )
	  {
	    h = Math.Pow( ( tol / ypk ), 0.2 );
	  }
	}
      }

      if ( toln <= 0.0 )
      {
	h = 0.0;
      }

      h = Math.Max ( h, 26.0 * eps * Math.Max ( Math.Abs( t ), Math.Abs( dt ) ) );

      if ( flag_return < 0 )
      {
	flag_save = -2;
      }
      else
      {
	flag_save = 2;
      }
    } // end of initialization

  /*
    Set stepsize for integration in the direction from T to TOUT.
  */
    h = r8_sign ( dt ) * Math.Abs( h );
  /*
    Test to see if too may output points are being requested.
  */
    if ( 2.0 * Math.Abs( dt ) <= Math.Abs( h ) )
    {
      kop = kop + 1;
    }
  /*
    Unnecessary frequency of output.
  */
    if ( kop == 100 )
    {
      kop = 0;
      // freed fjs
      flag_return = 7;
      throw new Exception("Warning: Output resolution very high");
      // return flag_return;
    }
  /*
    If we are too close to the output point, then simply extrapolate and return.
  */
    if ( Math.Abs( dt ) <= 26.0 * eps * Math.Abs( t ) )
    {
      t = tout;
      for (int i = 0; i < neqn; i++ )
      {
	y[i] = y[i] + dt * yp[i];
      }
      f ( t, y, yp );
      nfe = nfe + 1;

      // freed fjs
      flag_return = 2;
      return flag_return;
    }
  /*
    Initialize the output point indicator.
  */
    output = false;
  /*
    To avoid premature underflow in the error tolerance function,
    scale the error tolerances.
  */
    scale = 2.0 / relerr;
    ae = scale * abserr;
  /*
    Step by step integration.
  */
    for ( ; ; )
    {
      hfaild = false;
  /*
    Set the smallest allowable stepsize.
  */
      hmin = 26.0 * eps * Math.Abs( t );
  /*
    Adjust the stepsize if necessary to hit the output point.

    Look ahead two steps to avoid drastic changes in the stepsize and
    thus lessen the impact of output points on the code.
  */
      dt = tout - t;

      if ( 2.0 * Math.Abs( h ) <= Math.Abs( dt ) )
      {
      }
      else
  /*
    Will the next successful step complete the integration to the output point?
  */
      {
	if ( Math.Abs( dt ) <= Math.Abs( h ) )
	{
	  output = true;
	  h = dt;
	}
	else
	{
	  h = 0.5 * dt;
	}

      }
  /*
    Here begins the core integrator for taking a single step.

    The tolerances have been scaled to avoid premature underflow in
    computing the error tolerance function ET.
    To avoid problems with zero crossings, relative error is measured
    using the average of the magnitudes of the solution at the
    beginning and end of a step.
    The error estimate formula has been grouped to control loss of
    significance.

    To distinguish the various arguments, H is not permitted
    to become smaller than 26 units of roundoff in T.
    Practical limits on the change in the stepsize are enforced to
    smooth the stepsize selection process and to avoid excessive
    chattering on problems having discontinuities.
    To prevent unnecessary failures, the code uses 9/10 the stepsize
    it estimates will succeed.

    After a step failure, the stepsize is not allowed to increase for
    the next attempted step.  This makes the code more efficient on
    problems having discontinuities and more effective in general
    since local extrapolation is being used and extra caution seems
    warranted.

    Test the number of derivative function evaluations.
    If okay, try to advance the integration from T to T+H.
  */
      for ( ; ; )
      {
	/*
	  Have we done too much work?
	*/
	if ( MAXNFE < nfe )
	{
	  kflag = 4;
	  // freed fjs
	  flag_return = 4;
	  throw new Exception("Warning: Too many function evaluations");
	  // return flag_return;
	}
	/*
	  Advance an approximate solution over one step of length H.
	*/
	r8_fehl (y, t, h, yp, f1);
	nfe = nfe + 5;
	/*
	  Compute and test allowable tolerances versus local error estimates
	  and remove scaling of tolerances.  The relative error is
	  measured with respect to the average of the magnitudes of the
	  solution at the beginning and end of the step.
	*/
	eeoet = 0.0;

	for (int k = 0; k < neqn; k++ )
	{
	  et = Math.Abs( y[k] ) + Math.Abs( f1[k] ) + ae;

	  if ( et <= 0.0 )
	  {
	    // freed fjs
	    flag_return = 5;
	    throw new Exception("Warning: Pure relative error test impossible");
	    // return flag_return;
	  }

	  ee = Math.Abs 
	  ( ( -2090.0 * yp[k] 
	    + ( 21970.0 * f3[k] - 15048.0 * f4[k] ) 
	    ) 
	  + ( 22528.0 * f2[k] - 27360.0 * f5[k] ) 
	  );

	  eeoet = Math.Max ( eeoet, ee / et );

	}

	esttol = Math.Abs( h ) * eeoet * scale / 752400.0;

	if ( esttol <= 1.0 )
	{
	  break;
	}
  /*
    Unsuccessful step.  Reduce the stepsize, try again.
    The decrease is limited to a factor of 1/10.
  */
	hfaild = true;
	output = false;

	if ( esttol < 59049.0 )
	{
	  s = 0.9 / Math.Pow( esttol, 0.2 );
	}
	else
	{
	  s = 0.1;
	}

	h = s * h;

	if ( Math.Abs( h ) < hmin )
	{
	  kflag = 6;
	  // freed fjs
	  flag_return = 6;
	  throw new Exception("Warning: Requested accuracy could not be achieved");
	  // return flag_return;
	}

      }
  /*
    We exited the loop because we took a successful step.  
    Store the solution for T+H, and evaluate the derivative there.
  */
      t = t + h;
      for (int i = 0; i < neqn; i++ )
      {
	y[i] = f1[i];
      }
      f ( t, y, yp );
      nfe = nfe + 1;
  /*
    Choose the next stepsize.  The increase is limited to a factor of 5.
    If the step failed, the next stepsize is not allowed to increase.
  */
      if ( 0.0001889568 < esttol )
      {
	s = 0.9 / Math.Pow( esttol, 0.2 );
      }
      else
      {
	s = 5.0;
      }

      if ( hfaild )
      {
	s = Math.Min ( s, 1.0 );
      }

      h = r8_sign ( h ) * Math.Max ( s * Math.Abs( h ), hmin );
  /*
    End of core integrator

    Should we take another step?
  */
      if ( output )
      {
	t = tout;
	// freed fjs
	flag_return = 2;
	return flag_return;
      }

      if ( flag_return <= 0 )
      {
	// freed fjs
	flag_return = -2;
	return flag_return;
      }

    }
  }

  /******************************************************************************/

  /* Auxiliaries */

  static double r8_sign ( double x ) {
    return (double)(Math.Sign(x));
  }

  static double FindDoubleEpsilon() {
    double r = 1.0;
    while (1.0 < (1.0 + r))
      r = r / 2.0;
    return 2.0 * r;
  }

  /******************************************************************************/

  /* Main public solver functions for actuarial use */

  public static double[][] RKF45_n(Action<double, double[],double[]> dV, 
				   Action<double,double[]> bj_ii,
				   int a, int b, double err, double[] Va) {
    int n = Va.Length;
    double[][] result = new double[a-b+1][];
    for (int y=a; y>=b; y--) 
      result[y-b] = new double[n];
    Array.Copy(Va, result[a-b], Va.Length);

    double relerr = err, abserr = err;
    Rkf45 solver = new Rkf45(dV, n);
    double[] 
      v = new double[n], 
      tmp = new double[n];
    Array.Copy(Va, v, Va.Length);
    for (int y=a; y>b; y--) { 
      bj_ii(y, tmp);
      saxpy(1, tmp, v, v);  
      double t = y;
      // Integrate over [y,y+1]
      solver.r8_rkf45(v, ref t, t-1, ref relerr, abserr, 2);
      //      Console.WriteLine(y + " " + flag + " " + v[0]);
      Array.Copy(v, result[y-b-1], v.Length);
    }
    return result;
  }

  // saxpy = scalar a times x array plus y array, imperative version
  static void saxpy(double a, double[] x, double[] y, double[] res) {
    if (x.Length != y.Length)
      throw new Exception("saxpy: lengths of x and y differ");
    if (x.Length != res.Length)
      throw new Exception("saxpy: lengths of x and res differ");
    for (int i=0; i<x.Length; i++)
      res[i] = a * x[i] + y[i];
  }

}

public class Timer {
  private Stopwatch stopwatch;

  public Timer() {
    stopwatch = new Stopwatch();
    stopwatch.Reset();
    stopwatch.Start();
  }
  
  public double Check() {
    return stopwatch.ElapsedMilliseconds;
  }
}

/******************************************************************************/

class TestRkf45 {
  public static void Main(String[] args) {
    // Rkf45 solver = new Rkf45(DeferredTemporaryLifeAnnuity, 1);
    // double[] y = { 0.0 };
    // double relerr = 1E-10;

    // // Solve 
    // for (int year=50; year>=1; year-=1) {
    //   double t = year;
    //   int flag = solver.r8_rkf45(y, ref t, t-1, ref relerr, 1E-15, 2);
    //   Console.WriteLine(year + " " + flag + " " + y[0]);
    // }
    // Compute and print reserves

    //CalculationSpecifications.TimeAll(12288);
    CalculationSpecifications.ComputeAll();
  }

  public static void Zero(double t, double[] yv, /* for output: */ double[] yvp) {
    Debug.Assert(1 == yv.Length && yv.Length == yvp.Length);
    yvp[0] = 0;
  }

  public static void Two(double t, double[] yv, /* for output: */ double[] yvp) {
    Debug.Assert(1 == yv.Length && yv.Length == yvp.Length);
    yvp[0] = 2;
  }

  public static void Logistic(double t, double[] yv, /* for output: */ double[] yvp) {
    Debug.Assert(1 == yv.Length && yv.Length == yvp.Length);
    yvp[0] = 10 * yv[0] * (1 - yv[0]);
  }

  public static void DeferredTemporaryLifeAnnuity(double t, double[] yv, 
						  /* for output: */ double[] yvp) {
    Debug.Assert(1 == yv.Length && yv.Length == yvp.Length);
    yvp[0] = r(t) * yv[0] - b_0(t) - GM(t) * (0 - yv[0]);
  }

  static double b_0(double t) {
    return 1 * indicator(t > 35) * indicator(t < 35 + 10);
  }

  static double r(double t) { 
    return 0.05; 
  }

  static double indicator(bool b) {
    return b ? 1.0 : 0.0;
  }

  // Gompertz-Makeham mortality intensities for Danish women aged 30+t
  static double GM(double t) {
    return 0.0005 + Math.Pow(10, 5.728 - 10 + 0.038*(30 + t));
  }
}

// The weird encapsulation style (classes PureEndowment and so on) is
// motivated by the wish to stay relatively close to C, so no
// inheritance.  Also, while the function argument dV could be
// implemented in C using a function pointer, this is probably
// impossible or a bad idea on a GPU. 

class CalculationSpecifications {
  public static readonly double err = 1e-11;
  
  static readonly double age = 30,
    interestrate = 0.05,
    bpension = 1,
    pensiontime = 35;

  static double indicator(bool b) {
    return b ? 1.0 : 0.0;
  }

  // 

  public static void TimeAll(int customers) {
    Console.WriteLine("PureEndowment:               " + PureEndowment.Time(customers));
    Console.WriteLine("DeferredTemoraryLifeAnnuity: " + DeferredTemporaryLifeAnnuity.Time(customers));
    Console.WriteLine("TemporaryLifeAnnuityPremium: " + TemporaryLifeAnnuityPremium.Time(customers));
    Console.WriteLine("TermInsurance:               " + TermInsurance.Time(customers));
    Console.WriteLine("DisabilityAnnuity:           " + DisabilityAnnuity.Time(customers));
    Console.WriteLine("DisabilityTermInsurance:     " + DisabilityTermInsurance.Time(customers));
  }

  public static void ComputeAll() {
    // Compute and print reserves
    Print(PureEndowment.Compute());
    Print(DeferredTemporaryLifeAnnuity.Compute());
    Print(TemporaryLifeAnnuityPremium.Compute());
    Print(TermInsurance.Compute());
    Print(DisabilityAnnuity.Compute());
    Print(DisabilityTermInsurance.Compute());
  }

  // Gompertz-Makeham mortality intensities for Danish women
  static double GM(double t) {
    return 0.0005 + Math.Pow(10, 5.728 - 10 + 0.038*(age + t));
  }

  static double r(double t) { 
    return interestrate;    // Fixed interest rate
    //return rFsa(t);            // Finanstilsynet's rate curve
  }

  // The Danish FSA yield curve (Finanstilsynets rentekurve).
  // Data from 2011-11-16 
  static readonly double[] 
    ts = new double[] { 
      0.25, 0.5, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 
      15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30 },
    rs = new double[] { 
      1.146677033, 1.146677033, 1.146677033, 1.340669678, 1.571952911, 1.803236144, 
      2.034519377, 2.26580261, 2.497085843, 2.584085843, 2.710085843, 2.805085843, 
      2.871485843, 2.937885843, 3.004285843, 3.070685843, 3.137085843, 3.136485843, 
      3.135885843, 3.135285843, 3.134685843, 3.134085843, 3.113185843, 3.092285843, 
      3.071385843, 3.050485843, 3.029585843, 3.008685843, 2.987785843, 2.966885843, 
      2.945985843, 2.925085843
    };

  // Get discount rate at time t by linear interpolation into (ts,rs); then
  // compute the instantaneous forward rate as described in 
  // https://wiki.actulus.dk/Documentation-CalculationPlatform-YieldCurves.ashx

  // This method uses binary search, which is needlessly slow because
  // the t values are monotonically decreasing.  Hence it would be
  // faster to keep the most recent index m into the ts/rs arrays and
  // decrement m only when t < ts[m].  It would also be easier to get
  // wrong, because it relies on many assumptions, so for now we stick
  // to the binary search.
 
  static double rFsa(double t) { 
    // Requires ts non-empty and elements strictly increasing.
    int last = ts.Length-1;
    if (t <= ts[0])
      return Math.Log(1 + rs[0]/100);
    else if (t >= ts[last])
      return Math.Log(1 + rs[last]/100);
    else {
      int a = 0, b = last;
      // Now a < b (bcs. ts must have more than 1 element) and ts[a] < t < ts[b]
      while (a+1 < b) {
	// Now a < b and ts[a] <= t < ts[b]
	int i = (a+b)/2;
	if (ts[i] <= t)
	  a = i;
	else // t < ts[i]
	  b = i;
      }
      // Now a+1>=b and ts[a] <= t < ts[b]; so a!=b and hence a+1 == b <= last
      int m = a;
      double tm = ts[m], tm1 = ts[m+1];
      double rm = rs[m] / 100, rm1 = rs[m+1] / 100;
      double Rt = (rm * (tm1 - t) + rm1 * (t - tm)) / (tm1 - tm);
      return Math.Log(1 + Rt) + t / (tm1 - tm) * (rm1 - rm) / (1 + Rt);
    }
  }

  // Payment stream in state 0 (alive)
  static double b_0(double t) {
    return 0.0;
  }

  // Lump sum payments while in state 0 (alive) 
  static double bj_00(double t) {
    // This works only because t is known to be an integer
    return t == pensiontime ? bpension : 0.0;
  }

  // Lump sum payments while in state 1 (dead) 
  static double bj_11(double t) {
    // This works only because t is known to be an integer
    return 0.0;
  }

  // Transition intensity from state 0 (alive) to state 1 (dead)
  static double mu_01(double t) {
    return GM(t);
  }

  // Lump sum payment on transition from state 0 (alive) to state 1 (dead)
  static double bj_01(double t) {
    return 0.0;
  }

  // The dV0/dt function as used in Thiele:

  static double dV0(double t, double V0) {
    return r(t) * V0 - b_0(t) - mu_01(t) * (0 - V0 + bj_01(t));
  }

  static double dV0(double t, double V0, double V1) {
    return r(t) * V0 - b_0(t) - mu_01(t) * (0 - V0 + bj_01(t));
  }

  static double dV1(double t, double V0, double V1) {
    return 0.0;
  }

  static double[] dV(double t, double[] V) {
    double dV0 = r(t) * V[0] - b_0(t) - mu_01(t) * (V[1] - V[0] + bj_01(t));
    double dV1 = 0;
    return new double[] { dV0, dV1 };
  }

  static void Print(double[][] result) {
    for (int y=0; y<result.Length; y++) {
      Console.Write("{0,3}:", y);
      for (int i=0; i<result[y].Length; i++)
	Console.Write("  {0,20:F16}", result[y][i]);
      Console.WriteLine();
    }
  }

  // The two-state Actulus calculation kernel examples; 
  // really one-state because V1(t) = 0 for all t.

  public class PureEndowment {
    static double b_0(double t) {
      return 0.0;
    }
    
    static double mu_01(double t) {
      return GM(t);
    }
    
    static double bj_00(double t) {
      // This works only because t is known to be an integer
      return t == pensiontime ? bpension : 0.0;
    }

    static double bj_01(double t) {
      return 0.0;
    }

    public static double Time(int customers) {
      Timer timer = new Timer();
      double start_time = timer.Check();
      for(int i = 0;i<customers;i++)
        Compute();
      double end_time = timer.Check();
      return end_time - start_time;
    }

    public static double[][] Compute() {
      return Rkf45.RKF45_n((double t, double[] V, double[] res) =>
          { res[0] = r(t) * V[0] - b_0(t) - mu_01(t) * (0 - V[0] + bj_01(t)); },
          (double t, double[] res) => { res[0] = bj_00(t); },
          40, 0, err, new double[] { 0 });
    }
  }

  public class DeferredTemporaryLifeAnnuity {
    static int m = 35, n = 10;

    static double b_0(double t) {
      return bpension * indicator(t > m) * indicator(t < m + n);
    }

    static double mu_01(double t) {
      return GM(t);
    }

    static double bj_00(double t) {
      return 0.0;
    }

    static double bj_01(double t) {
      return 0.0;
    }

  public static double Time(int customers) {
    Timer timer = new Timer();
    double start_time = timer.Check();
    for(int i = 0;i<customers;i++)
      Compute();
    double end_time = timer.Check();
    return end_time - start_time;
  }
    public static double[][] Compute() {
      return Rkf45.RKF45_n((double t, double[] V, double[] res) =>
          { res[0] = r(t) * V[0] - b_0(t) - mu_01(t) * (0 - V[0] + bj_01(t)); },
          (double t, double[] res) => { res[0] = bj_00(t); },
          50, 0, err, new double[] { 0 });
    }
  }

  public class TemporaryLifeAnnuityPremium {
    static int n = 35;
    static double bpremium = 1;

    static double b_0(double t) {
      return -bpremium * indicator(t >= 0) * indicator(t < n);
    }

    static double mu_01(double t) {
      return GM(t);
    }

    static double bj_00(double t) {
      return 0.0;
    }

    static double bj_01(double t) {
      return 0.0;
    }
public static double Time(int customers) {
  Timer timer = new Timer();
  double start_time = timer.Check();
  for(int i = 0;i<customers;i++)
    Compute();
  double end_time = timer.Check();
  return end_time - start_time;
}

    public static double[][] Compute() {
      return Rkf45.RKF45_n((double t, double[] V, double[] res) =>
          { res[0] = r(t) * V[0] - b_0(t) - mu_01(t) * (0 - V[0] + bj_01(t)); },
          (double t, double[] res) => { res[0] = bj_00(t); },
          50, 0, err, new double[] { 0 });
    }
  }

  public class TermInsurance {
    static int n = 35;
    static double bdeath = 1;

    static double b_0(double t) {
      return 0.0;
    }

    static double mu_01(double t) {
      return GM(t);
    }

    static double bj_00(double t) {
      return 0.0;
    }

    static double bj_01(double t) {
      return bdeath * indicator(t > 0) * indicator(t < n);
    }

  public static double Time(int customers) {
    Timer timer = new Timer();
    double start_time = timer.Check();
    for(int i = 0;i<customers;i++)
      Compute();
    double end_time = timer.Check();
    return end_time - start_time;
  }
    public static double[][] Compute() {
      return Rkf45.RKF45_n((double t, double[] V, double[] res) =>
          { res[0] = r(t) * V[0] - b_0(t) - mu_01(t) * (0 - V[0] + bj_01(t)); },
          (double t, double[] res) => { res[0] = bj_00(t); },
          50, 0, err, new double[] { 0 });
    }
  }

  // The three-state Actulus calculation kernel examples; 
  // really two-state because V2(t) = 0 for all t.

  public class DisabilityAnnuity {
    static int n = 35;
    static double bdisabled = 1;

    static double b_0(double t) {
      return 0.0;
    }

    static double b_1(double t) {
      return bdisabled * indicator(t > 0) * indicator(t < n);
    }

    static double GM01(double t) {
      return 0.0006 + Math.Pow(10, 4.71609 - 10 + 0.06*(age + t));
    }

    static double GM02(double t) {
      return GM(t);
    }

    static double GM12(double t) {
      return GM(t);
    }

    static double mu_01(double t) {
      return GM01(t);
    }

    static double mu_02(double t) {
      return GM02(t);
    }

    static double mu_12(double t) {
      return GM12(t);
    }

    static double bj_00(double t) {
      return 0.0;
    }

    static double bj_01(double t) {
      return 0.0;
    }

    static double bj_02(double t) {
      return 0.0;
    }

    static double bj_11(double t) {
      return 0.0;
    }

    static double bj_12(double t) {
      return 0.0;
    }

  public static double Time(int customers) {
    Timer timer = new Timer();
    double start_time = timer.Check();
    for(int i = 0;i<customers;i++)
      Compute();
    double end_time = timer.Check();
    return end_time - start_time;
  }
    public static double[][] Compute() {
      return Rkf45.RKF45_n((double t, double[] V, double[] res) =>
          { res[0] = r(t) * V[0] - b_0(t) - mu_01(t) * (V[1] - V[0] + bj_01(t)) 
          - mu_02(t) * (0 - V[0] + bj_02(t));
          res[1] = r(t) * V[1] - b_1(t) - mu_12(t) * (0 - V[1] + bj_12(t)); },
          (double t, double[] res) => { res[0] = bj_00(t); res[1] = bj_11(t); },
          50, 0, err, new double[] { 0, 0 });
    }
  }

  public class DisabilityTermInsurance {
    static int n = 35;
    static double bdisabled = 1;

    static double b_0(double t) {
      return 0.0;
    }

    static double b_1(double t) {
      return 0.0;
    }

    static double GM01(double t) {
      return 0.0006 + Math.Pow(10, 4.71609 - 10 + 0.06*(age + t));
    }

    static double GM02(double t) {
      return GM(t);
    }

    static double GM12(double t) {
      return GM(t);
    }

    static double mu_01(double t) {
      return GM01(t);
    }

    static double mu_02(double t) {
      return GM02(t);
    }

    static double mu_12(double t) {
      return GM12(t);
    }

    static double bj_00(double t) {
      return 0.0;
    }

    static double bj_01(double t) {
      return bdisabled * indicator(t > 0) * indicator(t < n);
    }

    static double bj_02(double t) {
      return 0.0;
    }

    static double bj_11(double t) {
      return 0.0;
    }    

    static double bj_12(double t) {
      return 0.0;
    }

  public static double Time(int customers) {
    Timer timer = new Timer();
    double start_time = timer.Check();
    for(int i = 0;i<customers;i++)
      Compute();
    double end_time = timer.Check();
    return end_time - start_time;
  }
    public static double[][] Compute() {
      return Rkf45.RKF45_n((double t, double[] V, double[] res) =>
          { res[0] = r(t) * V[0] - b_0(t) - mu_01(t) * (V[1] - V[0] + bj_01(t)) 
          - mu_02(t) * (0 - V[0] + bj_02(t));
          res[1] = r(t) * V[1] - b_1(t) - mu_12(t) * (0 - V[1] + bj_12(t)); },
          (double t, double[] res) => { res[0] = bj_00(t); res[1] = bj_11(t); },
          50, 0, err, new double[] { 0, 0 });
    }
  }
}
