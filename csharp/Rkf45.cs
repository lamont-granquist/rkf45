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

    /*
     * TODO:
     * A Description of what is going on with each variable. 
     * Some Assertions?
     * Simplification?
     */

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
  {
  /* TODO:
   * Describe each variable. Delete some of them.
   *
   *
   *
   *
   */
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

  //Make a n*[a-b+1] matrix
  int n = Va.Length;
  double[][] result = new double[a-b+1][];
  for (int y=a; y>=b; y--) 
    result[y-b] = new double[n];
  Array.Copy(Va, result[a-b], Va.Length); // Insert start values

  double relerr = err, abserr = err; // Set errors

  Rkf45 solver = new Rkf45(dV, n); //Make solver

  // Va = v
  double[] v = new double[n], tmp = new double[n];
  Array.Copy(Va, v, Va.Length);

  //Solve for one year at a time
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
    return stopwatch.ElapsedMilliseconds / 1000.0;
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

  static bool IsEqual(double[][] a,double [][] b) {
    for (int i=0; i<a.Length; i++) {
    //Console.WriteLine(a[i].Length);  
    //Console.WriteLine(b[i].Length);  
      for (int j=0; j<a[i].Length; j++) {
        if (Math.Abs(a[i][j] - b[i][j]) > err)
          return false;
      }
    }
    return true;
  }

  static void Assert(bool b,string msg="Fail!") {
    if (b == false) {
      throw new Exception(msg);
    }
  }

  public static void ComputeAll() {
    // Compute and print reserves
    Assert(IsEqual(PureEndowment.Compute(),PureEndowment.test_values),"PureEndowment failed");
    Assert(IsEqual(DeferredTemporaryLifeAnnuity.Compute(),DeferredTemporaryLifeAnnuity.test_values),"DeferredTemporaryLifeAnnuity failed");
    Assert(IsEqual(TemporaryLifeAnnuityPremium.Compute(),TemporaryLifeAnnuityPremium.test_values),"TempLifeAnnuPrem failed");
    Assert(IsEqual(TermInsurance.Compute(),TermInsurance.test_values),"TempInsurance failed");
    Assert(IsEqual(DisabilityAnnuity.Compute(),DisabilityAnnuity.test_values),"DisAnnu failed");
    Assert(IsEqual(DisabilityTermInsurance.Compute(),DisabilityTermInsurance.test_values),"DisabilityTermInsurance failed");
    //Print(PureEndowment.Compute());
    //Print(DeferredTemporaryLifeAnnuity.Compute());
    //Print(TemporaryLifeAnnuityPremium.Compute());
    //Print(TermInsurance.Compute());
    //Print(DisabilityAnnuity.Compute());
    //Print(DisabilityTermInsurance.Compute());
  }

  // Gompertz-Makeham mortality intensities for Danish women
  static double GM(double t) {
    return 0.0005 + Math.Pow(10, 5.728 - 10 + 0.038*(age + t));
  }

  static double r(double t) { 
    return interestrate;    // Fixed interest rate
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

  //Print
  static void Print(double[][] result) {
    for (int y=0; y<result.Length; y++) {
      Console.Write("{0,3}:", y);
      for (int i=0; i<result[y].Length; i++)
        Console.Write("  {0,20:F16}", result[y][i]);
      Console.WriteLine();
    }
  }

  //Temp method
  static void Print_test_values(double[][] result) {
    for (int y=0; y<result.Length; y++) {
      //Console.Write("{0,3}:", y);
      Console.Write("new double[] {");
      for (int i=0; i<result[y].Length; i++) {
        Console.Write("  {0,20:F16}", result[y][i]);
        if (i<result[y].Length-1)
          Console.Write(",");
      }
      Console.Write("},");
      Console.WriteLine();
    }
  }

  // The two-state Actulus calculation kernel examples; 
  // really one-state because V1(t) = 0 for all t.

  public class PureEndowment {

    static public double[][] test_values = new double[][] {
new double[] {0.1437946974886250},
new double[] {0.1513594875720590},
new double[] {0.1593334830357050},
new double[] {0.1677404776222450},
new double[] {0.1766058889630630},
new double[] {0.1859569023966960},
new double[] {0.1958226315281160},
new double[] {0.2062342978884670},
new double[] {0.2172254324285080},
new double[] {0.2288321020171980},
new double[] {0.2410931646317720},
new double[] {0.2540505575320670},
new double[] {0.2677496234274150},
new double[] {0.2822394804908960},
new double[] {0.2975734430792770},
new double[] {0.3138095012097790},
new double[] {0.3310108682663330},
new double[] {0.3492466081065300},
new double[] {0.3685923547759590},
new double[] {0.3891311404830350},
new double[] {0.4109543504366150},
new double[] {0.4341628267168120},
new double[] {0.4588681476758340},
new double[] {0.4851941146396860},
new double[] {0.5132784841210240},
new double[] {0.5432749916555080},
new double[] {0.5753557231029890},
new double[] {0.6097139012822690},
new double[] {0.6465671707381610},
new double[] {0.6861614820516400},
new double[] {0.7287757004081980},
new double[] {0.7747270924521970},
new double[] {0.8243778824987340},
new double[] {0.8781431162166330},
new double[] {0.9365001299367070},
new double[] {0.0000000000000000},
new double[] {0.0000000000000000},
new double[] {0.0000000000000000},
new double[] {0.0000000000000000},
new double[] {0.0000000000000000},
new double[] {0.0000000000000000}
   };

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

    public static double[][] Compute() {
      //Console.WriteLine("\n PureEndowment");
      //Print(new double[][] { new double[] { 0.14379469738 } });
      //Console.WriteLine();
      return Rkf45.RKF45_n((double t, double[] V, double[] res) =>
          { res[0] = r(t) * V[0] - b_0(t) - mu_01(t) * (0 - V[0] + bj_01(t)); },
          (double t, double[] res) => { res[0] = bj_00(t); },
          40, 0, err, new double[] { 0 });
    }
  }

  public class DeferredTemporaryLifeAnnuity {

  static public double[][] test_values = new double[][] {
    new double[] {    1.0265607676014400},
    new double[] {    1.0805663523022000},
    new double[] {    1.1374932838714300},
    new double[] {    1.1975114275626800},
    new double[] {    1.2608022416886800},
    new double[] {    1.3275598043516300},
    new double[] {    1.3979919596875600},
    new double[] {    1.4723216004702400},
    new double[] {    1.5507881065874600},
    new double[] {    1.6336489620308100},
    new double[] {    1.7211815767162100},
    new double[] {    1.8136853437812100},
    new double[] {    1.9114839681141600},
    new double[] {    2.0149281079127500},
    new double[] {    2.1243983782341300},
    new double[] {    2.2403087740143500},
    new double[] {    2.3631105801829600},
    new double[] {    2.4932968486250200},
    new double[] {    2.6314075362739100},
    new double[] {    2.7780354160836500},
    new double[] {    2.9338328936850100},
    new double[] {    3.0995198879959800},
    new double[] {    3.2758929649601700},
    new double[] {    3.4638359512170500},
    new double[] {    3.6643323004950200},
    new double[] {    3.8784795419259000},
    new double[] {    4.1075062089359300},
    new double[] {    4.3527917332333000},
    new double[] {    4.6158898949982800},
    new double[] {    4.8985565532549200},
    new double[] {    5.2027825467745500},
    new double[] {    5.5308328651267400},
    new double[] {    5.8852934539509200},
    new double[] {    6.2691273543593000},
    new double[] {    6.6857423050148000},
    new double[] {    7.1390724853580000},
    new double[] {    6.5992994744162900},
    new double[] {    6.0319762072565800},
    new double[] {    5.4341959778729600},
    new double[] {    4.8025044157278400},
    new double[] {    4.1327786917263300},
    new double[] {    3.4200771644266400},
    new double[] {    2.6584512106845400},
    new double[] {    1.8407083534208200},
    new double[] {    0.9581122236639260},
    new double[] {    0.0000000000000000},
    new double[] {    0.0000000000000000},
    new double[] {    0.0000000000000000},
    new double[] {    0.0000000000000000},
    new double[] {    0.0000000000000000},
    new double[] {    0.0000000000000000}
    };

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

    public static double[][] Compute() {
      //Console.WriteLine("\n DeferredTemporaryLifeAnnuity");
      //Print(new double[][] { new double[] { 1.0265607675 } });
      //Console.WriteLine();
      return Rkf45.RKF45_n((double t, double[] V, double[] res) =>
          { res[0] = r(t) * V[0] - b_0(t) - mu_01(t) * (0 - V[0] + bj_01(t)); },
          (double t, double[] res) => { res[0] = bj_00(t); },
          50, 0, err, new double[] { 0 });
    }
  }

  public class TemporaryLifeAnnuityPremium {

  static public double[][] test_values = new double[][] {
    new double[] {  -15.9717676660001000},
    new double[] {  -15.7859295725898000},
    new double[] {  -15.5914495774420000},
    new double[] {  -15.3879467041606000},
    new double[] {  -15.1750230434772000},
    new double[] {  -14.9522626687192000},
    new double[] {  -14.7192304105538000},
    new double[] {  -14.4754704664848000},
    new double[] {  -14.2205048162637000},
    new double[] {  -13.9538314093213000},
    new double[] {  -13.6749220843743000},
    new double[] {  -13.3832201743607000},
    new double[] {  -13.0781377416068000},
    new double[] {  -12.7590523783886000},
    new double[] {  -12.4253034965330000},
    new double[] {  -12.0761880160465000},
    new double[] {  -11.7109553465719000},
    new double[] {  -11.3288015361432000},
    new double[] {  -10.9288624386922000},
    new double[] {  -10.5102057241661000},
    new double[] {  -10.0718215219905000},
    new double[] {   -9.6126114486954800},
    new double[] {   -9.1313757222581100},
    new double[] {   -8.6267980071471600},
    new double[] {   -8.0974275627022600},
    new double[] {   -7.5416581802140900},
    new double[] {   -6.9577032868934500},
    new double[] {   -6.3435664627206500},
    new double[] {   -5.6970064523888100},
    new double[] {   -5.0154955507238000},
    new double[] {   -4.2961699850952200},
    new double[] {   -3.5357705981019300},
    new double[] {   -2.7305717294876500},
    new double[] {   -1.8762956830915000},
    new double[] {   -0.9680095100191570},
    new double[] {    0.0000000000000000},
    new double[] {    0.0000000000000000},
    new double[] {    0.0000000000000000},
    new double[] {    0.0000000000000000},
    new double[] {    0.0000000000000000},
    new double[] {    0.0000000000000000},
    new double[] {    0.0000000000000000},
    new double[] {    0.0000000000000000},
    new double[] {    0.0000000000000000},
    new double[] {    0.0000000000000000},
    new double[] {    0.0000000000000000},
    new double[] {    0.0000000000000000},
    new double[] {    0.0000000000000000},
    new double[] {    0.0000000000000000},
    new double[] {    0.0000000000000000},
    new double[] {    0.0000000000000000}
  };

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

    public static double[][] Compute() {
      //Console.WriteLine("\n TemporaryLifeAnnuityPremium");
      //Print(new double[][] { new double[] { -15.971767666 } });
      //Console.WriteLine();
      return Rkf45.RKF45_n((double t, double[] V, double[] res) =>
          { res[0] = r(t) * V[0] - b_0(t) - mu_01(t) * (0 - V[0] + bj_01(t)); },
          (double t, double[] res) => { res[0] = bj_00(t); },
          50, 0, err, new double[] { 0 });
    }
  }

  public class TermInsurance {
    
  static public double[][] test_values = new double[][] {
        new double[] {    0.0576169193132673},
        new double[] {    0.0593440338503396},
        new double[] {    0.0610940381466785},
        new double[] {    0.0628621872268886},
        new double[] {    0.0646429589229975},
        new double[] {    0.0664299642301075},
        new double[] {    0.0682158480098718},
        new double[] {    0.0699921788564634},
        new double[] {    0.0717493268311646},
        new double[] {    0.0734763275934875},
        new double[] {    0.0751607312303756},
        new double[] {    0.0767884338351089},
        new double[] {    0.0783434895820515},
        new double[] {    0.0798079006843388},
        new double[] {    0.0811613821939618},
        new double[] {    0.0823810980933246},
        new double[] {    0.0834413645163828},
        new double[] {    0.0843133152038705},
        new double[] {    0.0849645234136384},
        new double[] {    0.0853585734399381},
        new double[] {    0.0854545736024947},
        new double[] {    0.0852066009948634},
        new double[] {    0.0845630663660149},
        new double[] {    0.0834659851665421},
        new double[] {    0.0818501379168519},
        new double[] {    0.0796420995167974},
        new double[] {    0.0767591127460391},
        new double[] {    0.0731077757869581},
        new double[] {    0.0685825068600615},
        new double[] {    0.0630637406431617},
        new double[] {    0.0564158005824136},
        new double[] {    0.0484843779035996},
        new double[] {    0.0390935313045591},
        new double[] {    0.0280420999246517},
        new double[] {    0.0150993948779494},
        new double[] {    0.0000000000000000},
        new double[] {    0.0000000000000000},
        new double[] {    0.0000000000000000},
        new double[] {    0.0000000000000000},
        new double[] {    0.0000000000000000},
        new double[] {    0.0000000000000000},
        new double[] {    0.0000000000000000},
        new double[] {    0.0000000000000000},
        new double[] {    0.0000000000000000},
        new double[] {    0.0000000000000000},
        new double[] {    0.0000000000000000},
        new double[] {    0.0000000000000000},
        new double[] {    0.0000000000000000},
        new double[] {    0.0000000000000000},
        new double[] {    0.0000000000000000},
        new double[] {    0.0000000000000000}
  };

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

    public static double[][] Compute() {
      //Console.WriteLine("\n TermInsurance");
      //Print(new double[][] { new double[] { 0.057616919318 } });
      //Console.WriteLine();
      return Rkf45.RKF45_n((double t, double[] V, double[] res) =>
          { res[0] = r(t) * V[0] - b_0(t) - mu_01(t) * (0 - V[0] + bj_01(t)); },
          (double t, double[] res) => { res[0] = bj_00(t); },
          50, 0, err, new double[] { 0 });
    }
  }

  // The three-state Actulus calculation kernel examples; 
  // really two-state because V2(t) = 0 for all t.

  public class DisabilityAnnuity {

  static public double[][] test_values = new double[][] {
    new double[] {    0.5555261079604120,   15.9717676673750000},
    new double[] {    0.5697939362470290,   15.7859295725873000},
    new double[] {    0.5842458860490700,   15.5914495774393000},
    new double[] {    0.5988112559939020,   15.3879467041578000},
    new double[] {    0.6134061668653150,   15.1750230434743000},
    new double[] {    0.6279320187147260,   14.9522626687161000},
    new double[] {    0.6422738476905750,   14.7192304105506000},
    new double[] {    0.6562985942478250,   14.4754704664814000},
    new double[] {    0.6698533011953220,   14.2205048162601000},
    new double[] {    0.6827632688462940,   13.9538314093175000},
    new double[] {    0.6948302058887720,   13.6749220843703000},
    new double[] {    0.7058304291697920,   13.3832201743564000},
    new double[] {    0.7155131842695250,   13.0781377416023000},
    new double[] {    0.7235991826741960,   12.7590523783855000},
    new double[] {    0.7297794820426480,   12.4253034965300000},
    new double[] {    0.7337148754958810,   12.0761880160438000},
    new double[] {    0.7350360067043140,   11.7109553465694000},
    new double[] {    0.7333444934112320,   11.3288015361409000},
    new double[] {    0.7282154278216300,   10.9288624386901000},
    new double[] {    0.7192017347676550,   10.5102057241642000},
    new double[] {    0.7058410171284570,   10.0718215219888000},
    new double[] {    0.6876657158018430,    9.6126114486939200},
    new double[] {    0.6642176772198330,    9.1313757222565700},
    new double[] {    0.6350685815395070,    8.6267980071455300},
    new double[] {    0.5998481774680660,    8.0974275627005400},
    new double[] {    0.5582829507417500,    7.5416581802122700},
    new double[] {    0.5102488039940220,    6.9577032868915200},
    new double[] {    0.4558426666237880,    6.3435664627186100},
    new double[] {    0.3954798644641650,    5.6970064523866500},
    new double[] {    0.3300268327704200,    5.0154955507223200},
    new double[] {    0.2609827682846260,    4.2961699850940100},
    new double[] {    0.1907297303473690,    3.5357705981010400},
    new double[] {    0.1228795253423570,    2.7305717294870500},
    new double[] {    0.0627590447157113,    1.8762956830912200},
    new double[] {    0.0180961562710709,    0.9680095100190450},
    new double[] {    0.0000000000000000,    0.0000000000000000},
    new double[] {    0.0000000000000000,    0.0000000000000000},
    new double[] {    0.0000000000000000,    0.0000000000000000},
    new double[] {    0.0000000000000000,    0.0000000000000000},
    new double[] {    0.0000000000000000,    0.0000000000000000},
    new double[] {    0.0000000000000000,    0.0000000000000000},
    new double[] {    0.0000000000000000,    0.0000000000000000},
    new double[] {    0.0000000000000000,    0.0000000000000000},
    new double[] {    0.0000000000000000,    0.0000000000000000},
    new double[] {    0.0000000000000000,    0.0000000000000000},
    new double[] {    0.0000000000000000,    0.0000000000000000},
    new double[] {    0.0000000000000000,    0.0000000000000000},
    new double[] {    0.0000000000000000,    0.0000000000000000},
    new double[] {    0.0000000000000000,    0.0000000000000000},
    new double[] {    0.0000000000000000,    0.0000000000000000},
    new double[] {    0.0000000000000000,    0.0000000000000000}
    };

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

    public static double[][] Compute() {
      //Console.WriteLine("\n DisabilityAnnuity");
      //Print(new double[][] { new double[] { 0.55552610797, 15.971767666 } });
      //Console.WriteLine();
      return Rkf45.RKF45_n((double t, double[] V, double[] res) =>
          { res[0] = r(t) * V[0] - b_0(t) - mu_01(t) * (V[1] - V[0] + bj_01(t)) 
          - mu_02(t) * (0 - V[0] + bj_02(t));
          res[1] = r(t) * V[1] - b_1(t) - mu_12(t) * (0 - V[1] + bj_12(t)); },
          (double t, double[] res) => { res[0] = bj_00(t); res[1] = bj_11(t); },
          50, 0, err, new double[] { 0, 0 });
    }
  }

  public class DisabilityTermInsurance {

  static public double[][] test_values = new double[][] {
    new double[] {    0.0714186989824431,    0.0000000000000000},
    new double[] {    0.0742705084048387,    0.0000000000000000},
    new double[] {    0.0772312485352858,    0.0000000000000000},
    new double[] {    0.0803007388769112,    0.0000000000000000},
    new double[] {    0.0834779380618836,    0.0000000000000000},
    new double[] {    0.0867607747788289,    0.0000000000000000},
    new double[] {    0.0901459511899878,    0.0000000000000000},
    new double[] {    0.0936287143026796,    0.0000000000000000},
    new double[] {    0.0972025899462083,    0.0000000000000000},
    new double[] {    0.1008590729874850,    0.0000000000000000},
    new double[] {    0.1045872661636680,    0.0000000000000000},
    new double[] {    0.1083734583544920,    0.0000000000000000},
    new double[] {    0.1122006311731980,    0.0000000000000000},
    new double[] {    0.1160478803013540,    0.0000000000000000},
    new double[] {    0.1198897349010070,    0.0000000000000000},
    new double[] {    0.1236953544561090,    0.0000000000000000},
    new double[] {    0.1274275773038030,    0.0000000000000000},
    new double[] {    0.1310417884996260,    0.0000000000000000},
    new double[] {    0.1344845660310780,    0.0000000000000000},
    new double[] {    0.1376920530492920,    0.0000000000000000},
    new double[] {    0.1405879887686990,    0.0000000000000000},
    new double[] {    0.1430813106541400,    0.0000000000000000},
    new double[] {    0.1450632136041590,    0.0000000000000000},
    new double[] {    0.1464035154086170,    0.0000000000000000},
    new double[] {    0.1469461280481600,    0.0000000000000000},
    new double[] {    0.1465033660147020,    0.0000000000000000},
    new double[] {    0.1448487279175560,    0.0000000000000000},
    new double[] {    0.1417076547118570,    0.0000000000000000},
    new double[] {    0.1367455798907500,    0.0000000000000000},
    new double[] {    0.1295523183539230,    0.0000000000000000},
    new double[] {    0.1196214525705190,    0.0000000000000000},
    new double[] {    0.1063228073533260,    0.0000000000000000},
    new double[] {    0.0888652648692060,    0.0000000000000000},
    new double[] {    0.0662459119488600,    0.0000000000000000},
    new double[] {    0.0371795952968856,    0.0000000000000000},
    new double[] {    0.0000000000000000,    0.0000000000000000},
    new double[] {    0.0000000000000000,    0.0000000000000000},
    new double[] {    0.0000000000000000,    0.0000000000000000},
    new double[] {    0.0000000000000000,    0.0000000000000000},
    new double[] {    0.0000000000000000,    0.0000000000000000},
    new double[] {    0.0000000000000000,    0.0000000000000000},
    new double[] {    0.0000000000000000,    0.0000000000000000},
    new double[] {    0.0000000000000000,    0.0000000000000000},
    new double[] {    0.0000000000000000,    0.0000000000000000},
    new double[] {    0.0000000000000000,    0.0000000000000000},
    new double[] {    0.0000000000000000,    0.0000000000000000},
    new double[] {    0.0000000000000000,    0.0000000000000000},
    new double[] {    0.0000000000000000,    0.0000000000000000},
    new double[] {    0.0000000000000000,    0.0000000000000000},
    new double[] {    0.0000000000000000,    0.0000000000000000},
    new double[] {    0.0000000000000000,    0.0000000000000000}
  };

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

    public static double[][] Compute() {
      //Console.WriteLine("\n DisabilityTermInsurance");
      //Print(new double[][] { new double[] { 0.071418699003, 0.000000000 } });
      //Console.WriteLine();
      return Rkf45.RKF45_n((double t, double[] V, double[] res) =>
          { res[0] = r(t) * V[0] - b_0(t) - mu_01(t) * (V[1] - V[0] + bj_01(t)) 
          - mu_02(t) * (0 - V[0] + bj_02(t));
          res[1] = r(t) * V[1] - b_1(t) - mu_12(t) * (0 - V[1] + bj_12(t)); },
          (double t, double[] res) => { res[0] = bj_00(t); res[1] = bj_11(t); },
          50, 0, err, new double[] { 0, 0 });
    }
  }
}
