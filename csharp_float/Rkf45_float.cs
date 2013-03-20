/*
 * Modern Rkf45. Made by simplification and refactoring of original code
 */

using System;
using System.Diagnostics;

class Estimator {

  //Public (to be set before estimation)
  private int neqn; //set in constructor
  public int start_year;
  public int end_year;
  public Action<float, float[], float[]> dy;
  public Action<float,float[]> bj_ii;
  public float relerr;
  public float abserr;
  public float [] end_year_y;

  //Private
  private float t;
  private float stepsize;
  private float[] f1;
  private float[] f2;
  private float[] f3;
  private float[] f4;
  private float[] f5;
  private float[] f_swap;
  private float[] y;
  private float[] y_diff;
  private float[] y_plus_one;
  private float[] y_plus_one_alternative;
  private int local_start_year; 
  private int local_end_year; 
  static float DoubleEpsilon; //Const used for calculations

  public int steps_taken_in_last_estimation;

  /***** TESTED *****/
  public int TESTED_NUMBER = 0;
  public int TESTED_MAX    = 5;

  public void TESTED() {
    TESTED_NUMBER++;
    if (TESTED_NUMBER > TESTED_MAX) {
      Environment.Exit(0);
    }
  }


  /************************** Constructor ***********************/

  /* Construct */
  public void construct(int this_neqn) {
    neqn = this_neqn;
    DoubleEpsilon = FindDoubleEpsilon();
    allocate_equation_space();
  }

  /* Allocate equation space */
  private void allocate_equation_space() {
    //Global for the class 
    y_plus_one = new float[neqn];
    end_year_y = new float[neqn];
    y = new float[neqn];
    y_diff = new float[neqn];

    //Temporary for the solve method
    f1 = new float[neqn];
    f2 = new float[neqn];
    f3 = new float[neqn];
    f4 = new float[neqn];
    f5 = new float[neqn];
    f_swap = new float[neqn];
    y_plus_one_alternative = new float[neqn];
  }

  /************************** Solve ***********************/

  /* Calculate the actual and the alternative solutions */
  //y_plus_one and y_plus_one_alternative will be set
  private void calculate_solutions() {

    float lcd_stepsize = stepsize / 4.0f; //lowest common denominator of stepsize

    //f1
    for (int i = 0; i < neqn; i++ )
      f_swap[i] = y[i] + lcd_stepsize * y_diff[i];
    dy ( t + lcd_stepsize, f_swap, f1 );

    //f2
    lcd_stepsize = 3.0f * stepsize / 32.0f;
    for (int i = 0; i < neqn; i++ )
      f_swap[i] = y[i] + lcd_stepsize * ( y_diff[i] + 3.0f * f1[i] );
    dy ( t + 3.0f * stepsize / 8.0f, f_swap, f2 );

    /*Console.WriteLine("f_swap!:             " + f_swap[0]);
    Console.WriteLine("!:             " + (t + 3.0f * stepsize / 8.0f));*/
    
    //f3
    lcd_stepsize = stepsize / 2197.0f;
    for (int i = 0; i < neqn; i++ )
      f_swap[i] = y[i] + lcd_stepsize * ( 1932.0f * y_diff[i] + ( 7296.0f * f2[i] - 7200.0f * f1[i] ) );
    dy ( t + 12.0f * stepsize / 13.0f, f_swap, f3 );

    //f4
    lcd_stepsize = stepsize / 4104.0f;
    for (int i = 0; i < neqn; i++ )
      f_swap[i] = y[i] + lcd_stepsize * ( ( 8341.0f * y_diff[i] - 845.0f * f3[i] ) + 
          ( 29440.0f * f2[i] - 32832.0f * f1[i] ) );
    dy ( t + stepsize, f_swap, f4 );

    //f5
    lcd_stepsize = stepsize / 20520.0f;
    for (int i = 0; i < neqn; i++ )
      f_swap[i] = y[i] + lcd_stepsize * ( ( -6080.0f * y_diff[i] + 
            ( 9295.0f * f3[i] - 5643.0f * f4[i] ) ) + ( 41040.0f * f1[i] - 28352.0f * f2[i] ) );
    dy ( t + stepsize / 2.0f, f_swap, f5 );

    //Calculate solution
    lcd_stepsize = stepsize / 7618050.0f;
    for (int i = 0; i < neqn; i++ )
      y_plus_one[i] = y[i] + lcd_stepsize * ( ( 902880.0f * y_diff[i] + 
            ( 3855735.0f * f3[i] - 1371249.0f * f4[i] ) ) + ( 3953664.0f * f2[i] + 277020.0f * f5[i] ) );

    //Calculate alternative solution
    for (int i = 0; i < neqn; i++ )
      y_plus_one_alternative[i] = ( -2090.0f * y_diff[i] + ( 21970.0f * f3[i] - 15048.0f * f4[i] ) ) + ( 22528.0f * f2[i] - 27360.0f * f5[i] );


      /*Console.WriteLine("y              " + y[0]);
      Console.WriteLine("y_diff         " + y_diff[0]);
      Console.WriteLine("t+lcd_stepsize " + (t+lcd_stepsize));
      Console.WriteLine("f1[0]:         " + f1[0]);
      Console.WriteLine("f2[0]:         " + f2[0]);
      Console.WriteLine("f3[0]:         " + f3[0]);
      Console.WriteLine("f4[0]:         " + f4[0]);
      Console.WriteLine("f5[0]:         " + f5[0]);
      Console.WriteLine("s [0]:         " + y_plus_one[0]);
      Console.WriteLine("sa[0]:         " + y_plus_one_alternative[0]);
      Console.WriteLine("");*/
      //TESTED();
  }

  /* Calculate the error of the solution */
  //Pure
  private float calculate_solution_error() {

    //Used in calculations
    float scale = 2.0f / relerr;

    //Calculate the biggest_difference
    float biggest_difference = 0.0f;
    for (int i = 0; i < neqn; i++ )
    {
      float et = Math.Abs( y[i] ) + Math.Abs( y_plus_one[i] ) + scale * abserr;
      float ee = Math.Abs( y_plus_one_alternative[i] );

      biggest_difference = Math.Max ( biggest_difference, ee / et );
    }

    //Return the error 
    return Math.Abs( stepsize ) * biggest_difference * scale / 752400.0f;
  }

  /******************* Local estimation ***********************/

  /* Move from current position to local_start_year, and update all values */
  // Updates y, h
  private void local_estimate() {
    t = (float)local_end_year;

    //Step by step integration.
    bool local_start_year_reached = false;
    while (!local_start_year_reached)
    {
      steps_taken_in_last_estimation++;

      //Variables used in calculations
      bool stepsize_decreased = false;
      float minimum_stepsize = 26.0f * DoubleEpsilon * Math.Abs( t );

      local_start_year_reached = local_start_year_to_be_reached(); //Has side effects for stepsize.

      //Try to solve
      calculate_solutions();
      float error = calculate_solution_error();

      //Integreate 1 step
      while(error > 1.0f)
      {
        local_start_year_reached = false;

        //Scale down.
        stepsize *= Math.Max(0.1f,0.9f / (float) Math.Round(Math.Pow( error, 0.2f ),7));
        stepsize_decreased = true;

        //Try again.
        calculate_solutions();
        error = calculate_solution_error();
      }
      
      //Apply solution
      for (int i = 0; i < neqn; i++ )
        y[i] = y_plus_one[i];

      //Advance in time
      t += stepsize; 

      //Update y_diff
      dy ( t, y, y_diff );

      //Apply scale to stepsize
      float scale = scale_from_error(error,stepsize_decreased);
      stepsize = Math.Sign ( stepsize ) * Math.Max ( scale * Math.Abs( stepsize ), minimum_stepsize );
    }
  }

  /******************* Local estimation help functions *******************/

  /* React if the "local start year" is about to be reached */
  //Effects stepsize, returns whether the start year is reached
  private bool local_start_year_to_be_reached() {
    float dt = local_start_year - t;
    if ( 2.0f * Math.Abs( stepsize ) > Math.Abs( dt ) )
    {
      if ( Math.Abs( dt ) <= Math.Abs( stepsize ) ) //Final step?
      {
        stepsize = dt;                   //Let stepsize hit output point
        return true;
      }
      else
      {
        stepsize = 0.5f * dt; // If not final step, set stepsize to be second final step. (evens out)
      }
    }
    return false;
  }

  /* Calculate stepsize's startvalue */
  public float calculate_initial_stepsize()
  {
    //Calculate the start value of stepsize
    float stepsize = Math.Abs( start_year - t );

    for (int k = 0; k < neqn; k++ )
    {
      float tol = relerr * Math.Abs( y[k] ) + abserr;
      if ( 0.0f < tol )
      {
        float y_diff_k = Math.Abs( y_diff[k] );
        if ( tol < y_diff_k * (float) Math.Round(Math.Pow( stepsize, 5 ),7) ) //We never really get in here...
        {
          stepsize = (float) Math.Round(Math.Pow( ( tol / y_diff_k ), 0.2f ),7);
          throw new Exception("initial stepsize calculation gone wrong");
        }
      }
    }

    return  Math.Max ( stepsize, 26.0f * DoubleEpsilon * Math.Max ( Math.Abs( t ), Math.Abs( start_year - t ) ) );
  }

  /* Scale from error calculations */
  public float scale_from_error(float error,bool stepsize_decreased) {
    float scale = Math.Min(5.0f,0.9f / (float) Math.Round(Math.Pow( error, 0.2f ),7));

    if (stepsize_decreased)
      scale = Math.Min ( scale, 1.0f );

    return scale;
  }

  /*************************** Estimate ***************************/

  public float[][] estimate()
  {
    steps_taken_in_last_estimation = 0;

    //Set the initial values
    Array.Copy(end_year_y,y,y.Length);       // y
    t = (float) end_year;                   // t
    dy ( t, y, y_diff );                     // y_diff
    stepsize = calculate_initial_stepsize(); // stepsize

    float[][] result = new float[end_year-start_year+1][];
    for (int year=end_year; year>=start_year; year--) 
      result[year-start_year] = new float[neqn];

    //Solve for one year at a time
    for (int year=end_year; year>start_year; year--) { 

      //Add this years benefit to y
      bj_ii(year, y); 
      
      // Integrate over [year,year-1]
      local_start_year = (int)year-1;
      local_end_year = (int)year;
      local_estimate();

      //Copy y to results
      Array.Copy(y, result[year-start_year-1], y.Length); 
    }

    return result;
  }

  /*************************** Auxiliaries ***************************/

  static float FindDoubleEpsilon() {
    float r = 1.0f;
    while (1.0f < (1.0f + r))
      r = r / 2.0f;
    return 2.0f * r;
  }

}

public class Timer {
  private Stopwatch stopwatch;

  public Timer() {
    stopwatch = new Stopwatch();
    stopwatch.Reset();
    stopwatch.Start();
  }

  public float Check() {
    return stopwatch.ElapsedMilliseconds;
  }
}

/******************************************************************************/

// The weird encapsulation style (classes PureEndowment and so on) is
// motivated by the wish to stay relatively close to C, so no
// inheritance.  Also, while the function argument dV could be
// implemented in C using a function pointer, this is probably
// impossible or a bad idea on a GPU. 

class CalculationSpecifications {

  public static void Main(String[] args) {
    //TimeAll(82288);
    //TestAll();
    PrintAll();
    //DeferredTemporaryLifeAnnuity.Print_dy();
    //PureEndowment.Print_dy();
  }

  public static readonly float err = 1e-7f;

  static readonly float age = 30.0f,
                  interestrate = 0.05f,
                  bpension = 1.0f,
                  pensiontime = 35.0f;

  static float indicator(bool b) {
    return b ? 1.0f : 0.0f;
  }

  static bool IsEqual(float[][] a,float [][] b) {
    for (int i=0; i<a.Length; i++) {
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

  public static void TimeAll(int customers) {
    //Console.WriteLine("PureEndowment:               " + PureEndowment.Time(customers));
    //Console.WriteLine("DeferredTemoraryLifeAnnuity: " + DeferredTemporaryLifeAnnuity.Time(customers));
    //Console.WriteLine("TemporaryLifeAnnuityPremium: " + TemporaryLifeAnnuityPremium.Time(customers));
    Console.WriteLine("TermInsurance:               " + TermInsurance.Time(customers));
    //Console.WriteLine("DisabilityAnnuity:           " + DisabilityAnnuity.Time(customers));
    //Console.WriteLine("DisabilityTermInsurance:     " + DisabilityTermInsurance.Time(customers));
  }


  public static void PrintAll() {
    //Print(PureEndowment.Compute());
    Print(DeferredTemporaryLifeAnnuity.Compute());
    /*Print(TemporaryLifeAnnuityPremium.Compute());
    Print(TermInsurance.Compute());
    Print(DisabilityAnnuity.Compute());
    Print(DisabilityTermInsurance.Compute());*/
  }

  public static void TestAll() {
    Assert(IsEqual(PureEndowment.Compute(),PureEndowment.test_values),"PureEndowment failed");
    Assert(IsEqual(DeferredTemporaryLifeAnnuity.Compute(),DeferredTemporaryLifeAnnuity.test_values),"DeferredTemporaryLifeAnnuity failed");
    Assert(IsEqual(TemporaryLifeAnnuityPremium.Compute(),TemporaryLifeAnnuityPremium.test_values),"TempLifeAnnuPrem failed");
    Assert(IsEqual(TermInsurance.Compute(),TermInsurance.test_values),"TempInsurance failed");
    Assert(IsEqual(DisabilityAnnuity.Compute(),DisabilityAnnuity.test_values),"DisAnnu failed");
    Assert(IsEqual(DisabilityTermInsurance.Compute(),DisabilityTermInsurance.test_values),"DisabilityTermInsurance failed");
    Console.WriteLine("Tests passed");
  }

  // Gompertz-Makeham mortality intensities for Danish women
  static float GM(float t) {
    //float r = 0.0005f + (float) Math.Round(Math.Pow(10.0f, 5.728f - 10.0f + 0.038f*(age + t));
    //Console.WriteLine("t: " + t + "GM(t): " + Math.Round(r,7)); 
    //return 0.5f;
    
    return 0.0005f + (float) Math.Round(Math.Pow(10.0f, 5.728f - 10.0f + 0.038f*(age + t)),7);

  }

  static float r(float t) { 
    return interestrate;    // Fixed interest rate
  }

  // Payment stream in state 0 (alive)
  static float b_0(float t) {
    return 0.0f;
  }

  // Lump sum payments while in state 0 (alive) 
  static float bj_00(float t) {
    // This works only because t is known to be an integer
    return t == pensiontime ? bpension : 0.0f;
  }

  // Lump sum payments while in state 1 (dead) 
  static float bj_11(float t) {
    // This works only because t is known to be an integer
    return 0.0f;
  }

  // Transition intensity from state 0 (alive) to state 1 (dead)
  static float mu_01(float t) {
    return GM(t);
  }

  // Lump sum payment on transition from state 0 (alive) to state 1 (dead)
  static float bj_01(float t) {
    return 0.0f;
  }

  //Print
  static void Print(float[][] result) {
    for (int y=0; y<result.Length; y++) {
      Console.Write("{0,3}:", y);
      for (int i=0; i<result[y].Length; i++)
        Console.Write("  {0,20:F7}", result[y][i]);
      Console.WriteLine();
    }
  }

  // The two-state Actulus calculation kernel examples; 
  // really one-state because V1(t) = 0 for all t.

  public class PureEndowment {

    static public float[][] test_values = new float[][] {
      new float[] {0.1437946974886250f},
          new float[] {0.1513594875720590f},
          new float[] {0.1593334830357050f},
          new float[] {0.1677404776222450f},
          new float[] {0.1766058889630630f},
          new float[] {0.1859569023966960f},
          new float[] {0.1958226315281160f},
          new float[] {0.2062342978884670f},
          new float[] {0.2172254324285080f},
          new float[] {0.2288321020171980f},
          new float[] {0.2410931646317720f},
          new float[] {0.2540505575320670f},
          new float[] {0.2677496234274150f},
          new float[] {0.2822394804908960f},
          new float[] {0.2975734430792770f},
          new float[] {0.3138095012097790f},
          new float[] {0.3310108682663330f},
          new float[] {0.3492466081065300f},
          new float[] {0.3685923547759590f},
          new float[] {0.3891311404830350f},
          new float[] {0.4109543504366150f},
          new float[] {0.4341628267168120f},
          new float[] {0.4588681476758340f},
          new float[] {0.4851941146396860f},
          new float[] {0.5132784841210240f},
          new float[] {0.5432749916555080f},
          new float[] {0.5753557231029890f},
          new float[] {0.6097139012822690f},
          new float[] {0.6465671707381610f},
          new float[] {0.6861614820516400f},
          new float[] {0.7287757004081980f},
          new float[] {0.7747270924521970f},
          new float[] {0.8243778824987340f},
          new float[] {0.8781431162166330f},
          new float[] {0.9365001299367070f},
          new float[] {0.0000000000000000f},
          new float[] {0.0000000000000000f},
          new float[] {0.0000000000000000f},
          new float[] {0.0000000000000000f},
          new float[] {0.0000000000000000f},
          new float[] {0.0000000000000000f}
    };

    static float b_0(float t) {
      return 0.0f;
    }

    static float mu_01(float t) {
      return GM(t);
    }

    static float bj_00(float t) {
      // This works only because t is known to be an integer
      return t == pensiontime ? bpension : 0.0f;
    }

    static float bj_01(float t) {
      return 0.0f;
    }

    public static void Print_dy() {
      Estimator estimator = new Estimator();
      estimator.dy =
        (float t, float[] V, float[] res) =>
        res[0] = r(t) * V[0] - b_0(t) - mu_01(t) * (0 - V[0] + bj_01(t));
      
      float[] result = new float[1];
      float[] W      = new float[1];
      W[0]            = 0;
      estimator.dy(0,W,result);
      Console.WriteLine(result[0]);

      W[0]            = 1;
      estimator.dy(0,W,result);
      Console.WriteLine(result[0]);

      W[0]            = 1;
      estimator.dy(1,W,result);
      Console.WriteLine(result[0]);

      W[0]            = -2;
      estimator.dy(1,W,result);
      Console.WriteLine(result[0]);
    }

    public static float Time(int customers) {
      Timer timer = new Timer();
      float start_time = timer.Check();
      for(int i = 0;i<customers;i++)
        Compute();
      float end_time = timer.Check();
      return end_time - start_time;
    }

    public static float[][] Compute() {
      //Construct the estimator
      Estimator estimator = new Estimator();
      estimator.construct(1);

      //Set estimator variables (Term insurrance)
      estimator.relerr = err;
      estimator.abserr = err;
      estimator.end_year_y = new float[] { 0 };
      estimator.start_year = 0;
      estimator.end_year = 40;
      estimator.bj_ii =
        (float t, float[] res) =>
        res[0] += bj_00(t);
      estimator.dy =
        (float t, float[] V, float[] res) =>
        res[0] = r(t) * V[0] - b_0(t) - mu_01(t) * (0 - V[0] + bj_01(t));

      return estimator.estimate();
    }
  }

  public class DeferredTemporaryLifeAnnuity {

    static public float[][] test_values = new float[][] {
      new float[] {    1.0265607676014400f},
          new float[] {    1.0805663523022000f},
          new float[] {    1.1374932838714300f},
          new float[] {    1.1975114275626800f},
          new float[] {    1.2608022416886800f},
          new float[] {    1.3275598043516300f},
          new float[] {    1.3979919596875600f},
          new float[] {    1.4723216004702400f},
          new float[] {    1.5507881065874600f},
          new float[] {    1.6336489620308100f},
          new float[] {    1.7211815767162100f},
          new float[] {    1.8136853437812100f},
          new float[] {    1.9114839681141600f},
          new float[] {    2.0149281079127500f},
          new float[] {    2.1243983782341300f},
          new float[] {    2.2403087740143500f},
          new float[] {    2.3631105801829600f},
          new float[] {    2.4932968486250200f},
          new float[] {    2.6314075362739100f},
          new float[] {    2.7780354160836500f},
          new float[] {    2.9338328936850100f},
          new float[] {    3.0995198879959800f},
          new float[] {    3.2758929649601700f},
          new float[] {    3.4638359512170500f},
          new float[] {    3.6643323004950200f},
          new float[] {    3.8784795419259000f},
          new float[] {    4.1075062089359300f},
          new float[] {    4.3527917332333000f},
          new float[] {    4.6158898949982800f},
          new float[] {    4.8985565532549200f},
          new float[] {    5.2027825467745500f},
          new float[] {    5.5308328651267400f},
          new float[] {    5.8852934539509200f},
          new float[] {    6.2691273543593000f},
          new float[] {    6.6857423050148000f},
          new float[] {    7.1390724853580000f},
          new float[] {    6.5992994744162900f},
          new float[] {    6.0319762072565800f},
          new float[] {    5.4341959778729600f},
          new float[] {    4.8025044157278400f},
          new float[] {    4.1327786917263300f},
          new float[] {    3.4200771644266400f},
          new float[] {    2.6584512106845400f},
          new float[] {    1.8407083534208200f},
          new float[] {    0.9581122236639260f},
          new float[] {    0.0000000000000000f},
          new float[] {    0.0000000000000000f},
          new float[] {    0.0000000000000000f},
          new float[] {    0.0000000000000000f},
          new float[] {    0.0000000000000000f},
          new float[] {    0.0000000000000000f}
    };

    static int m = 35, n = 10;

    static float b_0(float t) {
      return bpension * indicator(t > m) * indicator(t < m + n);
    }

    static float mu_01(float t) {
      return GM(t);
    }

    static float bj_00(float t) {
      return 0.0f;
    }

    static float bj_01(float t) {
      return 0.0f;
    }

    public static float Time(int customers) {
      Timer timer = new Timer();
      float start_time = timer.Check();
      for(int i = 0;i<customers;i++)
        Compute();
      float end_time = timer.Check();
      return end_time - start_time;
    }
    public static void Print_dy() {
      Estimator estimator = new Estimator();
      estimator.dy =
        (float t, float[] V, float[] res) =>
        res[0] = r(t) * V[0] - b_0(t) - mu_01(t) * (0 - V[0] + bj_01(t));
      
      float[] result = new float[1];
      float[] W      = new float[1];

      W[0]            = 0;
      estimator.dy(0,W,result);
      Console.WriteLine(result[0]);

      W[0]            = 1;
      estimator.dy(0,W,result);
      Console.WriteLine(result[0]);

      W[0]            = 1;
      estimator.dy(1,W,result);
      Console.WriteLine(result[0]);

      W[0]            = -2;
      estimator.dy(1,W,result);
      Console.WriteLine(result[0]);
    }

    public static float[][] Compute() {
      Estimator estimator = new Estimator();
      estimator.construct(1);
      estimator.relerr = err;
      estimator.abserr = err;
      estimator.start_year = 0;
      estimator.end_year = 50;
      estimator.end_year_y = new float[] { 0 };
      estimator.bj_ii =
        (float t, float[] res) => res[0] += bj_00(t);
      estimator.dy =
        (float t, float[] V, float[] res) =>
        res[0] = r(t) * V[0] - b_0(t) - mu_01(t) * (0 - V[0] + bj_01(t));

      return estimator.estimate();
    }
  }

  public class TemporaryLifeAnnuityPremium {

    static public float[][] test_values = new float[][] {
      new float[] {  -15.9717676660001000f},
          new float[] {  -15.7859295725898000f},
          new float[] {  -15.5914495774420000f},
          new float[] {  -15.3879467041606000f},
          new float[] {  -15.1750230434772000f},
          new float[] {  -14.9522626687192000f},
          new float[] {  -14.7192304105538000f},
          new float[] {  -14.4754704664848000f},
          new float[] {  -14.2205048162637000f},
          new float[] {  -13.9538314093213000f},
          new float[] {  -13.6749220843743000f},
          new float[] {  -13.3832201743607000f},
          new float[] {  -13.0781377416068000f},
          new float[] {  -12.7590523783886000f},
          new float[] {  -12.4253034965330000f},
          new float[] {  -12.0761880160465000f},
          new float[] {  -11.7109553465719000f},
          new float[] {  -11.3288015361432000f},
          new float[] {  -10.9288624386922000f},
          new float[] {  -10.5102057241661000f},
          new float[] {  -10.0718215219905000f},
          new float[] {   -9.6126114486954800f},
          new float[] {   -9.1313757222581100f},
          new float[] {   -8.6267980071471600f},
          new float[] {   -8.0974275627022600f},
          new float[] {   -7.5416581802140900f},
          new float[] {   -6.9577032868934500f},
          new float[] {   -6.3435664627206500f},
          new float[] {   -5.6970064523888100f},
          new float[] {   -5.0154955507238000f},
          new float[] {   -4.2961699850952200f},
          new float[] {   -3.5357705981019300f},
          new float[] {   -2.7305717294876500f},
          new float[] {   -1.8762956830915000f},
          new float[] {   -0.9680095100191570f},
          new float[] {    0.0000000000000000f},
          new float[] {    0.0000000000000000f},
          new float[] {    0.0000000000000000f},
          new float[] {    0.0000000000000000f},
          new float[] {    0.0000000000000000f},
          new float[] {    0.0000000000000000f},
          new float[] {    0.0000000000000000f},
          new float[] {    0.0000000000000000f},
          new float[] {    0.0000000000000000f},
          new float[] {    0.0000000000000000f},
          new float[] {    0.0000000000000000f},
          new float[] {    0.0000000000000000f},
          new float[] {    0.0000000000000000f},
          new float[] {    0.0000000000000000f},
          new float[] {    0.0000000000000000f},
          new float[] {    0.0000000000000000f}
    };

    static int n = 35;
    static float bpremium = 1;

    static float b_0(float t) {
      return -bpremium * indicator(t >= 0) * indicator(t < n);
    }

    static float mu_01(float t) {
      return GM(t);
    }

    static float bj_00(float t) {
      return 0.0f;
    }

    static float bj_01(float t) {
      return 0.0f;
    }

    public static float Time(int customers) {
      Timer timer = new Timer();
      float start_time = timer.Check();
      for(int i = 0;i<customers;i++)
        Compute();
      float end_time = timer.Check();
      return end_time - start_time;
    }

    public static float[][] Compute() {
      Estimator estimator = new Estimator();
      estimator.construct(1);
      estimator.relerr = err;
      estimator.abserr = err;
      estimator.start_year = 0;
      estimator.end_year = 50;
      estimator.end_year_y = new float[] { 0 };
      estimator.bj_ii =
        (float t, float[] res) => res[0] += bj_00(t);
      estimator.dy =
        (float t, float[] V, float[] res) =>
        res[0] = r(t) * V[0] - b_0(t) - mu_01(t) * (0 - V[0] + bj_01(t));

      return estimator.estimate();
    }
  }

  public class TermInsurance {


    static public float[][] test_values = new float[][] {
      new float[] {    0.0576169193132673f},
          new float[] {    0.0593440338503396f},
          new float[] {    0.0610940381466785f},
          new float[] {    0.0628621872268886f},
          new float[] {    0.0646429589229975f},
          new float[] {    0.0664299642301075f},
          new float[] {    0.0682158480098718f},
          new float[] {    0.0699921788564634f},
          new float[] {    0.0717493268311646f},
          new float[] {    0.0734763275934875f},
          new float[] {    0.0751607312303756f},
          new float[] {    0.0767884338351089f},
          new float[] {    0.0783434895820515f},
          new float[] {    0.0798079006843388f},
          new float[] {    0.0811613821939618f},
          new float[] {    0.0823810980933246f},
          new float[] {    0.0834413645163828f},
          new float[] {    0.0843133152038705f},
          new float[] {    0.0849645234136384f},
          new float[] {    0.0853585734399381f},
          new float[] {    0.0854545736024947f},
          new float[] {    0.0852066009948634f},
          new float[] {    0.0845630663660149f},
          new float[] {    0.0834659851665421f},
          new float[] {    0.0818501379168519f},
          new float[] {    0.0796420995167974f},
          new float[] {    0.0767591127460391f},
          new float[] {    0.0731077757869581f},
          new float[] {    0.0685825068600615f},
          new float[] {    0.0630637406431617f},
          new float[] {    0.0564158005824136f},
          new float[] {    0.0484843779035996f},
          new float[] {    0.0390935313045591f},
          new float[] {    0.0280420999246517f},
          new float[] {    0.0150993948779494f},
          new float[] {    0.0000000000000000f},
          new float[] {    0.0000000000000000f},
          new float[] {    0.0000000000000000f},
          new float[] {    0.0000000000000000f},
          new float[] {    0.0000000000000000f},
          new float[] {    0.0000000000000000f},
          new float[] {    0.0000000000000000f},
          new float[] {    0.0000000000000000f},
          new float[] {    0.0000000000000000f},
          new float[] {    0.0000000000000000f},
          new float[] {    0.0000000000000000f},
          new float[] {    0.0000000000000000f},
          new float[] {    0.0000000000000000f},
          new float[] {    0.0000000000000000f},
          new float[] {    0.0000000000000000f},
          new float[] {    0.0000000000000000f}
    };

    static int n = 35;
    static float bdeath = 1;

    static float b_0(float t) {
      return 0.0f;
    }

    static float mu_01(float t) {
      return GM(t);
    }

    static float bj_00(float t) {
      return 0.0f;
    }

    static float bj_01(float t) {
      return bdeath * indicator(t > 0) * indicator(t < n);
    }

    public static float Time(int customers) {
      Timer timer = new Timer();
      float start_time = timer.Check();
      for(int i = 0;i<customers;i++)
        Compute();
      float end_time = timer.Check();
      return end_time - start_time;
    }

    public static float[][] Compute() {
      Estimator estimator = new Estimator();
      estimator.construct(1);
      estimator.relerr = err;
      estimator.abserr = err;
      estimator.start_year = 0;
      estimator.end_year = 50;
      estimator.end_year_y = new float[] { 0 };
      estimator.bj_ii =
        (float t, float[] res) => res[0] += bj_00(t);
      estimator.dy =
        (float t, float[] V, float[] res) =>
        res[0] = r(t) * V[0] - b_0(t) - mu_01(t) * (0 - V[0] + bj_01(t));

      return estimator.estimate();
    }
  }

  // The three-state Actulus calculation kernel examples; 
  // really two-state because V2(t) = 0 for all t.

  public class DisabilityAnnuity {

    static public float[][] test_values = new float[][] {
          new float[] {    0.5555261079604120f,   15.9717676673750000f},
          new float[] {    0.5697939362470290f,   15.7859295725873000f},
          new float[] {    0.5842458860490700f,   15.5914495774393000f},
          new float[] {    0.5988112559939020f,   15.3879467041578000f},
          new float[] {    0.6134061668653150f,   15.1750230434743000f},
          new float[] {    0.6279320187147260f,   14.9522626687161000f},
          new float[] {    0.6422738476905750f,   14.7192304105506000f},
          new float[] {    0.6562985942478250f,   14.4754704664814000f},
          new float[] {    0.6698533011953220f,   14.2205048162601000f},
          new float[] {    0.6827632688462940f,   13.9538314093175000f},
          new float[] {    0.6948302058887720f,   13.6749220843703000f},
          new float[] {    0.7058304291697920f,   13.3832201743564000f},
          new float[] {    0.7155131842695250f,   13.0781377416023000f},
          new float[] {    0.7235991826741960f,   12.7590523783855000f},
          new float[] {    0.7297794820426480f,   12.4253034965300000f},
          new float[] {    0.7337148754958810f,   12.0761880160438000f},
          new float[] {    0.7350360067043140f,   11.7109553465694000f},
          new float[] {    0.7333444934112320f,   11.3288015361409000f},
          new float[] {    0.7282154278216300f,   10.9288624386901000f},
          new float[] {    0.7192017347676550f,   10.5102057241642000f},
          new float[] {    0.7058410171284570f,   10.0718215219888000f},
          new float[] {    0.6876657158018430f,    9.6126114486939200f},
          new float[] {    0.6642176772198330f,    9.1313757222565700f},
          new float[] {    0.6350685815395070f,    8.6267980071455300f},
          new float[] {    0.5998481774680660f,    8.0974275627005400f},
          new float[] {    0.5582829507417500f,    7.5416581802122700f},
          new float[] {    0.5102488039940220f,    6.9577032868915200f},
          new float[] {    0.4558426666237880f,    6.3435664627186100f},
          new float[] {    0.3954798644641650f,    5.6970064523866500f},
          new float[] {    0.3300268327704200f,    5.0154955507223200f},
          new float[] {    0.2609827682846260f,    4.2961699850940100f},
          new float[] {    0.1907297303473690f,    3.5357705981010400f},
          new float[] {    0.1228795253423570f,    2.7305717294870500f},
          new float[] {    0.0627590447157113f,    1.8762956830912200f},
          new float[] {    0.0180961562710709f,    0.9680095100190450f},
          new float[] {    0.0000000000000000f,    0.0000000000000000f},
          new float[] {    0.0000000000000000f,    0.0000000000000000f},
          new float[] {    0.0000000000000000f,    0.0000000000000000f},
          new float[] {    0.0000000000000000f,    0.0000000000000000f},
          new float[] {    0.0000000000000000f,    0.0000000000000000f},
          new float[] {    0.0000000000000000f,    0.0000000000000000f},
          new float[] {    0.0000000000000000f,    0.0000000000000000f},
          new float[] {    0.0000000000000000f,    0.0000000000000000f},
          new float[] {    0.0000000000000000f,    0.0000000000000000f},
          new float[] {    0.0000000000000000f,    0.0000000000000000f},
          new float[] {    0.0000000000000000f,    0.0000000000000000f},
          new float[] {    0.0000000000000000f,    0.0000000000000000f},
          new float[] {    0.0000000000000000f,    0.0000000000000000f},
          new float[] {    0.0000000000000000f,    0.0000000000000000f},
          new float[] {    0.0000000000000000f,    0.0000000000000000f},
          new float[] {    0.0000000000000000f,    0.0000000000000000f}
    };

    static int n = 35;
    static float bdisabled = 1;

    static float b_0(float t) {
      return 0.0f;
    }

    static float b_1(float t) {
      return bdisabled * indicator(t > 0) * indicator(t < n);
    }

    static float GM01(float t) {
      return 0.0006f + (float) Math.Round(Math.Pow(10, 4.71609f - 10 + 0.06f*(age + t)),7);
    }

    static float GM02(float t) {
      return GM(t);
    }

    static float GM12(float t) {
      return GM(t);
    }

    static float mu_01(float t) {
      return GM01(t);
    }

    static float mu_02(float t) {
      return GM02(t);
    }

    static float mu_12(float t) {
      return GM12(t);
    }

    static float bj_00(float t) {
      return 0.0f;
    }

    static float bj_01(float t) {
      return 0.0f;
    }

    static float bj_02(float t) {
      return 0.0f;
    }

    static float bj_11(float t) {
      return 0.0f;
    }

    static float bj_12(float t) {
      return 0.0f;
    }

    public static void Print_dy() {
      Estimator estimator = new Estimator();
      estimator.dy =
        (float t, float[] V, float[] res) => {
          res[0] = r(t) * V[0] - b_0(t) - mu_01(t) * (V[1] - V[0] + bj_01(t)) - mu_02(t) * (0 - V[0] + bj_02(t));
          res[1] = r(t) * V[1] - b_1(t) - mu_12(t) * (0 - V[1] + bj_12(t)); 
        };
      
      float[] result = new float[2];
      float[] W      = new float[2];
      W[0]            = 0;
      estimator.dy(0,W,result);
      Console.WriteLine(result[0]);

      W[0]            = 1;
      estimator.dy(0,W,result);
      Console.WriteLine(result[0]);

      W[0]            = 1;
      estimator.dy(1,W,result);
      Console.WriteLine(result[0]);

      W[0]            = -2;
      estimator.dy(1,W,result);
      Console.WriteLine(result[0]);

      W[0]            = 0;
      estimator.dy(34.625f,W,result);
      Console.WriteLine("result: " + result[1]);
    }

    public static float Time(int customers) {
      Timer timer = new Timer();
      float start_time = timer.Check();
      for(int i = 0;i<customers;i++)
        Compute();
      float end_time = timer.Check();
      return end_time - start_time;
    }

    public static float[][] Compute() {
      Estimator estimator = new Estimator();
      estimator.construct(2);
      estimator.relerr = err;
      estimator.abserr = err;
      estimator.start_year = 0;
      estimator.end_year = 50;
      estimator.end_year_y = new float[] { 0,0 };
      estimator.bj_ii =
        (float t, float[] res) => {res[0] += bj_00(t); res[1] += bj_11(t);};
      estimator.dy =
        (float t, float[] V, float[] res) => {
          res[0] = r(t) * V[0] - b_0(t) - mu_01(t) * (V[1] - V[0] + bj_01(t)) - mu_02(t) * (0 - V[0] + bj_02(t));
          res[1] = r(t) * V[1] - b_1(t) - mu_12(t) * (0 - V[1] + bj_12(t)); 
        };

      return estimator.estimate();
    }
  }

  public class DisabilityTermInsurance {

    static public float[][] test_values = new float[][] {
      new float[] {    0.0714186989824431f,    0.0000000000000000f},
          new float[] {    0.0742705084048387f,    0.0000000000000000f},
          new float[] {    0.0772312485352858f,    0.0000000000000000f},
          new float[] {    0.0803007388769112f,    0.0000000000000000f},
          new float[] {    0.0834779380618836f,    0.0000000000000000f},
          new float[] {    0.0867607747788289f,    0.0000000000000000f},
          new float[] {    0.0901459511899878f,    0.0000000000000000f},
          new float[] {    0.0936287143026796f,    0.0000000000000000f},
          new float[] {    0.0972025899462083f,    0.0000000000000000f},
          new float[] {    0.1008590729874850f,    0.0000000000000000f},
          new float[] {    0.1045872661636680f,    0.0000000000000000f},
          new float[] {    0.1083734583544920f,    0.0000000000000000f},
          new float[] {    0.1122006311731980f,    0.0000000000000000f},
          new float[] {    0.1160478803013540f,    0.0000000000000000f},
          new float[] {    0.1198897349010070f,    0.0000000000000000f},
          new float[] {    0.1236953544561090f,    0.0000000000000000f},
          new float[] {    0.1274275773038030f,    0.0000000000000000f},
          new float[] {    0.1310417884996260f,    0.0000000000000000f},
          new float[] {    0.1344845660310780f,    0.0000000000000000f},
          new float[] {    0.1376920530492920f,    0.0000000000000000f},
          new float[] {    0.1405879887686990f,    0.0000000000000000f},
          new float[] {    0.1430813106541400f,    0.0000000000000000f},
          new float[] {    0.1450632136041590f,    0.0000000000000000f},
          new float[] {    0.1464035154086170f,    0.0000000000000000f},
          new float[] {    0.1469461280481600f,    0.0000000000000000f},
          new float[] {    0.1465033660147020f,    0.0000000000000000f},
          new float[] {    0.1448487279175560f,    0.0000000000000000f},
          new float[] {    0.1417076547118570f,    0.0000000000000000f},
          new float[] {    0.1367455798907500f,    0.0000000000000000f},
          new float[] {    0.1295523183539230f,    0.0000000000000000f},
          new float[] {    0.1196214525705190f,    0.0000000000000000f},
          new float[] {    0.1063228073533260f,    0.0000000000000000f},
          new float[] {    0.0888652648692060f,    0.0000000000000000f},
          new float[] {    0.0662459119488600f,    0.0000000000000000f},
          new float[] {    0.0371795952968856f,    0.0000000000000000f},
          new float[] {    0.0000000000000000f,    0.0000000000000000f},
          new float[] {    0.0000000000000000f,    0.0000000000000000f},
          new float[] {    0.0000000000000000f,    0.0000000000000000f},
          new float[] {    0.0000000000000000f,    0.0000000000000000f},
          new float[] {    0.0000000000000000f,    0.0000000000000000f},
          new float[] {    0.0000000000000000f,    0.0000000000000000f},
          new float[] {    0.0000000000000000f,    0.0000000000000000f},
          new float[] {    0.0000000000000000f,    0.0000000000000000f},
          new float[] {    0.0000000000000000f,    0.0000000000000000f},
          new float[] {    0.0000000000000000f,    0.0000000000000000f},
          new float[] {    0.0000000000000000f,    0.0000000000000000f},
          new float[] {    0.0000000000000000f,    0.0000000000000000f},
          new float[] {    0.0000000000000000f,    0.0000000000000000f},
          new float[] {    0.0000000000000000f,    0.0000000000000000f},
          new float[] {    0.0000000000000000f,    0.0000000000000000f},
          new float[] {    0.0000000000000000f,    0.0000000000000000f}
    };

    static int n = 35;
    static float bdisabled = 1;

    static float b_0(float t) {
      return 0.0f;
    }

    static float b_1(float t) {
      return 0.0f;
    }

    static float GM01(float t) {
      return 0.0006f + (float) Math.Round(Math.Pow(10, 4.71609f - 10 + 0.06f*(age + t)),7);
    }

    static float GM02(float t) {
      return GM(t);
    }

    static float GM12(float t) {
      return GM(t);
    }

    static float mu_01(float t) {
      return GM01(t);
    }

    static float mu_02(float t) {
      return GM02(t);
    }

    static float mu_12(float t) {
      return GM12(t);
    }

    static float bj_00(float t) {
      return 0.0f;
    }

    static float bj_01(float t) {
      return bdisabled * indicator(t > 0) * indicator(t < n);
    }

    static float bj_02(float t) {
      return 0.0f;
    }

    static float bj_11(float t) {
      return 0.0f;
    }    

    static float bj_12(float t) {
      return 0.0f;
    }

    public static float Time(int customers) {
      Timer timer = new Timer();
      float start_time = timer.Check();
      for(int i = 0;i<customers;i++)
        Compute();
      float end_time = timer.Check();
      return end_time - start_time;
    }

    public static float[][] Compute() {
      Estimator estimator = new Estimator();
      estimator.construct(2);
      estimator.relerr = err;
      estimator.abserr = err;
      estimator.start_year = 0;
      estimator.end_year = 50;
      estimator.end_year_y = new float[] { 0,0 };
      estimator.bj_ii =
        (float t, float[] res) => { res[0] += bj_00(t); res[1] += bj_11(t); };
      estimator.dy =
        (float t, float[] V, float[] res) => {
          res[0] = r(t) * V[0] - b_0(t) - mu_01(t) * (V[1] - V[0] + bj_01(t)) - mu_02(t) * (0 - V[0] + bj_02(t));
          res[1] = r(t) * V[1] - b_1(t) - mu_12(t) * (0 - V[1] + bj_12(t));
        };

      return estimator.estimate();
    }
  }
}
