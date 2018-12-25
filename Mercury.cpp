/*
  Mercury: A 2nd-order Finite Volume code for
  the Euler equations of inviscid, compressible
  flow of an ideal gas in 2D. 

  Boundary Conditions: Code's configuration
  is zero-penetration walls on top and 
  bottom. 
  Right and Left edges are programmed
  as farfield outflows.
*/

#include<stdio.h>
#include<math.h>
#include<vector>
#include<iostream>

//resolution
const int Mx = 256;
const int My = Mx/2;

//Geometry
double LEFT = 0.0;
double RIGHT = 2.0;
double BASE = 0.0;
double TOP = 1.0;

//screen update rate
const int interval_out = 100;

//quadrature for initialization
const int GQdeg = 6;

//fluid parameters
const double rsh = 1.4; //ratio of specific heats
const double Rg = 287.15; //gas constant
const double cv = Rg/(rsh-1);
const double cp = rsh*Rg / (rsh-1.0);

//Some constants (These are for the exact Riemann solver) 
const double G1 = (rsh-1.0)/(2.0*rsh);
const double G2 = (rsh + 1.0) / (2.0*rsh);
const double G3 = 2.0*rsh/(rsh-1.0);
const double G4 = 2.0 / (rsh-1.0);
const double G5 = 2.0 / (rsh+1.0);
const double G6 = (rsh-1.0) / (rsh+1.0);
const double G7 = (rsh-1.0) / 2.0;
const double G8 = rsh-1.0;
const double TOL_NR = pow(10,-10.0);
const double Chi = (rsh-1.0)/(rsh+1.0);

//Steady atmospheric parameters, can be used for initial condition (IC):
const double a_free = 1.0;
const double Mach_free = 0.5;
const double p_free = 100000.0;
const double T_free = 300.0;
const double rho_free = p_free / (Rg*T_free);
const double u_free = Mach_free * sqrt( rsh*Rg*T_free );
const double v_free = 0;

const double shock_loc = 0.2;

const double V_shock = 1.1 * a_free;
const int initialize = 0; //initialization routine: standard if 0, continue if 1
const double t_start = 0.0;

//Rault shock-vortex parameters
const double t_final = 15.0/30.0;
const double Mshock = 1.5;
const double Mvort = 0.5;
const double rad_a = 0.075;
const double rad_b = 0.175;
const double rhoRault_0 = 1.0; //pre-shock density
const double TRault_0 = 1.0; //pre-shock temperature
const double chiRault = (rsh-1.0) / (Rg*rsh);
const double epiX_Rault = 0.25;
const double epiY_Rault = 0.5;

const double eps = pow(10,-15.0); //limiter function safety factor

//Math constants, for general use:
const double pi = 3.141592653589793238;
const double e = 2.71828182845904523536;

//CFL number: in the case of an unstable calculation, lowering this number can help.
//Setting it higher than 1.0 is inadviseable.
const double nu = 0.5;

//The cells are rectangular; delbox_x and delbox_y are the
//cell dimensions.
const double delbox_x = (RIGHT-LEFT) / (Mx+0.0);
const double delbox_y = (TOP-BASE) / (My+0.0);


double weights[GQdeg];

double t_check[11]; //Array of the times where I pause code to OUTPUT SOLUTION TO FILE
int t_mark = 0;
double t_target;

typedef struct
{
  double ma;
  double mx;
  double my;
  double en;
} vec4;
typedef struct
{
  double loc;
  double weight;
} gwn;
class timer 
{
  //Simple timer class
 public:
  clock_t start;
  std::string name;
  double elapsed;
  
  timer(std::string n) : name(n), elapsed(0.0f) {}
  ~timer(){}
  
  void startBlock() 
  {
    start = clock();
  }
  
  void endBlock() 
  {
    elapsed += (double)(clock() - start)/CLOCKS_PER_SEC;
  }
  
  void print() {
    printf("%s total time: %fs\n", name.c_str(), elapsed);
  }
};
typedef std::vector<std::vector<double> > D2_double;
typedef std::vector<std::vector<vec4> > D2_vec4;
typedef std::vector<std::vector<std::vector<std::vector<double> > > > D4_double;
typedef std::vector<std::vector<std::vector<std::vector<vec4> > > > D4_vec4;

std::vector<timer*> timers;
/*timer Roe_Timer0L("Roe section 0Left");
timer Roe_Timer0R("Roe section 0Right");
timer Roe_Timer1("Roe section 1");
timer Roe_Timer2("Roe section 2");
timer Roe_Timer3("Roe section 3");
timer Roe_Timer4("Roe section 4");
timer Roe_Timer5("Roe section 5");

timer Flux_E_Timer("E flux function");
timer Flux_F_Timer("F flux function");*/
//Output array
D2_vec4 U_out;
D2_double x;
D2_double y;
//std::vector<std::vector<vec4> > U_out;
//std::vector<std::vector<double> > x;
//std::vector<std::vector<double> > y;

const vec4 eps_vec = {eps,eps,eps,eps};

double get_pressure(vec4 input)
{
  //Function to extract pressure from U vector
  double rho = input.ma;
  double u = input.mx / rho;
  double v = input.my / rho;
  double sq_Vmag = u*u + v*v;
  double p = (rsh-1.0) * (input.en - 0.5*rho*sq_Vmag);
  return p;
}
void write_output(int station)
{
  //Function for writing solution to file
  //when needed.
  double rho;
  double u;
  double v;
  double Eg;
  double p;
  double Vel;
  FILE*fileA;
  switch (station)
    {
    case 0:
      {
	fileA = fopen("xU_t0.csv","w");
	break;
      }	
    case 1:
      {
	fileA = fopen("xU_t1.csv","w");
	break;
      }
    case 2:
      {
	fileA = fopen("xU_t2.csv","w");
	break;
      }
    case 3:
      {
	fileA = fopen("xU_t3.csv","w");
	break;
      }
    case 4:
      {
	fileA = fopen("xU_t4.csv","w");
	break;
      }
    case 5:
      {
	fileA = fopen("xU_t5.csv","w");
	break;
      }
    case 6:
      {
	fileA = fopen("xU_t6.csv","w");
	break;
      }
    case 7:
      {
	fileA = fopen("xU_t7.csv","w");
	break;
      }
    case 8:
      {
	fileA = fopen("xU_t8.csv","w");
	break;
      }
    case 9:
      {
	fileA = fopen("xU_t9.csv","w");
	break;
      }
    case 10:
      {
	fileA = fopen("xU_t10.csv","w");
	break;
      }
    default:
      {
	break;
      }
    }
  //Loop through all the cells, output each cell's current solution.
  for (int I = 0; I < Mx; I = I + 1)
    {
      for (int J = 0; J < My; J = J + 1)
	{
	  rho = U_out[I][J].ma;
	  u = U_out[I][J].mx/rho;
	  v = U_out[I][J].my/rho;
	  Eg = U_out[I][J].en/rho;
	  p = get_pressure(U_out[I][J]);
	  Vel = sqrt(u*u + v*v);
	  
	  fprintf(fileA,"%8.7f, %8.7f, %16.10f,%16.10f,%16.10f,%16.10f,%16.10f,%16.10f\n", x[I][J], y[I][J], rho, u, v, Eg, p, Vel);
	}
    }
  fclose(fileA);
}
D2_double Allocate_2D_double(int d1, int d2)
{
  std::vector<std::vector<double> > output;
  output.resize(d1);
  for (int A = 0; A < d1; A = A+1)
    {
      output[A].resize(d2);
    }
  return output;
}
D2_vec4 Allocate_2D_vec4(int d1, int d2) 
{
  std::vector<std::vector<vec4> > output;
  output.resize(d1);
  for (int A = 0; A < d1; A = A+1)
    {
      output[A].resize(d2);
    }
  return output;
}
D4_double Allocate_4D_double(int d1, int d2, int d3, int d4)
{
  std::vector<std::vector<std::vector<std::vector<double> > > > output;
  output.resize(d1);
  for (int A = 0; A < d1; A = A + 1)
    {
      output[A].resize(d2);
      for (int B = 0; B < d2; B = B + 1)
	{
	  output[A][B].resize(d3);
	  for (int C = 0; C < d3; C = C + 1)
	    {
	      output[A][B][C].resize(d4);
	    }
	}
    }
  return output;
}
D4_vec4 Allocate_4D_vec4(int d1, int d2, int d3, int d4)
{
  std::vector<std::vector<std::vector<std::vector<vec4> > > > output;
  output.resize(d1);
  for (int A = 0; A < d1; A = A + 1)
    {
      output[A].resize(d2);
      for (int B = 0; B < d2; B = B + 1)
	{
	  output[A][B].resize(d3);
	  for (int C = 0; C < d3; C = C + 1)
	    {
	      output[A][B][C].resize(d4);
	    }
	}
    }
  return output;
}
gwn gwn_return(int dim, int spot)
{
  //Returns Gaussian weights and nodes for a particular quadrature precision.:
  gwn output;
 
  double weight20[20];	//nodal weights
  double GQ20[20];	//nodal points
  
  double weight10[10];
  double GQ10[10];
  
  double weight6[6];
  double GQ6[6];

  double weight5[5];
  double GQ5[5];
  
  double weight4[4];
  double GQ4[4];
  
  double weight3[3];
  double GQ_loc_3[3];
  
  double weight2[2];
  double GQ_loc_2[2];

  
  weight20[9] = 0.1527533871307258506980843;
  weight20[8] = 0.1491729864726037467878287;
  weight20[7] = 0.1420961093183820513292983;
  weight20[6] = 0.1316886384491766268984945;
  weight20[5] = 0.1181945319615184173123774;
  weight20[4] = 0.1019301198172404350367501;
  weight20[3] = 0.0832767415767047487247581;
  weight20[2] = 0.0626720483341090635695065;
  weight20[1] = 0.0406014298003869413310400;
  weight20[0] = 0.0176140071391521183118620;
  GQ20[9] = -0.0765265211334973337546404;
  GQ20[8] = -0.2277858511416450780804962;
  GQ20[7] = -0.3737060887154195606725482;
  GQ20[6] = -0.5108670019508270980043641;
  GQ20[5] = -0.6360536807265150254528367;
  GQ20[4] = -0.7463319064601507926143051;
  GQ20[3] = -0.8391169718222188233945291;
  GQ20[2] = -0.9122344282513259058677524;
  GQ20[1] = -0.9639719272779137912676661;
  GQ20[0] = -0.9931285991850949247861224;
  
  weight10[4] = 0.2955242247147528701738930;
  weight10[3] = 0.2692667193099963550912269;
  weight10[2] = 0.2190863625159820439955349;
  weight10[1] = 0.1494513491505805931457763;
  weight10[0] = 0.0666713443086881375935688;
  GQ10[4] = -0.1488743389816312108848260;
  GQ10[3] = -0.4333953941292471907992659;
  GQ10[2] = -0.6794095682990244062343274;
  GQ10[1] = -0.8650633666889845107320967;
  GQ10[0] = -0.9739065285171717200779640;
  
  weight6[2] = 0.4679139345726910473898703;
  weight6[1] = 0.3607615730481386075698335;
  weight6[0] = 0.1713244923791703450402961;
  GQ6[2] = -0.2386191860831969086305017;
  GQ6[1] = -0.6612093864662645136613996;
  GQ6[0] = -0.9324695142031520278123016;

  weight4[1] = 0.6521451548625461;
  weight4[0] = 0.3478548451374538;
  GQ4[1] = -0.3399810435848563;
  GQ4[0] = -0.8611363115940526;

  GQ_loc_2[0] = -0.5773502691896257;
  weight2[0] = 1.0000000000000000;

  weight3[2] = 0.5555555555555556;
  weight3[1] = 0.8888888888888888;
  weight3[0] = 0.5555555555555556;
  GQ_loc_3[2] = 0.7745966692414834;
  GQ_loc_3[1] = 0.0000000000000000;
  GQ_loc_3[0] = -0.7745966692414834;
  
 

  GQ5[0] = -0.9061798459386640;
  GQ5[1] = -0.5384693101056831;
  GQ5[2] = 0.0000000000000000;
  GQ5[3] = 0.5384693101056831;
  GQ5[4] = 0.9061798459386640;
  weight5[0] = 0.2369268850561891;
  weight5[1] = 0.4786286704993665;
  weight5[2] = 0.5688888888888889;
  weight5[3] = 0.4786286704993665;
  weight5[4] = 0.2369268850561891;
  
  
  //More Gaussian weights and nodes:
  for (int i = 1; i <= 10; i = i + 1)
    {
      weight20[i + 9] = weight20[10 - i];
      GQ20[i + 9] = -1 * GQ20[10 - i];
    }
  for (int i = 1; i <= 5; i = i + 1)
    {
      weight10[i + 4] = weight10[5 - i];
      GQ10[i + 4] = -1 * GQ10[5 - i];
    }
  for (int i = 1; i <= 3; i = i + 1)
    {
      weight6[i + 2] = weight6[3 - i];
      GQ6[i + 2] = -1 * GQ6[3 - i];
    }

for (int i = 1; i <= 2; i = i + 1)
    {
      weight4[i + 1] = weight4[2 - i];
      GQ4[i + 1] = -1 * GQ4[2 - i];
    }

for (int i = 1; i <= 1; i = i + 1)
    {
      weight2[i + 0] = weight2[1 - i];
      GQ_loc_2[i + 0] = -1 * GQ_loc_2[1 - i];
    }
  switch (dim)
    {
    case 2: 
      {
	output.loc = GQ_loc_2[spot];
	output.weight = weight2[spot];
	return output;
	break;
      }
    case 3: 
      {
	output.loc = GQ_loc_3[spot];
	output.weight = weight3[spot];
	return output;
	break;
      }
    case 4: 
      {
	output.loc = GQ4[spot];
	output.weight = weight4[spot];
	return output;
	break;
      }
    case 5: 
      {
	output.loc = GQ5[spot];
	output.weight = weight5[spot];
	return output;
	break;
      }
    case 6: 
      {
	output.loc = GQ6[spot];
	output.weight = weight6[spot];
	return output;
	break;
      }
    case 10:
      {
	output.loc = GQ10[spot];
	output.weight = weight10[spot];
	return output;
	break;
      }
    case 20:
      {
	output.loc = GQ20[spot];
	output.weight = weight20[spot];
	return output;
	break;
      }
    }
}
double gquad(double f_data[GQdeg]) //quadrature for unit domain
{
  double sum = f_data[0] * weights[0];
  for (int j = 1; j < GQdeg - 1; j = j + 1)
    {
      sum = sum + weights[j] * f_data[j];
    }
  return 0.5 * (sum + f_data[GQdeg - 1] * weights[GQdeg - 1]);
}
vec4 vec_gquad(vec4 values[GQdeg]) //vectored quadrature, unit domain
{
  vec4 output;
  double mass[GQdeg];
  double mo_x[GQdeg];
  double mo_y[GQdeg];
  double energy[GQdeg];
  for (int i = 0; i < GQdeg; i = i + 1)
    {
      mass[i] = values[i].ma;
      mo_x[i] = values[i].mx;
      mo_y[i] = values[i].my;
      energy[i] = values[i].en;
    }
  output.ma = gquad(mass);
  output.mx = gquad(mo_x);
  output.my = gquad(mo_y);
  output.en = gquad(energy);
  return output;
}
double harmonic(double uL, double uM, double uR) //limiter function for scalar variable
{
  double dL = uM-uL;
  double dR = uR-uM;
  double magL = fabs(dL);
  double magR = fabs(dR);
  return (dL*magR + dR*magL)/(magL+magR+eps);
}
vec4 vec_harmonic(vec4 UL, vec4 UM, vec4 UR) //limiter function for vector.
{
  vec4 output;
  output.ma = harmonic(UL.ma, UM.ma, UR.ma);
  output.mx = harmonic(UL.mx, UM.mx, UR.mx);
  output.my = harmonic(UL.my, UM.my, UR.my);
  output.en = harmonic(UL.en, UM.en, UR.en);
  return output;
}
vec4 set_vec(double input) //sets a vector to the input scalar
{
  vec4 output;
  output.ma = input;
  output.mx = input;
  output.my = input;
  output.en = input;
  return output;
}
vec4 f2vec(double mass, double mo_x, double mo_y, double energy) //sets a vector to a vector input
{
  vec4 output;
  output.ma = mass;
  output.mx = mo_x;
  output.my = mo_y;
  output.en = energy;
  return output;
}
//General math operations on the vec4 structure:
vec4 vec2vec(vec4 input)
{
  vec4 output;
  output.ma = input.ma;
  output.mx = input.mx;
  output.my = input.my;
  output.en = input.en;
  return output;
}
vec4 vecdiff(vec4 lead, vec4 trail)
{
  vec4 output;
  output.ma = lead.ma - trail.ma;
  output.mx = lead.mx - trail.mx;
  output.my = lead.my - trail.my;
  output.en = lead.en - trail.en;
  return output;
}
vec4 vecdot(vec4 lead, vec4 trail)
{
  vec4 output;
  output.ma = lead.ma * trail.ma;
  output.mx = lead.mx * trail.mx;
  output.my = lead.my * trail.my;
  output.en = lead.en * trail.en;
  return output;
}
vec4 scavec(double scalar, vec4 trail)
{
  vec4 output;
  output.ma = scalar*trail.ma;
  output.mx = scalar*trail.mx;
  output.my = scalar*trail.my;
  output.en = scalar*trail.en;
  return output;
}
vec4 vecadd(vec4 lead, vec4 trail)
{
  vec4 output;
  output.ma = lead.ma + trail.ma;
  output.mx = lead.mx + trail.mx;
  output.my = lead.my + trail.my;
  output.en = lead.en + trail.en;
  return output;
}
vec4 vecdiv(vec4 lead, vec4 trail)
{
  vec4 output;
  output.ma = lead.ma / trail.ma;
  output.mx = lead.mx / trail.mx;
  output.my = lead.my / trail.my;
  output.en = lead.en / trail.en;
  return output;
}
double trans_up2(double in, double a, double b) //sends reference in [-1,1] to physical coordinate in interval [a,b]
{
  return (in + 1) * (b - a) / (2.0) + a;
}
double trans_up1(double in, double a, double b) //sends reference in [0,1] to physical coordinate in interval [a,b]
{
  return in * (b - a) + a;
}
//Some items specific to the Euler equations in 2D:
double get_prop(vec4 input) //Fastest propogation speed at a point where U=input.
{
  double rho = input.ma;
  double u = input.mx / rho;
  double v = input.my / rho;
  double sq_Vmag = u*u + v*v;
  double p = (rsh-1.0) * (input.en - 0.5*rho*sq_Vmag);
  return sqrt(rsh*p/rho) + sqrt(sq_Vmag);
}
double g_of_p(double AK, double BK, double pres) //Handy function for the HLLE Riemann solver.
{
  return sqrt( AK/(pres+BK));
}
vec4 Egov(vec4 U) //Egov is the function in the dE/dx portion of the GDE
{
  vec4 E;
  double rho = U.ma;
  double u = U.mx/rho;
  double v = U.my/rho;
  double sq_Vmag = u*u + v*v;
  double p = (rsh-1.0) * (U.en - 0.5*rho*sq_Vmag);
  E.ma = U.mx;
  E.mx = u*U.mx + p;
  E.my = v*U.mx;
  E.en = u*(U.en + p);
  return E;
}
vec4 Fgov(vec4 U) //Fgov is the function in the dF/dy portion of the GDE
{
  vec4 F;
  double rho = U.ma;
  double u = U.mx/rho;
  double v = U.my/rho;
  double sq_Vmag = u*u + v*v;
  double p = (rsh-1.0) * (U.en - 0.5*rho*sq_Vmag);
  F.ma = U.my;
  F.mx = u*U.my;
  F.my = v*U.my + p;
  F.en = v*(U.en + p);
  return F;
}
vec4 Ugov(double density, double vel_X, double vel_Y, double pressure) //Gets vector of conserved variables from primitives.
{
  double sq_Vmag = vel_X * vel_X + vel_Y * vel_Y;
  vec4 output;
  output.ma = density;
  output.mx = density*vel_X;
  output.my = density*vel_Y;
  output.en = pressure/(rsh-1.0) + 0.5 * density * sq_Vmag;
  return output;
}

vec4 E_flux_Riemann(vec4 UL, vec4 UR) //Exact Riemann solver, flux in x direction
{
  int c_limit;
  vec4 output;
  double change_NR = 1.0;
  double sol_NR;
  double dsol_NR;
  double sol_uL;
  double sol_uR;
  double dsol_uL;
  double dsol_uR;
  double p_old;
  //Interface primitives:
  double u_flux;
  double p_flux;
  double rho_flux;
  double v_flux;
  double Eg_flux;
  //Left state primitives:
  double rhoL = UL.ma;
  double uL = UL.mx/rhoL;
  double vL = UL.my/rhoL;
  double sq_VmagL = uL*uL + vL*vL;
  double pL = (rsh-1.0) * (UL.en - 0.5*rhoL*sq_VmagL);
  double aL = sqrt(rsh*pL/rhoL);
  //Right state primitives
  double rhoR = UR.ma;
  double uR = UR.mx/rhoR;
  double vR = UR.my/rhoR;
  double sq_VmagR = uR*uR + vR*vR;
  double pR = (rsh-1.0) * (UR.en - 0.5*rhoR*sq_VmagR);
  double aR = sqrt(rsh*pR/rhoR);
  //Star State Variables
  double u_star;
  double p_star;
  int path;
  //Identify star state
  double sol_p = 0.5*(pL+pR); //initial pressure guess
  while (change_NR > TOL_NR)
    {
      c_limit = c_limit + 1;
      //Get get u_star of right side based on (u+a) wave from left
      if (sol_p > pR)  //shock
	{
	  //sol_uR = uR + (sol_p - pR) / pow( rhoR*(Chi*pR + sol_p) / (1.0-Chi) , 0.5 );
     	  //dsol_uR = pow((1.0-Chi)/rhoR, 0.5) / (Chi*pR+sol_p) * ( pow(Chi*pR+sol_p, 0.5) - 0.5*(sol_p-pR)/pow(Chi*pR+sol_p, 0.5) );
	  sol_uR = uR + (sol_p - pR) / sqrt(rhoR*(Chi*pR + sol_p) / (1.0-Chi));
     	  dsol_uR = sqrt((1.0-Chi)/rhoR) / (Chi*pR+sol_p) * ( sqrt(Chi*pR+sol_p) - 0.5*(sol_p-pR)/sqrt(Chi*pR+sol_p) );
	}
      else //rarefaction
	{
	  sol_uR = uR + 2.0*aR/(rsh-1.0) * ( pow((sol_p/pR), ((rsh-1.0)/(2.0*rsh)) ) - 1.0);
	  dsol_uR = 2.0*aR/(rsh-1.0) * (rsh-1.0)/(2.0*rsh* pow(pR,((rsh-1.0)/(2.0*rsh)))) * pow(sol_p, ((-rsh-1.0)/(2.0*rsh)) );
	}
      
      //Get u_star of left side based on (u-a) wave from right
      if (sol_p > pL)  //shock
	{
	  //sol_uL = uL + (pL - sol_p) / sqrt(rhoL*(Chi*pL + sol_p) / (1.0-Chi));
	  //dsol_uL = sqrt((1.0-Chi)/rhoL) / (Chi*pL+sol_p) * ( -sqrt(Chi*pL+sol_p) - 0.5*(pL-sol_p)/sqrt(Chi*pL+sol_p) );
	   sol_uL = uL + (pL - sol_p) / sqrt(rhoL*(Chi*pL + sol_p) / (1.0-Chi));
	   dsol_uL = sqrt((1.0-Chi)/rhoL) / (Chi*pL+sol_p) * ( -sqrt(Chi*pL+sol_p) - 0.5*(pL-sol_p)/sqrt(Chi*pL+sol_p) );
	}
      else //rarefaction
	{
	  sol_uL = uL + 2.0*aL/(rsh-1.0) * ( 1.0 - pow(sol_p/pL, (rsh-1.0)/(2.0*rsh)) );
	  dsol_uL = 2.0*aL/(rsh-1.0) * (1.0-rsh)/(2.0*rsh*pow(pL,((rsh-1.0)/(2.0*rsh)))) *  pow(sol_p, ((-rsh-1.0)/(2.0*rsh)));
        }
      
      sol_NR = sol_uL - sol_uR;
      dsol_NR = dsol_uL - dsol_uR;
         
      p_old = sol_p;
      sol_p = fabs(p_old - sol_NR/dsol_NR);
      change_NR = 2.0 * fabs(sol_p - p_old) / fabs(p_old + sol_p);
    } //End of NewtonRhapson iterations
      
  p_star = sol_p;
  //Grab velocity from final density solution
  if (p_star > pR)  //shock
    {
      u_star = uR + (p_star - pR) / sqrt(rhoR*(Chi*pR + p_star) / (1.0-Chi));
    }
  else //rarefaction
    {
      u_star = uR + 2.0*aR/(rsh-1.0) * (pow(p_star/pR, ((rsh-1.0)/(2.0*rsh))) -1.0);
    }
  //p_star and u_star are now known
  //Determine path:
  //1: positive u_star, left wave is shock
  //2: positive u_star, left wave is rarefaction
  //3: negative u_star, right wave is shock
  //4: negative u_star, right wave is rarefaction
  if (u_star > 0) //contact is right of interfac
    {
      path = 1; //default decision is left shock
      if (p_star < pL) //left rarefaction
	{
	  path = 2;
	}
    }
  else //Contact is left of interface, or ustar=0
    {
      path = 3; //default decision is right shock
      if (p_star < pR) //right rarefaction
	{
	  path = 4;
	}
    }
  //Refine analysis based on path
  switch (path)
    {
    case 1: //need to analyze left shock
      {
	double rho_star_L = rhoL * (p_star/pL + Chi) / (Chi*p_star/pL + 1);
	//double vel_L = uL - aL * pow( (rsh+1.0)/(2.0*rsh) * p_star/pL + (rsh-1.0)/(2.0*rsh) , 0.5 );
	double vel_L = uL - aL * sqrt( (rsh+1.0)/(2.0*rsh) * p_star/pL + (rsh-1.0)/(2.0*rsh) );
	if (vel_L > 0.0) //left shock is right of interface
	  {
	    rho_flux = rhoL;
	    u_flux = uL;
	    v_flux = vL;
	    p_flux = pL;
	  }
	else //interface lies between left shock and contact
	  {
	    rho_flux = rho_star_L;
	    u_flux = u_star;
	    v_flux = vL;
	    p_flux = p_star;
	  }
	break;
      }
    case 2: //u-a wave is rarefaction
      {
	double rho_star_L = rhoL * pow(p_star/pL, 1.0/rsh);
	double aL_star = aL*pow(p_star/pL , (rsh-1.0)/(2.0*rsh) );
	double vel_HL = uL - aL;
	double vel_TL = u_star - aL_star;
	if (vel_HL > 0.0) //left head is right of interface
	  {
	    rho_flux = rhoL;
	    u_flux = uL;
	    v_flux = vL;
	    p_flux = pL;
	  }
	else  //left head is left of interface
	  {
	    if (vel_TL > 0.0) //left rarefaction straddles interface
	      {
		double ratio_base = 2.0/(rsh+1.0) + (rsh-1.0)/(aL*(rsh+1.0))*uL;
		rho_flux = rhoL * pow(ratio_base, 2.0/(rsh-1.0));
		u_flux = 2.0/(rsh+1.0) * (aL + 0.5*uL*(rsh-1.0) );
		v_flux = vL;
		p_flux = pL * pow( ratio_base , 2.0*rsh/(rsh-1.0) );
	      }
	    else //star state straddles interface
	      {
		rho_flux = rho_star_L;
		u_flux = u_star;
		v_flux = vL;
		p_flux = p_star;
	      }
	  }
	break;
      }
    case 3: //u+a wave is shock, contact is left of interface
      {
	double rho_star_R = rhoR * (p_star/pR + Chi) / (Chi*p_star/pR + 1);
	//double vel_R = uR + aR * pow( (rsh+1.0)/(2.0*rsh) * p_star/pR + (rsh-1.0)/(2.0*rsh) , 0.5);
	double vel_R = uR + aR * sqrt( (rsh+1.0)/(2.0*rsh) * p_star/pR + (rsh-1.0)/(2.0*rsh) );
	if (vel_R < 0.0) //right shock is left of interface
	  {
	    rho_flux = rhoR;
	    u_flux = uR;
	    v_flux = vR;
	    p_flux = pR;
	  }
	else //interface lies between right shock and contact
	  {
	    rho_flux = rho_star_R;
	    u_flux = u_star;
	    v_flux = vR;
	    p_flux = p_star;
	  }
	break;
      }
    case 4:
      {
	double rho_star_R = rhoR * pow(p_star/pR , 1.0/rsh);
	double aR_star = aR * pow(p_star/pR , (rsh-1.0)/(2.0*rsh));
	double vel_HR = uR + aR;
	double vel_TR = u_star + aR_star;
	if (vel_HR < 0.0) //right head is left of interface
	  {
	    rho_flux = rhoR;
	    u_flux = uR;
	    v_flux = vR;
	    p_flux = pR;
	  }
	else  //right head is right of interface
	  {
	    if (vel_TR < 0.0) //right rarefaction straddles interface
	      {
		double ratio_base = 2.0/(rsh+1.0) - (rsh-1.0)/(aR*(rsh+1.0))*uR;
		rho_flux = rhoR * pow(ratio_base, 2.0/(rsh-1.0));
		u_flux = 2.0/(rsh+1.0) * (-aR + 0.5*uR*(rsh-1.0) );
		v_flux = vR;
		p_flux = pR * pow(ratio_base, 2.0*rsh/(rsh-1.0));
	      }
	    else //star state straddles interface
	      {
		rho_flux = rho_star_R;
		u_flux = u_star;
		v_flux = vR;
		p_flux = p_star;
	      }
	  }
	break;
      }
    }
  //Use interface values to define interface flux E:
  Eg_flux = 1.0/rho_flux * ( p_flux/(rsh-1.0) + 0.5*rho_flux*(u_flux*u_flux + v_flux*v_flux) );
  output.ma = rho_flux*u_flux;
  output.mx = rho_flux*u_flux*u_flux + p_flux;
  output.my = rho_flux*u_flux*v_flux;
  output.en = u_flux * (rho_flux * Eg_flux + p_flux);
  return output;
}
vec4 F_flux_Riemann(vec4 UL, vec4 UR) //Exact Riemann solver, flux in y direction
{
  //L is for low, R is for roof. 
  int c_limit;
  vec4 output;
  double change_NR = 1.0;
  double sol_NR;
  double dsol_NR;
  double sol_vL;
  double sol_vR;
  double dsol_vL;
  double dsol_vR;
  double p_old;
  //Interface primitives:
  double u_flux;
  double p_flux;
  double rho_flux;
  double v_flux;
  double Eg_flux;
  //Left state primitives:
  double rhoL = UL.ma;
  double uL = UL.mx/rhoL;
  double vL = UL.my/rhoL;
  double sq_VmagL = uL*uL + vL*vL;
  double pL = (rsh-1.0) * (UL.en - 0.5*rhoL*sq_VmagL);
  double aL = sqrt(rsh*pL/rhoL);
  //Right state primitives
  double rhoR = UR.ma;
  double uR = UR.mx/rhoR;
  double vR = UR.my/rhoR;
  double sq_VmagR = uR*uR + vR*vR;
  double pR = (rsh-1.0) * (UR.en - 0.5*rhoR*sq_VmagR);
  double aR = sqrt(rsh*pR/rhoR);
  //Star State Variables
  double v_star;
  double p_star;
  int path;
  //Identify star state
  double sol_p = 0.5*(pL+pR); //initial pressure guess
  while (change_NR > TOL_NR)
    {
      c_limit = c_limit + 1;
      //Get get v_star of right side based on (u+a) wave from floor
      if (sol_p > pR)  //shock
	{
	  //sol_vR = vR + (sol_p - pR) / pow( rhoR*(Chi*pR + sol_p) / (1.0-Chi) , 0.5 );
     	  //dsol_vR = pow((1.0-Chi)/rhoR, 0.5) / (Chi*pR+sol_p) * ( pow(Chi*pR+sol_p, 0.5) - 0.5*(sol_p-pR)/pow(Chi*pR+sol_p, 0.5) );
	  sol_vR = vR + (sol_p - pR) / sqrt( rhoR*(Chi*pR + sol_p) / (1.0-Chi) );
     	  dsol_vR = sqrt((1.0-Chi)/rhoR) / (Chi*pR+sol_p) * ( sqrt(Chi*pR+sol_p) - 0.5*(sol_p-pR)/sqrt(Chi*pR+sol_p) );
	}
      else //rarefaction
	{
	  sol_vR = vR + 2.0*aR/(rsh-1.0) * ( pow((sol_p/pR), ((rsh-1.0)/(2.0*rsh)) ) - 1.0);
	  dsol_vR = 2.0*aR/(rsh-1.0) * (rsh-1.0)/(2.0*rsh* pow(pR,((rsh-1.0)/(2.0*rsh)))) * pow(sol_p, ((-rsh-1.0)/(2.0*rsh)) );
	}
      
      //Get v_star of left side based on (v-a) wave from roof
      if (sol_p > pL)  //shock
	{
	  sol_vL = vL + (pL - sol_p) / sqrt(rhoL*(Chi*pL + sol_p) / (1.0-Chi));
	  dsol_vL = sqrt((1.0-Chi)/rhoL) / (Chi*pL+sol_p) * ( -sqrt(Chi*pL+sol_p) - 0.5*(pL-sol_p)/sqrt(Chi*pL+sol_p) );
	}
      else //rarefaction
	{
	  sol_vL = vL + 2.0*aL/(rsh-1.0) * ( 1.0 - pow(sol_p/pL, (rsh-1.0)/(2.0*rsh)) );
	  dsol_vL = 2.0*aL/(rsh-1.0) * (1.0-rsh)/(2.0*rsh*pow(pL,((rsh-1.0)/(2.0*rsh)))) *  pow(sol_p, ((-rsh-1.0)/(2.0*rsh)));
        }
      
      sol_NR = sol_vL - sol_vR;
      dsol_NR = dsol_vL - dsol_vR;
         
      p_old = sol_p;
      sol_p = fabs(p_old - sol_NR/dsol_NR);
      change_NR = 2.0 * fabs(sol_p - p_old) / fabs(p_old + sol_p);
    } //End of NewtonRhapson iterations
      
  p_star = sol_p;
  //Grab velocity from final density solution
  if (p_star > pR)  //shock
    {
      v_star = vR + (p_star - pR) / sqrt(rhoR*(Chi*pR + p_star) / (1.0-Chi));
    }
  else //rarefaction
    {
      v_star = vR + 2.0*aR/(rsh-1.0) * (pow(p_star/pR, ((rsh-1.0)/(2.0*rsh))) -1.0);
    }
  //p_star and u_star are now known
  //Determine path:
  //1: positive v_star, low wave is shock
  //2: positive v_star, low wave is rarefaction
  //3: negative v_star, high wave is shock
  //4: negative v_star, high wave is rarefaction
  if (v_star > 0) //contact is high of interfac
    {
      path = 1; //default decision is low shock
      if (p_star < pL) //low rarefaction
	{
	  path = 2;
	}
    }
  else //Contact is low of interface, or ustar=0
    {
      path = 3; //default decision is right shock
      if (p_star < pR) //high rarefaction
	{
	  path = 4;
	}
    }
  //Refine analysis based on path
  switch (path)
    {
    case 1: //need to analyze low shock
      {
	double rho_star_L = rhoL * (p_star/pL + Chi) / (Chi*p_star/pL + 1);
	//double vel_L = vL - aL * pow( (rsh+1.0)/(2.0*rsh) * p_star/pL + (rsh-1.0)/(2.0*rsh) , 0.5 );
	double vel_L = vL - aL * sqrt( (rsh+1.0)/(2.0*rsh) * p_star/pL + (rsh-1.0)/(2.0*rsh) );
	if (vel_L > 0.0) //low shock is high of interface
	  {
	    rho_flux = rhoL;
	    u_flux = uL;
	    v_flux = vL;
	    p_flux = pL;
	  }
	else //interface lies between low shock and contact
	  {
	    rho_flux = rho_star_L;
	    u_flux = uL;
	    v_flux = v_star;
	    p_flux = p_star;
	  }
	break;
      }
    case 2: //v-a wave is rarefaction
      {
	double rho_star_L = rhoL * pow(p_star/pL, 1.0/rsh);
	double aL_star = aL*pow(p_star/pL , (rsh-1.0)/(2.0*rsh) );
	double vel_HL = vL - aL;
	double vel_TL = v_star - aL_star;
	if (vel_HL > 0.0) //low head is high of interface
	  {
	    rho_flux = rhoL;
	    u_flux = uL;
	    v_flux = vL;
	    p_flux = pL;
	  }
	else  //low head is low of interface
	  {
	    if (vel_TL > 0.0) //low rarefaction straddles interface
	      {
		double ratio_base = 2.0/(rsh+1.0) + (rsh-1.0)/(aL*(rsh+1.0))*vL;
		rho_flux = rhoL * pow(ratio_base, 2.0/(rsh-1.0));
		v_flux = 2.0/(rsh+1.0) * (aL + 0.5*vL*(rsh-1.0) );
		u_flux = uL;
		p_flux = pL * pow( ratio_base , 2.0*rsh/(rsh-1.0) );
	      }
	    else //star state straddles interface
	      {
		rho_flux = rho_star_L;
		u_flux = uL;
		v_flux = v_star;
		p_flux = p_star;
	      }
	  }
	break;
      }
    case 3: //u+a wave is shock, contact is low of interface
      {
	double rho_star_R = rhoR * (p_star/pR + Chi) / (Chi*p_star/pR + 1);
	//double vel_R = vR + aR * pow( (rsh+1.0)/(2.0*rsh) * p_star/pR + (rsh-1.0)/(2.0*rsh) , 0.5);
	double vel_R = vR + aR * sqrt( (rsh+1.0)/(2.0*rsh) * p_star/pR + (rsh-1.0)/(2.0*rsh) );
	if (vel_R < 0.0) //high shock is low of interface
	  {
	    rho_flux = rhoR;
	    u_flux = uR;
	    v_flux = vR;
	    p_flux = pR;
	  }
	else //interface lies between high shock and contact
	  {
	    rho_flux = rho_star_R;
	    u_flux = uR;
	    v_flux = v_star;
	    p_flux = p_star;
	  }
	break;
      }
    case 4:
      {
	double rho_star_R = rhoR * pow(p_star/pR , 1.0/rsh);
	double aR_star = aR * pow(p_star/pR , (rsh-1.0)/(2.0*rsh));
	double vel_HR = vR + aR;
	double vel_TR = v_star + aR_star;
	if (vel_HR < 0.0) //high head is low of interface
	  {
	    rho_flux = rhoR;
	    u_flux = uR;
	    v_flux = vR;
	    p_flux = pR;
	  }
	else  //high head is high of interface
	  {
	    if (vel_TR < 0.0) //high rarefaction straddles interface
	      {
		double ratio_base = 2.0/(rsh+1.0) - (rsh-1.0)/(aR*(rsh+1.0))*vR;
		rho_flux = rhoR * pow(ratio_base, 2.0/(rsh-1.0));
		v_flux = 2.0/(rsh+1.0) * (-aR + 0.5*vR*(rsh-1.0) );
		u_flux = uR;
		p_flux = pR * pow(ratio_base, 2.0*rsh/(rsh-1.0));
	      }
	    else //star state straddles interface
	      {
		rho_flux = rho_star_R;
		u_flux = uR;
		v_flux = v_star;
		p_flux = p_star;
	      }
	  }
	break;
      }
    }
  //Use interface values to define interface flux F:
  Eg_flux = 1.0/rho_flux * ( p_flux/(rsh-1.0) + 0.5*rho_flux*(u_flux*u_flux + v_flux*v_flux) );
  output.ma = rho_flux * v_flux;
  output.mx = rho_flux * u_flux * v_flux;
  output.my = rho_flux * v_flux * v_flux + p_flux;
  output.en = v_flux * (rho_flux * Eg_flux + p_flux);
  return output;
}
//These next four functions are for the Roe Riemann solver
vec4 populate_K0_E(double u, double v, double Vsq, double H, double sos)
{
  vec4 output;
  output.ma = 1.0;
  output.mx = u-sos;
  output.my = v;
  output.en = H - u*sos;
  return output;
}
vec4 populate_K1_E(double u, double v, double Vsq, double H, double sos)
{
  vec4 output;
  output.ma = 1.0;
  output.mx = u;
  output.my = v;
  output.en = 0.5*Vsq;
  return output;
}
vec4 populate_K2_E(double u, double v, double Vsq, double H, double sos)
{
  vec4 output;
  output.ma = 0.0;
  output.mx = 0.0;
  output.my = 1.0;
  output.en = v;
  return output;
}
vec4 populate_K3_E(double u, double v, double Vsq, double H, double sos)
{
  vec4 output;
  output.ma = 1.0;
  output.mx = u+sos;
  output.my = v;
  output.en = H + u*sos;
  return output;
}
//Roe solvers for x-direction flux. RoeHart uses entropy correction
vec4 E_flux_Roe(vec4 UL, vec4 UR)
{
  //Flux_E_Timer.startBlock();
  //Step 0: Extract information for UL, UR---------------- (very costly)
  double rho_L = UL.ma;
  double u_L = UL.mx / rho_L;
  double v_L = UL.my / rho_L;
  double Vsq_L = u_L*u_L + v_L*v_L;
  double p_L = (rsh-1.0) * (UL.en - 0.5*rho_L*Vsq_L);
  double sos_L = sqrt(rsh*p_L/rho_L);
  double H_L = (UL.en + p_L) / rho_L;
  
  double rho_R = UR.ma;
  double u_R = UR.mx / rho_R;
  double v_R = UR.my / rho_R;
  double Vsq_R = u_R*u_R + v_R*v_R;
  double p_R = (rsh-1.0) * (UR.en - 0.5*rho_R*Vsq_R);
  double sos_R = sqrt(rsh*p_R/rho_R);
  double H_R = (UR.en + p_R) / rho_R;
  

  
  //Step 1: Get Roe-averaged variables------------------
  double rho_avg = sqrt(rho_L * rho_R);
  double sqrt_rho_L = sqrt(rho_L);
  double sqrt_rho_R = sqrt(rho_R);
  double div = sqrt_rho_L + sqrt_rho_R;
  double u_avg = (sqrt_rho_L*u_L + sqrt_rho_R*u_R) / div;
  double v_avg = (sqrt_rho_L*v_L + sqrt_rho_R*v_R) / div;
  double Vsq_avg = u_avg*u_avg + v_avg*v_avg;
  double H_avg = (sqrt_rho_L*H_L + sqrt_rho_R*H_R) / div;
  double sos_avg = sqrt( (rsh-1.0) * (H_avg-0.5*Vsq_avg) );
  

  
  //Step 2:Define wave speeds in X direction based on average state----------
  double lam[4]; //wave speeds
  lam[0] = u_avg - sos_avg;
  lam[1] = sos_avg;
  lam[2] = sos_avg; 
  lam[3] = u_avg + sos_avg;
  
  //Step 3: Populate Right Eigenvectors------------------------
  vec4 K_eigen[4];
  K_eigen[0] = populate_K0_E(u_avg, v_avg, Vsq_avg, H_avg, sos_avg); 
  K_eigen[1] = populate_K1_E(u_avg, v_avg, Vsq_avg, H_avg, sos_avg); 
  K_eigen[2] = populate_K2_E(u_avg, v_avg, Vsq_avg, H_avg, sos_avg); 
  K_eigen[3] = populate_K3_E(u_avg, v_avg, Vsq_avg, H_avg, sos_avg);
  

  
  //Step 4: Get the wave strengths-------------------------------
  double del_p = p_R - p_L;
  double del_u = u_R - u_L;
  double del_v = v_R - v_L;
  double del_rho = rho_R - rho_L;
  double alpha[4]; //wave strengths
  double DA_du = rho_avg * sos_avg * del_u;
  double sos_sq = sos_avg * sos_avg;
  alpha[0] = 0.5/sos_sq * (del_p - DA_du);
  alpha[1] = del_rho - del_p/sos_sq;
  alpha[2] = rho_avg * del_v;
  alpha[3] = 0.5/sos_sq * (del_p + DA_du);
  

  
  //Step 5: Calculate wave-based perturbations, add to average flux------
  double scalar[4]; //speed*strength
  for (int i = 0; i < 4; i = i + 1)
    {
      scalar[i] = alpha[i] * fabs(lam[i]);
    }

  vec4 perts = set_vec(0.0); //perturbation, calculated by summation
  for (int i = 0; i < 4; i = i + 1)
    {
      perts = vecadd( perts, scavec(scalar[i], K_eigen[i]) );
    }
  

  
  vec4 flux_avg = scavec( 0.5 , vecadd(Egov(UL), Egov(UR)) );
  vec4 output;
  output = vecdiff( flux_avg, scavec(0.5, perts) );
  
  //Flux_E_Timer.endBlock();
  return output;
}
vec4 E_flux_RoeHart(vec4 UL, vec4 UR)
{
  
  vec4 output;
  //Step 0: Extract information for UL, UR----------------
  double rho_L = UL.ma;
  double u_L = UL.mx / rho_L;
  double v_L = UL.my / rho_L;
  double Vsq_L = u_L*u_L + v_L*v_L;
  double p_L = (rsh-1.0) * (UL.en - 0.5*rho_L*Vsq_L);
  double sos_L = sqrt(rsh*p_L/rho_L);
  double H_L = (UL.en + p_L) / rho_L;
  
  double rho_R = UR.ma;
  double u_R = UR.mx / rho_R;
  double v_R = UR.my / rho_R;
  double Vsq_R = u_R*u_R + v_R*v_R;
  double p_R = (rsh-1.0) * (UR.en - 0.5*rho_R*Vsq_R);
  double sos_R = sqrt(rsh*p_R/rho_R);
  double H_R = (UR.en + p_R) / rho_R;

  //Step 1: Get Roe-averaged variables------------------
  double rho_avg = sqrt(rho_L * rho_R);
  double sqrt_rho_L = sqrt(rho_L);
  double sqrt_rho_R = sqrt(rho_R);
  double div = sqrt_rho_L + sqrt_rho_R;
  double u_avg = (sqrt_rho_L*u_L + sqrt_rho_R*u_R) / div;
  double v_avg = (sqrt_rho_L*v_L + sqrt_rho_R*v_R) / div;
  double Vsq_avg = u_avg*u_avg + v_avg*v_avg;
  double H_avg = (sqrt_rho_L*H_L + sqrt_rho_R*H_R) / div;
  double sos_avg = sqrt( (rsh-1.0) * (H_avg-0.5*Vsq_avg) );

  //Check sonic indicators
  double CL = rho_L*sos_L; //rho*sos at the foot of each cell's characteristics
  double CR = rho_R*sos_R;
  double pvrs = (CR*p_L + CL*p_R + CL*CR*(u_L - u_R)) / (CL+CR);
  double p_pre = fmax(pvrs,0.0); //predictor pressure for TSRS solution
  
  //TSRS solution, Toro absolute page 322
  double A_L = 2.0/(rho_L*(rsh+1.0));
  double A_R = 2.0/(rho_R*(rsh+1.0));
  double B_L = p_L * (rsh-1.0) / (rsh+1.0);
  double B_R = p_R * (rsh-1.0) / (rsh+1.0);
  double g_L = g_of_p(A_L, B_L, p_pre);
  double g_R = g_of_p(A_R, B_R, p_pre);
  double p_star = (g_L*p_L + g_R*p_R - u_R + u_L) / (g_L + g_R);
  
  double u_star = 0.5 * (u_L+u_R + (p_star-p_R)*g_R - (p_star-p_L)*g_L);
  double rhoL_star = rho_L *  (p_star/p_L + G6) / (G6*p_star/p_L + 1) ;
  double rhoR_star = rho_R *  (p_star/p_R + G6) / (G6*p_star/p_R + 1) ;
  double aL_star = sqrt(rsh*p_L / rhoL_star);
  double aR_star = sqrt(rsh*p_R / rhoR_star);
  
  //Wavespeeds for Harten decision (First interpretation of aL_star, aR_star was wrong)
  double S0_L = u_L - sos_L;
  double S3_R = u_R + sos_R;
  double S0_Rstar = u_star - aR_star;
  double S3_Lstar = u_star + aL_star;

  if (S0_L < 0 && 0 < S0_Rstar)
    {
      //Left sonic fix
      //u-a characteristic perturbs left side
      double lam_0 = u_avg - sos_avg;

      //wave strength and eigenvector:
      double del_p = p_R - p_L;
      double del_u = u_R - u_L;
      double del_v = v_R - v_L;
      double del_rho = rho_R - rho_L;
      
      double DA_du = rho_avg * sos_avg * del_u;
      double sos_sq = sos_avg * sos_avg;
      double alpha_0 = 0.5/sos_sq * (del_p - DA_du);

      vec4 K_eigen_0 = populate_K0_E(u_avg, v_avg, Vsq_avg, H_avg, sos_avg);

      //Get fix for the lam_0 wavespeed, which is currently based on Roe averages
      double lam_fix = S0_L * (S0_Rstar - lam_0) / (S0_Rstar - S0_L);
      
      //get perturbation
      vec4 perts = scavec( lam_fix*alpha_0, K_eigen_0);
      output = vecadd(Egov(UL) , perts);
    }
  else if (S3_Lstar < 0 && 0 < S3_R )
    {
      //Right sonic fix
      //u+a characteristic perturbs right side
      double lam_3 = u_avg + sos_avg;

      //wave strength and eigenvector:
      double del_p = p_R - p_L;
      double del_u = u_R - u_L;
      double del_v = v_R - v_L;
      double del_rho = rho_R - rho_L;
      
      double DA_du = rho_avg * sos_avg * del_u;
      double sos_sq = sos_avg * sos_avg;
      double alpha_3 = 0.5/sos_sq * (del_p - DA_du);

      vec4 K_eigen_3 = populate_K3_E(u_avg, v_avg, Vsq_avg, H_avg, sos_avg);

      //Get fix for the lam_3 wavespeed
      double lam_fix = S3_R * (lam_3 - S3_Lstar) / (S3_R - S3_Lstar);

      vec4 perts = scavec(lam_fix*alpha_3 , K_eigen_3);
      output = vecdiff(Egov(UR) , perts);
    }
  else
    {
      //Regular RoePike procedure
      //Step 2:Define wave speeds in X direction based on average state----------
      double lam[4]; //wave speeds
      lam[0] = u_avg - sos_avg;
      lam[1] = sos_avg;
      lam[2] = sos_avg; 
      lam[3] = u_avg + sos_avg;
      
      //Step 3: Populate Right Eigenvectors------------------------
      vec4 K_eigen[4];
      K_eigen[0] = populate_K0_E(u_avg, v_avg, Vsq_avg, H_avg, sos_avg); 
      K_eigen[1] = populate_K1_E(u_avg, v_avg, Vsq_avg, H_avg, sos_avg); 
      K_eigen[2] = populate_K2_E(u_avg, v_avg, Vsq_avg, H_avg, sos_avg); 
      K_eigen[3] = populate_K3_E(u_avg, v_avg, Vsq_avg, H_avg, sos_avg);
      
      //Step 4: Get the wave strengths-------------------------------
      double del_p = p_R - p_L;
      double del_u = u_R - u_L;
      double del_v = v_R - v_L;
      double del_rho = rho_R - rho_L;
      double alpha[4]; //wave strengths
      double DA_du = rho_avg * sos_avg * del_u;
      double sos_sq = sos_avg * sos_avg;
      alpha[0] = 0.5/sos_sq * (del_p - DA_du);
      alpha[1] = del_rho - del_p/sos_sq;
      alpha[2] = rho_avg * del_v;
      alpha[3] = 0.5/sos_sq * (del_p + DA_du);
      
      //Step 5: Calculate wave-based perturbations, add to average flux------
      double scalar[4]; //speed*strength
      for (int i = 0; i < 4; i = i + 1)
	{
	  scalar[i] = alpha[i] * fabs(lam[i]);
	}
      
      vec4 perts = set_vec(0.0); //perturbation, calculated by summation
      for (int i = 0; i < 4; i = i + 1)
	{
	  perts = vecadd( perts, scavec(scalar[i], K_eigen[i]) );
	}
      vec4 flux_avg = scavec( 0.5 , vecadd(Egov(UL), Egov(UR)) );
      output = vecdiff( flux_avg, scavec(0.5, perts) );
    }
  
  return output;
  
}
//These next four functions are for the Roe Riemann solver
vec4 populate_K0_F(double u, double v, double Vsq, double H, double sos)
{
  vec4 output;
  output.ma = 1.0;
  output.mx = u;
  output.my = v-sos;
  output.en = H - v*sos;
  return output;
}
vec4 populate_K1_F(double u, double v, double Vsq, double H, double sos)
{
  vec4 output;
  output.ma = 0.0;
  output.mx = 1.0;
  output.my = 0.0;
  output.en = u;
  return output;
}
vec4 populate_K2_F(double u, double v, double Vsq, double H, double sos)
{
  vec4 output;
  output.ma = 1.0;
  output.mx = u;
  output.my = v;
  output.en = 0.5*Vsq;
  return output;
}
vec4 populate_K3_F(double u, double v, double Vsq, double H, double sos)
{
  vec4 output;
  output.ma = 1.0;
  output.mx = u;
  output.my = v+sos;
  output.en = H + v*sos;
  return output;
}
//Roe solvers for x-direction flux. RoeHart uses entropy correction
vec4 F_flux_Roe(vec4 UL, vec4 UR)
{
  //Flux_F_Timer.startBlock();
  //Step 0: Extract information for UL, UR----------------
  double rho_L = UL.ma;
  double u_L = UL.mx / rho_L;
  double v_L = UL.my / rho_L;
  double Vsq_L = u_L*u_L + v_L*v_L;
  double p_L = (rsh-1.0) * (UL.en - 0.5*rho_L*Vsq_L);
  double sos_L = sqrt(rsh*p_L/rho_L);
  double H_L = (UL.en + p_L) / rho_L;
  
  double rho_R = UR.ma;
  double u_R = UR.mx / rho_R;
  double v_R = UR.my / rho_R;
  double Vsq_R = u_R*u_R + v_R*v_R;
  double p_R = (rsh-1.0) * (UR.en - 0.5*rho_R*Vsq_R);
  double sos_R = sqrt(rsh*p_R/rho_R);
  double H_R = (UR.en + p_R) / rho_R;

  //Step 1: Get Roe-averaged variables------------------
  double rho_avg = sqrt(rho_L * rho_R);
  double sqrt_rho_L = sqrt(rho_L);
  double sqrt_rho_R = sqrt(rho_R);
  double div = sqrt_rho_L + sqrt_rho_R;
  double u_avg = (sqrt_rho_L*u_L + sqrt_rho_R*u_R) / div;
  double v_avg = (sqrt_rho_L*v_L + sqrt_rho_R*v_R) / div;
  double Vsq_avg = u_avg*u_avg + v_avg*v_avg;
  double H_avg = (sqrt_rho_L*H_L + sqrt_rho_R*H_R) / div;
  double sos_avg = sqrt( (rsh-1.0) * (H_avg-0.5*Vsq_avg) );

  //Step 2:Define wave speeds in Y direction based on average state----------
  double lam[4]; //wave speeds
  lam[0] = v_avg - sos_avg;
  lam[1] = sos_avg;
  lam[2] = sos_avg; 
  lam[3] = v_avg + sos_avg;
  
  //Step 3: Populate Right Eigenvectors------------------------
  vec4 K_eigen[4];
  K_eigen[0] = populate_K0_F(u_avg, v_avg, Vsq_avg, H_avg, sos_avg); 
  K_eigen[1] = populate_K1_F(u_avg, v_avg, Vsq_avg, H_avg, sos_avg); 
  K_eigen[2] = populate_K2_F(u_avg, v_avg, Vsq_avg, H_avg, sos_avg); 
  K_eigen[3] = populate_K3_F(u_avg, v_avg, Vsq_avg, H_avg, sos_avg);

  //Step 4: Get the wave strengths-------------------------------
  double del_p = p_R - p_L;
  double del_u = u_R - u_L;
  double del_v = v_R - v_L;
  double del_rho = rho_R - rho_L;
  double alpha[4]; //wave strengths
  double DA_dv = rho_avg * sos_avg * del_v;
  double sos_sq = sos_avg * sos_avg;
  alpha[0] = 0.5/sos_sq * (del_p - DA_dv);
  alpha[1] = rho_avg * del_u;
  alpha[2] = del_rho - del_p/sos_sq;
  alpha[3] = 0.5/sos_sq * (del_p + DA_dv);

  //Step 5: Calculate wave-based perturbations, add to average flux------
  double scalar[4]; //speed*strength
  for (int i = 0; i < 4; i = i + 1)
    {
      scalar[i] = alpha[i] * fabs(lam[i]);
    }

  vec4 perts = set_vec(0.0); //perturbation, calculated by summation
  for (int i = 0; i < 4; i = i + 1)
    {
      perts = vecadd( perts, scavec(scalar[i], K_eigen[i]) );
    }
  vec4 flux_avg = scavec( 0.5 , vecadd(Fgov(UL), Fgov(UR)) );
  vec4 output;
  output = vecdiff( flux_avg, scavec(0.5, perts) );
  //Flux_F_Timer.endBlock();
  return output; 
}
vec4 F_flux_RoeHart(vec4 UL, vec4 UR)
{
  
  vec4 output;
  //~ indicates a change from X orientation to Y orientation
  //Step 0: Extract information for UL, UR----------------
  
  double rho_L = UL.ma;
  double u_L = UL.mx / rho_L;
  double v_L = UL.my / rho_L;
  double Vsq_L = u_L*u_L + v_L*v_L;
  double p_L = (rsh-1.0) * (UL.en - 0.5*rho_L*Vsq_L);
  double sos_L = sqrt(rsh*p_L/rho_L);
  double H_L = (UL.en + p_L) / rho_L;
    
  double rho_R = UR.ma;
  double u_R = UR.mx / rho_R;
  double v_R = UR.my / rho_R;
  double Vsq_R = u_R*u_R + v_R*v_R;
  double p_R = (rsh-1.0) * (UR.en - 0.5*rho_R*Vsq_R);
  double sos_R = sqrt(rsh*p_R/rho_R);
  double H_R = (UR.en + p_R) / rho_R;
  
  //Step 1: Get Roe-averaged variables------------------
  double rho_avg = sqrt(rho_L * rho_R);
  double sqrt_rho_L = sqrt(rho_L);
  double sqrt_rho_R = sqrt(rho_R);
  double div = sqrt_rho_L + sqrt_rho_R;
  double u_avg = (sqrt_rho_L*u_L + sqrt_rho_R*u_R) / div;
  double v_avg = (sqrt_rho_L*v_L + sqrt_rho_R*v_R) / div;
  double Vsq_avg = u_avg*u_avg + v_avg*v_avg;
  double H_avg = (sqrt_rho_L*H_L + sqrt_rho_R*H_R) / div;
  double sos_avg = sqrt( (rsh-1.0) * (H_avg-0.5*Vsq_avg) );
  
  //Check sonic indicators
  double CL = rho_L*sos_L; //rho*sos at the foot of each cell's characteristics
  double CR = rho_R*sos_R;
  double pvrs = (CR*p_L + CL*p_R + CL*CR*(v_L - v_R)) / (CL+CR); //~
  double p_pre = fmax(pvrs,0.0); //predictor pressure for TSRS solution
  
  //TSRS solution, Toro absolute page 322
  double A_L = 2.0/(rho_L*(rsh+1.0));
  double A_R = 2.0/(rho_R*(rsh+1.0));
  double B_L = p_L * (rsh-1.0) / (rsh+1.0);
  double B_R = p_R * (rsh-1.0) / (rsh+1.0);
  
  double g_L = g_of_p(A_L, B_L, p_pre);
  double g_R = g_of_p(A_R, B_R, p_pre);
  double p_star = (g_L*p_L + g_R*p_R - v_R + v_L) / (g_L + g_R); //~
  
  double v_star = 0.5 * (v_L+v_R + (p_star-p_R)*g_R - (p_star-p_L)*g_L); //~
  double rhoL_star = rho_L *  (p_star/p_L + G6) / (G6*p_star/p_L + 1) ;
  double rhoR_star = rho_R *  (p_star/p_R + G6) / (G6*p_star/p_R + 1) ;
  double aL_star = sqrt(rsh*p_L / rhoL_star);
  double aR_star = sqrt(rsh*p_R / rhoR_star);
  
  //Wavespeeds for Harten decision (First interpretation of aL_star, aR_star was wrong)
  double S0_L = v_L - sos_L; //~
  double S3_R = v_R + sos_R; //~
  double S0_Rstar = v_star - aR_star; //~
  double S3_Lstar = v_star + aL_star; //~
  
  if (S0_L < 0 && 0 < S0_Rstar)
    {
      //Low sonic fix
      //v-a characteristic perturbs low side
      double lam_0 = v_avg - sos_avg; //~

      //wave strength and eigenvector:
      double del_p = p_R - p_L;
      double del_u = u_R - u_L;
      double del_v = v_R - v_L;
      double del_rho = rho_R - rho_L;
      
      double DA_dv = rho_avg * sos_avg * del_v; //~
      double sos_sq = sos_avg * sos_avg;
      double alpha_0 = 0.5/sos_sq * (del_p - DA_dv); //~

      vec4 K_eigen_0 = populate_K0_F(u_avg, v_avg, Vsq_avg, H_avg, sos_avg); //~

      //Get fix for the lam_0 wavespeed, which is currently based on Roe averages
      double lam_fix = S0_L * (S0_Rstar - lam_0) / (S0_Rstar - S0_L);
      
      //get perturbation
      vec4 perts = scavec( lam_fix*alpha_0, K_eigen_0);
      output = vecadd(Fgov(UL) , perts); //~
    }
  
  
  else if (S3_Lstar < 0 && 0 < S3_R )
    {
      //High sonic fix
      //u+a characteristic perturbs high side
      double lam_3 = v_avg + sos_avg; //~

      //wave strength and eigenvector:
      double del_p = p_R - p_L;
      double del_u = u_R - u_L;
      double del_v = v_R - v_L;
      double del_rho = rho_R - rho_L;
      
      double DA_dv = rho_avg * sos_avg * del_v; //~
      double sos_sq = sos_avg * sos_avg;
      double alpha_3 = 0.5/sos_sq * (del_p - DA_dv); //~

      vec4 K_eigen_3 = populate_K3_F(u_avg, v_avg, Vsq_avg, H_avg, sos_avg); //~

      //Get fix for the lam_3 wavespeed
      double lam_fix = S3_R * (lam_3 - S3_Lstar) / (S3_R - S3_Lstar);

      vec4 perts = scavec(lam_fix*alpha_3 , K_eigen_3);
      output = vecdiff(Fgov(UR) , perts); //~
    }
  
  else
    {
      //Regular RoePike procedure
      //Step 2:Define wave speeds in Y direction based on average state----------
      double lam[4]; //wave speeds
      lam[0] = v_avg - sos_avg; //~
      lam[1] = sos_avg;
      lam[2] = sos_avg; 
      lam[3] = v_avg + sos_avg; //~
      
      //Step 3: Populate Right Eigenvectors------------------------
      vec4 K_eigen[4];
      K_eigen[0] = populate_K0_F(u_avg, v_avg, Vsq_avg, H_avg, sos_avg); //~
      K_eigen[1] = populate_K1_F(u_avg, v_avg, Vsq_avg, H_avg, sos_avg); //~
      K_eigen[2] = populate_K2_F(u_avg, v_avg, Vsq_avg, H_avg, sos_avg); //~
      K_eigen[3] = populate_K3_F(u_avg, v_avg, Vsq_avg, H_avg, sos_avg); //~
      
      //Step 4: Get the wave strengths-------------------------------
      double del_p = p_R - p_L;
      double del_u = u_R - u_L;
      double del_v = v_R - v_L;
      double del_rho = rho_R - rho_L;
      double alpha[4]; //wave strengths
      double DA_dv = rho_avg * sos_avg * del_v; //~
      double sos_sq = sos_avg * sos_avg;
      alpha[0] = 0.5/sos_sq * (del_p - DA_dv); //~
      alpha[1] = rho_avg * del_u; //~
      alpha[2] = del_rho - del_p/sos_sq; //~
      alpha[3] = 0.5/sos_sq * (del_p + DA_dv); //~
      
      //Step 5: Calculate wave-based perturbations, add to average flux------
      double scalar[4]; //speed*strength
      for (int i = 0; i < 4; i = i + 1)
	{
	  scalar[i] = alpha[i] * fabs(lam[i]);
	}
      
      vec4 perts = set_vec(0.0); //perturbation, calculated by summation
      for (int i = 0; i < 4; i = i + 1)
	{
	  perts = vecadd( perts, scavec(scalar[i], K_eigen[i]) );
	}
      vec4 flux_avg = scavec( 0.5 , vecadd(Fgov(UL), Fgov(UR)) ); //~
      output = vecdiff( flux_avg, scavec(0.5, perts) );
    }
  
  return output;
}
//IMPORTANTL E_flux_HQ and F_flux_HG choose the Riemann solver.
//I would recommend using same scheme for E, F, fluxes.
vec4 E_flux_HQ(vec4 UL, vec4 UR)
{
  //return E_flux_Riemann(UL, UR);
  //return E_flux_Roe(UL, UR);
  return E_flux_RoeHart(UL, UR);
}
vec4 F_flux_HQ(vec4 UL, vec4 UR)
{
  //return F_flux_Riemann(UL, UR);
  //return F_flux_Roe(UL, UR);
  return F_flux_RoeHart(UL, UR);
}

//A bunch of routines for building
//the initial condition; the IC itself
//is set by U_init.
vec4 explode(double x_loc, double y_loc)
{
  //2D explosion problem
  vec4 output;
  double epiX = 1.0;
  double epiY = 0.5;
  double RAD = 0.25;
  double dist = sqrt( pow(x_loc-epiX,2)+pow(y_loc-epiY,2) );
  double rho = 1.0;
  double p = 1.0;
  double u = 0.0;
  double v = 0.0;
  if (dist > RAD)
    {
      rho = 0.1;
      p = 0.125;
    }
  output = Ugov(rho,u,v,p);
  return output;
}
vec4 explode_mean(double x_loc, double y_loc)
{
  //More violent 2D explosion problem.
  vec4 output;
  double epiX = 1.0;
  double epiY = 0.5;
  double RAD = 0.25;
  double dist = sqrt( pow(x_loc-epiX,2)+pow(y_loc-epiY,2) );
  double rho = 5.0;
  double p = 5.0;
  double u = 0.0;
  double v = 0.0;
  if (dist > RAD)
    {
      rho = 0.1;
      p = 0.125;
    }
  output = Ugov(rho,u,v,p);
  return output;
}
vec4 init_calm(double x_loc, double y_loc) //quiescent IC
{
  vec4 output;
  double rho = rho_free;
  double u = 0;
  double v = 0;
  double p = p_free;
  return Ugov(rho,u,v,p);
}
vec4 Shock_start(double Tau) //shock wave initialization
{
  //currently set to Mach 1.1 attacking wave
  vec4 output;
  double p = 1.245*p_free;
  double T = 1.064983*T_free;
  double rho = p / (Rg*T);
  double v = 0.0;
  double sonic = pow(rsh*p/rho, 0.5);
  double u = 1.1 * a_free - 0.91177*sonic;
  output = Ugov(rho,u,v,p);
  return output;
}
vec4 train_init(double x_loc, double y_loc) //Initializes moving shock at x=shock_loc.
{
  vec4 output;
  output = init_calm(0.0, 0.0); //this function is actually independent of input
  if (x_loc <= shock_loc) //shock location 
    {
      output = Shock_start(0.0);
    }
  return output;
}
vec4 Rault_shockVort(double x_loc, double y_loc) //Shock-Vortex interaction. Very cool.
{
  /*
    Rault's shock-vortex problem from 2003 paper.
    Features standing shock at x=0.5,
    mobile (isentropic?) vortex at x=0.25
   */
  /*
    //My original setup, which runs to time smaller than one
  double rho;
  double u;
  double v;
  double p;
  double p0 = rhoRault_0 * Rg * TRault_0; //pre-shock pressure
  double rho0 = rhoRault_0; //pre-shock density
  double T0 = TRault_0; //pre-shock temperature
  if (x_loc >= 0.5)
    {
      //post-shock conditions
      double ruby = rsh - 1.0;
      double chip = rsh + 1.0;
      double M0sq = Mshock*Mshock;
      double M1 = pow((ruby*M0sq + 2.0) / (2.0*rsh*M0sq - ruby) , 0.5);
      double p1 = p0 * (2.0*rsh*M0sq - ruby) / (chip);
      double rho1 = rho0 * (chip*M0sq) / (ruby*M0sq + 2.0);
      double sos1 = pow(rsh*p1/rho1 , 0.5);
      double u1 = M1 * sos1;
      rho = rho1;
      u = u1;
      v = 0;
      p = p1;
    }
  else
    {
      //pre-shock freestream conditions, plus vortex considerations
      double xsq = pow(x_loc - epiX_Rault , 2);
      double ysq = pow(y_loc - epiY_Rault , 2);
      double rad = pow(xsq + ysq , 0.5);
      double sos0 = pow(rsh*p0 / rho0, 0.5);
      double u0 = Mshock * sos0;
      double vmax = Mvort * sos0; //max azimuthal velocity of vortex
      //Some nasty constants for the temperature profile:
      double C3 = T0;
      double C2 = C3 / (chiRault*pow(rad_a*vmax / (rad_a*rad_a - rad_b*rad_b) , 2) )+ 2.0*rad_b*rad_b*log(rad_b);
      double asq = rad_a*rad_a;
      double bsq = rad_b*rad_b;
      double term1 = -0.5*asq;
      double term2 = pow(rad_a/vmax , 2);
      double term3 = pow(rad_a*vmax / (rad_a*rad_a - rad_b*rad_b) , 2);
      double term4 = 0.5*asq - 2.0*bsq*log(rad_a) - bsq*bsq/(2.0*asq) + C2;
      double C1 = term1 + term2*term3*term4;
      //Start ith freestram conditions, correct if inside vortex
      rho = rho0;
      u = u0;
      v = 0;
      p = p0;
      if (rad <= rad_a) //vortex inner core
	{
	  double Tv = chiRault * pow(vmax/rad_a , 2) * (C1 + 0.5*rad*rad);
	  double vtheta = vmax * rad / rad_a;
	  double theta = atan2(y_loc - epiY_Rault, x_loc - epiX_Rault);
	  double uvort = vtheta * (-sin(theta));
	  double vvort = vtheta * (cos(theta));
	  u = u + uvort;
	  v = v + vvort;
	  p = p0 * pow(Tv/T0 , rsh/(rsh-1.0));
	  rho = rho0 * pow(Tv/T0 , 1/(rsh-1.0));
	}
      else if (rad <= rad_b) //outer vortex core
	{
	  
	  double Tv = chiRault * pow(rad_a*vmax / (asq-bsq) , 2) * (0.5*rad*rad - 2.0*bsq*log(rad) - bsq*bsq/(2.0*rad*rad) + C2);
	  double vtheta = vmax * rad_a / (asq - bsq) * (rad - bsq/rad);
	  double theta = atan2(y_loc - epiY_Rault, x_loc - epiX_Rault);
	  double uvort = vtheta * (-sin(theta));
	  double vvort = vtheta * (cos(theta));
	  u = u + uvort;
	  v = v + vvort;
	  p = p0 * pow(Tv/T0 , rsh/(rsh-1.0));
	  rho = rho0 * pow(Tv/T0 , 1/(rsh-1.0));
	}
      else
	{
	  //do nothing
	}
    }
  */
  //My adjusted setup, runs to t=15 to get vortex near exit
  double sos0 = 1.0/(10.0*Mshock) * (0.25 + 0.75*((rsh+1.0)*Mshock*Mshock) / ((rsh-1.0)*Mshock*Mshock + 2.0)); //pre-shock speed of sound
  double rho0 = 1.0; //pre-shock density
  double p0 = sos0*sos0 * rho0 / rsh; //pre-shock pressure
  double T0 = p0 / (Rg * rho0); //pre-shock temperature
  double rho;
  double u;
  double v;
  double p;

  if (x_loc >= 0.5)
    {
      //post-shock conditions
      double ruby = rsh - 1.0;
      double chip = rsh + 1.0;
      double M0sq = Mshock*Mshock;
      double M1 = pow((ruby*M0sq + 2.0) / (2.0*rsh*M0sq - ruby) , 0.5);
      double p1 = p0 * (2.0*rsh*M0sq - ruby) / (chip);
      double rho1 = rho0 * (chip*M0sq) / (ruby*M0sq + 2.0);
      double sos1 = pow(rsh*p1/rho1 , 0.5);
      double u1 = M1 * sos1;
      rho = rho1;
      u = u1;
      v = 0;
      p = p1;
    }
  else
    {
      //pre-shock freestream conditions, plus vortex considerations
      double xsq = pow(x_loc - epiX_Rault , 2);
      double ysq = pow(y_loc - epiY_Rault , 2);
      double rad = pow(xsq + ysq , 0.5);
      double sos0 = pow(rsh*p0 / rho0, 0.5);
      double u0 = Mshock * sos0;
      double vmax = Mvort * sos0; //max azimuthal velocity of vortex
      //Some nasty constants for the temperature profile:
      double C3 = T0;
      double C2 = C3 / (chiRault*pow(rad_a*vmax / (rad_a*rad_a - rad_b*rad_b) , 2) )+ 2.0*rad_b*rad_b*log(rad_b);
      double asq = rad_a*rad_a;
      double bsq = rad_b*rad_b;
      double term1 = -0.5*asq;
      double term2 = pow(rad_a/vmax , 2);
      double term3 = pow(rad_a*vmax / (rad_a*rad_a - rad_b*rad_b) , 2);
      double term4 = 0.5*asq - 2.0*bsq*log(rad_a) - bsq*bsq/(2.0*asq) + C2;
      double C1 = term1 + term2*term3*term4;
      //Start ith freestram conditions, correct if inside vortex
      rho = rho0;
      u = u0;
      v = 0;
      p = p0;
      if (rad <= rad_a) //vortex inner core
	{
	  double Tv = chiRault * pow(vmax/rad_a , 2) * (C1 + 0.5*rad*rad);
	  double vtheta = vmax * rad / rad_a;
	  double theta = atan2(y_loc - epiY_Rault, x_loc - epiX_Rault);
	  double uvort = vtheta * (-sin(theta));
	  double vvort = vtheta * (cos(theta));
	  u = u + uvort;
	  v = v + vvort;
	  p = p0 * pow(Tv/T0 , rsh/(rsh-1.0));
	  rho = rho0 * pow(Tv/T0 , 1/(rsh-1.0));
	}
      else if (rad <= rad_b) //outer vortex core
	{
	  
	  double Tv = chiRault * pow(rad_a*vmax / (asq-bsq) , 2) * (0.5*rad*rad - 2.0*bsq*log(rad) - bsq*bsq/(2.0*rad*rad) + C2);
	  double vtheta = vmax * rad_a / (asq - bsq) * (rad - bsq/rad);
	  double theta = atan2(y_loc - epiY_Rault, x_loc - epiX_Rault);
	  double uvort = vtheta * (-sin(theta));
	  double vvort = vtheta * (cos(theta));
	  u = u + uvort;
	  v = v + vvort;
	  p = p0 * pow(Tv/T0 , rsh/(rsh-1.0));
	  rho = rho0 * pow(Tv/T0 , 1/(rsh-1.0));
	}
      else
	{
	  //do nothing
	}
    }
  return Ugov(rho,u,v,p);
}
vec4 W_from_U(vec4 U) //transformation from convervatives to primitives
{
  vec4 output;
  double rho = U.ma;
  double u = U.mx/rho;
  double v = U.my/rho;
  double sq_Vmag = u*u + v*v;
  double p = (rsh-1.0) * (U.en - 0.5*rho*sq_Vmag);
  output.ma = rho;
  output.mx = u;
  output.my = v;
  output.en = p;
  return output;
}
vec4 U_from_W(vec4 W) //transformation from primitives to conservatives
{
  vec4 U;
  U.ma = W.ma; //rho
  U.mx = W.ma * W.mx; //rho*u
  U.my = W.ma * W.my; //rho*v
  //U.en = p/(rsh-1.0) + 0.5*rho*(u*u+v*v);
  U.en = W.en/(rsh-1.0) + 0.5*W.ma*(W.mx*W.mx + W.my*W.my); 
  return U;
}
//Set the IC in this subroutine:
vec4 U_init(double x_loc, double y_loc)
{
  return explode_mean(x_loc,y_loc); //use t_final = 0.5
  //return explode(x_loc,y_loc); //Use t_final = 0.5
  //return train_init(x_loc,y_loc);
  //return Rault_shockVort(x_loc,y_loc); //Use t_final = 15
}
vec4 W_init(double x_loc, double y_loc)
{
  return W_from_U(U_init(x_loc,y_loc));
}
vec4 E_reflect(vec4 image) //reflective (no-penetration) slip wall Boundary Procedure
{
  vec4 output;
  output.ma = image.ma;
  output.mx = -image.mx;
  output.my = image.my;
  output.en = image.en;
  return output;
}
vec4 F_reflect(vec4 image) //reflective (no-penetration) slip wall Boundary Procedure
{
  vec4 output;
  output.ma = image.ma;
  output.mx = image.mx;
  output.my = -image.my;
  output.en = image.en;
  return output;
}
vec4 outflow(vec4 image)
{
  return image; //supersonic outflow. Heavyhanded simplification
}


int main()
{
  printf("\n***| Mercury deployed with %d x %d cells |***\n\n", Mx, My);
  
  //The x, y cell coordinates are stored on Global level, outside main.
  x = Allocate_2D_double(Mx,My);
  y = Allocate_2D_double(Mx,My);
  
  //time-dependent output system
  U_out = Allocate_2D_vec4(Mx,My);
  double t_gap = (t_final-t_start) / (10.0);
  t_check[0] = t_start;
  t_check[1] = t_start + t_gap;
  printf("Solution output marks:\n");
  printf("0th check at t=%f, (t%d)\n",t_check[0],0);
  printf("1st check at t=%f, (t%d)\n",t_check[1],1);
  for (int K = 2; K < 11; K = K + 1)
    {
      t_check[K] = t_check[K-1] + t_gap;
      printf("Next check at t=%f (t%d)\n",t_check[K],K); 
    }
  t_target = t_check[0];
  printf("Last Solution check should be t=%f (t%d)\n", t_final, 10);
  
  //Approximate (numerical) Solution. Evolves in time after being set to initial condition.
  D2_vec4 Uh = Allocate_2D_vec4(Mx,My);
  D2_vec4 Uh0 = Allocate_2D_vec4(Mx,My);

  printf("\nCell Geometry is rectangular. Dimensions are:\n");
  printf("delbox_x = %f\n", delbox_x);
  printf("delbox_y = %f\n", delbox_y);
  
  //Quadrature node locations, soon to be initialized:
  double gam[GQdeg];
  double zeta[GQdeg];
  
  double dmin;
  double time = t_start;
  int count = 0;
  double delt;
  double prop_max;

  //***************************************************
  /*
    "Overwritten" variables. These items are calculated
    during every timestep as part of the solution
    evolution process. They are all pieces
    of the commonly applied Muscl-Hancock
    space/time discretization for the Euler equations.
  */
  int BOW;
  int MAST;
  
  //Reconstructed solution gradients. Overwritten during every timestep
  D2_vec4 grad_LR = Allocate_2D_vec4(Mx,My);  
  D2_vec4 grad_BT = Allocate_2D_vec4(Mx,My);  

  //Primitive Limiting Considerations:
  D2_vec4 Wh = Allocate_2D_vec4(Mx,My);  

  vec4 W_local;

  vec4 halfdel_LR;
  vec4 halfdel_BT;
  vec4 U_local;
  vec4 E_L_local;
  vec4 E_R_local;
  vec4 F_B_local;
  vec4 F_T_local;
  vec4 Eterm;
  vec4 Fterm;

  //Predicted interface values
  D2_vec4 UcL = Allocate_2D_vec4(Mx,My); 
  D2_vec4 UcR = Allocate_2D_vec4(Mx,My); 
  D2_vec4 UcB = Allocate_2D_vec4(Mx,My); 
  D2_vec4 UcT = Allocate_2D_vec4(Mx,My); 
  
  //Interface flux values
  D2_vec4 E_flux_L = Allocate_2D_vec4(Mx,My); 
  D2_vec4 E_flux_R = Allocate_2D_vec4(Mx,My); 
  D2_vec4 F_flux_B = Allocate_2D_vec4(Mx,My); 
  D2_vec4 F_flux_T = Allocate_2D_vec4(Mx,My); 

  vec4 change;
  vec4 du_dt;

  //*********************************************************
  //Now, some considerations for setting the initial condition:
  
  //integration holders
  vec4 garage[GQdeg];
  vec4 foyer[GQdeg];

  
  double messenger;
  double garbage;
  
  //*******************************************************
  //output tools:
  vec4 delivery;
  vec4 delivery0;
  double xrec;
  double yrec;

  double rho;
  double rho0;
  double u;
  double v;
  double theta;
  double p_loc;
  double Vel;
  double safety = pow(10,-12);
  double stream;
  double max_Vmag0 = 0;
  double max_u0 = 0;
  double max_v0 = 0;
  double p0;
  double Vel0;

  //*****************************************************
  //Necessary arrays have been built, now
  //deal with some geometry issues and
  //build the initial condition.
  
  dmin = fmin(delbox_x, delbox_y);
  //Length 2 quadrature nodes
  for (int i = 0; i < GQdeg; i = i + 1)
    {
      gam[i] = gwn_return(GQdeg,i).loc;
    }
  //quadrature weights
  for (int i = 0; i < GQdeg; i = i + 1)
    {
      weights[i] = gwn_return(GQdeg,i).weight;
    }
  //Establish Reference Coordinates on [0,1] reference domain
  //These are distributed in prportion to Gaussian quadrature points
  for (int i = 0; i < GQdeg; i = i + 1)
    {
      zeta[i] = trans_up2(gam[i] , 0.0 , 1.0);
    }

  //Establish ce;; centroid coordinates
  for (int I = 0; I < Mx; I = I + 1)
    {
      for (int J = 0; J < My; J = J + 1)
	{
	  x[I][J] = LEFT + delbox_x*(I + 0.5);
	  y[I][J] = BASE + delbox_y*(J + 0.5); //cell centroid coordinates
	}
    }

  //Initialize field variables via Gaussian quadrature for cell average.
  //This is an element-wise operation, done in a way to prevent
  //full storage of the solution at all quadrature points, which
  //would substantially amplify code's memory consumption.
  for (int I = 0; I < Mx; I = I + 1)
    {
      for (int J = 0; J < My; J = J + 1)
	{
	  double xquad;
	  double yquad;
	  D2_vec4 W0_cell = Allocate_2D_vec4(GQdeg,GQdeg);
	  D2_vec4 U0_cell = Allocate_2D_vec4(GQdeg,GQdeg);
	  for (int i = 0; i < GQdeg; i = i + 1)
	    {
	      for (int j = 0; j < GQdeg; j = j + 1)
		{
		  xquad = trans_up1(zeta[i], x[I][J]-0.5*delbox_x, x[I][J]+0.5*delbox_x);
		  yquad = trans_up1(zeta[j], y[I][J]-0.5*delbox_y, y[I][J]+0.5*delbox_y);
		  W0_cell[i][j] = W_init(xquad, yquad);
		  U0_cell[i][j] = U_init(xquad, yquad);
		}
	    }
	  //Combine the solutions at quadrature points to get local cell average.
	  for (int i = 0; i < GQdeg; i = i + 1)
	    {
	      for (int j = 0; j < GQdeg; j = j + 1)
		{
		  garage[j] = W0_cell[i][j];
		}
	      foyer[i] = vec_gquad(garage);
	    }
	  Uh0[I][J] = U_from_W(vec_gquad(foyer));
	  Uh[I][J] = Uh0[I][J]; //This line sets the initial Uh values that will be evolved during time integration.
	  u = Uh0[I][J].mx / Uh0[I][J].ma;
	  v = Uh0[I][J].my / Uh0[I][J].ma;
	  Vel = sqrt(u*u+v*v);
	  max_u0 = fmax(max_u0, u);
	  max_v0 = fmax(max_v0, v);
	  max_Vmag0 = fmax(max_Vmag0, Vel);
	}
    }
 
  if (initialize == 1) //we have an initial values file to use. This is for restarting the code from a stored solution.
    {
      max_Vmag0 = 0;
      max_u0 = 0;
      max_v0 = 0;
      
      FILE*file100;
      file100 = fopen("xU_relay.csv","r");
      printf("Grabbing data from last run\n");
      //trhe file contains primitives rho,u,v,specific energy
      for (int I = 0; I < Mx; I = I + 1)
	{
	  for (int J = 0; J < My; J = J + 1)
	    {
	      //get xrec and yrec entries
	      fscanf(file100,"%lf,", &messenger);
	      xrec = messenger;
	      fscanf(file100,"%lf,", &messenger);
	      yrec = messenger;

	      //Now, get rho,u,v,specific energy
	      fscanf(file100,"%lf,",&messenger);
	      Uh[I][J].ma = messenger;
	      fscanf(file100,"%lf,", &messenger);
	      Uh[I][J].mx = messenger * Uh[I][J].ma;
	      fscanf(file100,"%lf,", &messenger);
	      Uh[I][J].my = messenger *  Uh[I][J].ma;
	      fscanf(file100,"%lf,", &messenger);
	      Uh[I][J].en = messenger * Uh[I][J].ma;
	    }
	}

      fclose(file100);
    
      //Set initial cell averages
      for (int I = 0; I < Mx; I = I + 1)
	{
	  for (int J = 0; J < My; J = J + 1)
	    {
	      Uh0[I][J] = Uh[I][J];
	      u = Uh0[I][J].mx/Uh0[I][J].ma;
	      v = Uh0[I][J].my/Uh0[I][J].ma;
	      Vel = sqrt(u*u+v*v);
	      max_u0 = fmax(max_u0, u);
	      max_v0 = fmax(max_v0, v);
	      max_Vmag0 = fmax(max_Vmag0, Vel);
	    }
	}
    }

  //Have a look at some very general information
  //about the initial condition:
  printf("\nSummary of Initial Condition:\n");
  printf("max axial velocity = %f\n", max_u0);
  printf("max normal velocity = %f\n", max_v0);
  printf("max speed = %f\n", max_Vmag0);

  //Send the t0 solution (initial condition) to file.
  printf("Outputting state for tmark = %d\n", t_mark);
  for (int I = 0; I < Mx; I = I + 1)
    {
      for (int J = 0; J < My; J = J + 1)
	{
	  U_out[I][J] = Uh[I][J];
	}
    }
  write_output(t_mark);
  t_mark = t_mark + 1;
  t_target = t_check[t_mark];
  
  //***************************************************************
  //Initialize timers.

  timer subTimer0("Cumulative Marching");
  timer subTimer1("First Reconstruction");
  timer subTimer2("Second Reconstruction");
  timer subTimer3("Riemann + Update Loop");
  timer subTimer4("Embedded Output");
  
  timers.push_back(&subTimer0);
  timers.push_back(&subTimer1);
  timers.push_back(&subTimer2);
  timers.push_back(&subTimer3);
  timers.push_back(&subTimer4);


  
  //***************************************************************
  //Everything is ready to go. Now, evolve the initial condition
  //forward in time with the Muscl-Hancock discretization.
  double fac_x;
  double fac_y;
  printf("\n---| Setup Complete. Entering explicit time integration loop to t=%f |---\n\n",t_final);
  subTimer0.startBlock();
 
  while (time < t_final)
    { 
      count = count + 1;
      //Need some information to set the timestep size
      //Get maximum wavespeed in physical domain, transform U to W
      prop_max = 0;
      for (int I = 0; I < Mx; I = I + 1)
	{
	  for (int J = 0; J < My; J = J + 1)
	    {
	      prop_max = fmax(prop_max, get_prop(Uh[I][J]));
	      Wh[I][J] = W_from_U(Uh[I][J]);
	    }
	}
      
      delt = nu*dmin / prop_max;
      fac_x = 0.5*delt / delbox_x;
      fac_y = 0.5*delt / delbox_y;
      time = time + delt;     
      if (time > t_final)
	{
	  time = time-delt;
	  delt = t_final-time;
	  printf("Corrected aggressive timestep to finish time loop on the dot\n");
	  time = time+delt;
	}
      if (count % interval_out == 0)
	{
	  printf("current time=%f, steps=%d, prop_max = %f, delt=%f, percent complete= %4.1f\n", time, count, prop_max, delt, time/t_final*100);
	}
      
      //CONSTRUCT GRADIENTS--------------------------------------------------
      subTimer1.startBlock();
      
      //Get dUdx gradients along vertical interfaces:
      for (int J = 0; J < My; J = J + 1) //J means we run through the cell structure in y direction.
	{
	  //boundary cells:
	  grad_LR[0][J] = vec_harmonic(outflow(Uh[0][J]) , Uh[0][J] , Uh[1][J]);
	  grad_LR[Mx-1][J] = vec_harmonic(Uh[Mx-2][J] , Uh[Mx-1][J] , outflow(Uh[Mx-1][J]));
	  
	  //Interior cells:
	  for (int I = 1; I < Mx-1; I = I + 1)
	    {
	      grad_LR[I][J] = vec_harmonic(Wh[I-1][J], Wh[I][J], Wh[I+1][J]);
	    }
	}

      //Get dUdy gradients along horizontal interfaces:
      for (int I = 0; I < Mx; I = I + 1) //I means we run the cell structure in x direction.
	{
	  //boundary cells:

	  //The outflow BC:
       	  //grad_BT[I][0] = vec_harmonic(outflow(Uh[I][0]) , Uh[I][0] , Uh[I][1]);
	  //grad_BT[I][My-1] = vec_harmonic(Uh[I][My-2] , Uh[I][My-1] , outflow(Uh[I][My-1]));

	  //The Slip-wall BC:
	  grad_BT[I][0] = vec_harmonic(F_reflect(Uh[I][0]) , Uh[I][0] , Uh[I][1]);
	  grad_BT[I][My-1] = vec_harmonic(Uh[I][My-2] , Uh[I][My-1] , F_reflect(Uh[I][My-1]));
	  
	  //interior cells:
	  for (int J = 1; J < My-1; J = J + 1)
	    {
	      grad_BT[I][J] = vec_harmonic(Wh[I][J-1], Wh[I][J], Wh[I][J+1]);
	    }
	}
     
      //GET PREDICTED CELL INTERFACE VALUES---------------------------------------------
  
      for (int I = 0; I < Mx; I =I + 1)
	{
	  for (int J = 0; J < My; J = J + 1)
	    {
	      halfdel_LR = scavec(0.5, grad_LR[I][J]);
	      halfdel_BT = scavec(0.5, grad_BT[I][J]); 
	      W_local = Wh[I][J];
	      
	      UcL[I][J] = U_from_W( vecdiff(W_local, halfdel_LR) );
	      UcR[I][J] = U_from_W( vecadd(W_local, halfdel_LR) );
	      UcB[I][J] = U_from_W( vecdiff(W_local, halfdel_BT) );
	      UcT[I][J] = U_from_W( vecadd(W_local, halfdel_BT) );
	    }
	}

      //USE GRADIENTS TO UPDATE TO U PREDICTOR, W PREDICTOR AT HALF OF TIMESTEP------
  
      for (int I = 0; I < Mx; I = I + 1)
	{
	  for (int J = 0; J < My; J = J + 1)
	    {
	      Eterm = vecdiff( Egov(UcL[I][J]), Egov(UcR[I][J]) );
	      Fterm = vecdiff( Fgov(UcB[I][J]), Fgov(UcT[I][J]) );
	      change = vecadd(scavec(fac_x,Eterm), scavec(fac_y,Fterm));
	      Wh[I][J] = W_from_U( vecadd(Uh[I][J], change) );
	    }
	}
      subTimer1.endBlock();
      //USE PREDICTOR W TO RECONSTRUCT CELL GRADIENTS-------------------------------
      subTimer2.startBlock();
      for (int J = 0; J < My; J = J + 1) //vertical interfaces
	{
	  //boundary cells:
	  grad_LR[0][J] = vec_harmonic(outflow(Wh[0][J]) ,Wh[0][J] , Wh[1][J]);
	  grad_LR[Mx-1][J] = vec_harmonic(Wh[Mx-2][J], Wh[Mx-1][J], outflow(Wh[Mx-1][J]));
	  
	  //interior left/right cells
	  for (int I = 1; I < Mx-1; I = I + 1)
	    {
	      grad_LR[I][J] = vec_harmonic(Wh[I-1][J] , Wh[I][J] , Wh[I+1][J]);
	    }
	}
      for (int I = 0; I < Mx; I = I + 1) //horizontal interfaces
	{
	  //boundary cells

	  //outflow BC:
	  //grad_BT[I][0] = vec_harmonic(outflow(Wh[I][0]) , Wh[I][0] , Wh[I][1]);
	  //grad_BT[I][My-1] = vec_harmonic(Wh[I][My-2] , Wh[I][My-1] , outflow(Wh[I][My-1]));

	  //Slip-Wall BC:
	  grad_BT[I][0] = vec_harmonic(F_reflect(Wh[I][0]) , Wh[I][0] , Wh[I][1]);
	  grad_BT[I][My-1] = vec_harmonic(Wh[I][My-2] , Wh[I][My-1] , F_reflect(Wh[I][My-1]));
	  
	  //interior base/top cells
	  for (int J = 1; J < My-1; J = J + 1)
	    {
	      grad_BT[I][J] = vec_harmonic(Wh[I][J-1] , Wh[I][J] , Wh[I][J+1]);
	    }
	}

      //USE RECONSTRUCTED GRADIENTS TO REFINE EDGE PREDICTIONS--------------------
     
      for (int I = 0; I < Mx; I =I + 1)
	{
	  for (int J = 0; J < My; J = J + 1)
	    {
	      halfdel_LR = scavec(0.5, grad_LR[I][J]);
	      halfdel_BT = scavec(0.5, grad_BT[I][J]); 
	      W_local = Wh[I][J];
	      UcL[I][J] = U_from_W( vecdiff(W_local, halfdel_LR) );
	      UcR[I][J] = U_from_W( vecadd(W_local, halfdel_LR) );
	      UcB[I][J] = U_from_W( vecdiff(W_local, halfdel_BT) );
	      UcT[I][J] = U_from_W( vecadd(W_local, halfdel_BT) );
	    }
	}
      subTimer2.endBlock();
      //INVOKE RIEMANN SOLVER TO UPDATE CELL-AVERAGED SOLUTION-----------------
      subTimer3.startBlock();
      //Get F terms
      for (int I = 0; I < Mx; I = I + 1) //horizontal interfaces
	{
	  //boundary interfaces:

	  //Outflow BC:
	  //F_flux_B[I][0] = F_flux_HQ(outflow(UcB[I][0]) , UcB[I][0]);
	  //F_flux_T[I][My-1] = F_flux_HQ(UcT[I][My-1] , outflow(UcT[I][My-1]));

	  //Slip-Wall BC:
	  F_flux_B[I][0] = F_flux_HQ(F_reflect(UcB[I][0]) , UcB[I][0]);
	  F_flux_T[I][My-1] = F_flux_HQ(UcT[I][My-1] , F_reflect(UcT[I][My-1]));
	  
	  //interior base/top interfaces
	  for (int J = 0; J < My-1; J = J + 1)
	    {
	      MAST = (J+1);
	      //top of current cell, base of mast cell
	      F_flux_T[I][J] = F_flux_HQ(UcT[I][J], UcB[I][MAST]);
	      F_flux_B[I][MAST] = F_flux_T[I][J];
	    }
	}
    
      //Get E terms
      for (int J = 0; J < My; J = J + 1) //vertical interfaces
	{
	  //boundary interfaces
	  E_flux_L[0][J] = E_flux_HQ(outflow(UcL[0][J]) , UcL[0][J]);
	  E_flux_R[Mx-1][J] = E_flux_HQ(UcR[Mx-1][J] , outflow(UcR[Mx-1][J]));
	  
	  //interior left/right interfaces
	  for (int I = 0; I < Mx-1; I = I + 1)
	    {
	      BOW = I+1;
	      //right of current cell, left of BOW cell
	      E_flux_R[I][J] = E_flux_HQ(UcR[I][J], UcL[BOW][J]);
	      E_flux_L[BOW][J] = E_flux_R[I][J];
	    }
	}
      
      
      //UPDATE ACTUAL SOLUTION---------------------------------------------------------
      for (int I = 0; I < Mx; I = I + 1)
	{
	  for (int J = 0; J < My; J = J + 1)
	    {
	      du_dt = vecadd( scavec(1.0/delbox_x, vecdiff(E_flux_L[I][J], E_flux_R[I][J])), scavec(1.0/delbox_y, vecdiff(F_flux_B[I][J], F_flux_T[I][J])) );
	      change = scavec(delt, du_dt);
	      Uh[I][J] = vecadd(Uh[I][J], change);
	    }
	}
      subTimer3.endBlock();
      subTimer4.startBlock();
      //IF at predetermined output station, output present distribution
      if (time >= t_target && time != t_final)
	{
	  printf("Outputting state for tmark = %d\n", t_mark);
	  for (int I = 0; I < Mx; I = I + 1)
	    {
	      for (int J = 0; J < My; J = J + 1)
		{
		  U_out[I][J] = Uh[I][J];
		}
	    }
	  write_output(t_mark);
	  t_mark = t_mark + 1;
	  t_target = t_check[t_mark];
	}
      
      subTimer4.endBlock();
      
    } //End master timestepping loop
  subTimer0.endBlock();

  printf("\nOutputting state for tmark = %d, should be end of time loop\n\n", t_mark);
  for (int I = 0; I < Mx; I = I + 1)
    {
      for (int J = 0; J < My; J = J + 1)
	{
	  U_out[I][J] = Uh[I][J];
	}
    }
  write_output(t_mark);
  t_mark = t_mark + 1;
  
  //********************************************************************
  //Simulation Complete.
  //Output:
  FILE*file1;
  file1 = fopen("xU_Final.csv","w");
  for (int I = 0; I < Mx; I = I + 1)
    {
      for (int J = 0; J < My; J = J + 1)
	{
	  rho0 = Uh0[I][J].ma;
	  rho = Uh[I][J].ma;
	  u = Uh[I][J].mx/rho;
	  v = Uh[I][J].my/rho;
	  p_loc = get_pressure(Uh[I][J]);
	  p0 = get_pressure(Uh0[I][J]);
	  Vel = sqrt(u*u + v*v);
	  Vel0 = sqrt(pow(Uh0[I][J].mx/Uh0[I][J].ma,2) + pow(Uh0[I][J].my/Uh0[I][J].ma,2));
	  fprintf(file1,"%8.7f, %8.7f, %16.10f,%16.10f,%16.10f,%16.10f,%16.10f, %16.10f,%16.10f,%16.10f,%16.10f,%16.10f,%16.10f, %16.10f\n", x[I][J], y[I][J], rho0, Uh0[I][J].mx/rho0, Uh0[I][J].my/rho0, Uh0[I][J].en/rho0, p0, Vel0, rho, u, v, Uh[I][J].en/rho, p_loc, Vel);
	}
    }
  fclose(file1);

  FILE*file2;
  file2 = fopen("xU_relay.csv","w");
  for (int I = 0; I < Mx; I = I + 1)
    {
      for (int J = 0; J < My; J = J + 1)
	{
	  rho = Uh[I][J].ma;
	  u = Uh[I][J].mx/rho;
	  v = Uh[I][J].my/rho;
	  fprintf(file2,"%8.7f, %8.7f, %16.10f,%16.10f,%16.10f,%16.10f\n", x[I][J], y[I][J], rho, u, v, Uh[I][J].en/rho);
	}
    }
  fclose(file2);

  FILE*file3;
  file3 = fopen("run_summary.csv","a");
  fprintf(file3,"%d,%d,%8.7f,%8.7f,%f,%f,%f,\n", Mx,My,delbox_x,delbox_y,subTimer0.elapsed, t_start,t_final);
  fclose(file3);

  FILE*file5;
  file5 = fopen("size_params.csv","w");
  fprintf(file5, "%d,%d\n",Mx,My);
  fclose(file5);

  printf("\n---| Timer Information: |---\n");
  for (int t=0; t < timers.size(); t++) 
    {
      timers[t]->print();
    }

  printf("\n***| Mercury execution successful, now exiting. |***\n\n");
  
  return 0;
}
