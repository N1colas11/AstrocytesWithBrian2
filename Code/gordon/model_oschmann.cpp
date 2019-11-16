#include <math.h>
#include <stdlib.h>
#include "typedefs.h"
#include "model_heat.h"
#include "model_oschmann.h"
#include "branch.h"
#include "normal.h"
#include "time_correlation.h"
#include "globals.h"
#include <math.h>

//----------------------------------------------------------------------
ModelOschmann::ModelOschmann() : ModelBase()
{
    printf("Model Oschmann constructor\n");

    // Equation 0: c
    // Equation 1: ce
    n_equ = 4;
    diff.resize(n_equ);
    Globals g;

    if (g.reaction_model == "linear") {
        ;
    } else {
        printf("ModelOschmann::reaction model %s not implemented\n", g.reaction_model.c_str());
        exit(0);
    }

    //-------------- Parameters --------------//
    
    //INCX
    I_NCX_max   = 0.0001; // amp/meter**2
    K_mN         = 87.5; //mole/meter**3
    K_mC         = 1.380; //mole/meter**3
    eta          = 0.35;
    k_sat        = 0.1;


    //-- Li-Rinzel CICR model parameters --//
    gamma=5;

    cai0=100.; 
    cai0=g.cai0; 

    c1=0.185; // ER surface to cyt volume ratio

   
    //Jchan 
    rc=6.;
    d1=0.13;     //IP3 binding affinity
    d5=0.05; 
    d5=0.08234; 

    //Jleak
    rl=0.55;  // comp glio for oschmann model
    rl=0.11; 

    //Jpump
    ver=0.9;
    ver=4.; //Oschmann 2017 incorporates effect of c1 or rho_A
    ver=6.; //ullah et al 2006 incorporates effect of c1 or rho_A
    ver=4.4; // COmp glio for Oschmann 2017  says 44. but this must be a typo

    ker=0.075; // .1 for AM .05 for FM
    ker=0.1; // .1 for AM .05 for FM
    ker=0.05; // .1 for AM .05 for FM

    d2=1.049;     //inactivating Ca2+ binding
    d3=0.9434;     //IP3 binding affinity (with Ca2+ inactivation)
    a2=0.2;

    //-- IP3 Degradation --/
    //IP-5p
    r5p_bar=.86; // CompGlio pg130
    r5p_bar=0.05; // Oschmann2009: .04 for AM .05 for FM 
    //IP3-3K
    O3k = .86; // CompGlio pg130
    O3k = 2.; // Oschmann2009
    kd=0.7; // Oschmann2009 .7 | CompGlio pg130 .5
    k3=1.;

    //-- IP3 Production --/
    
    //-- Poisson IP3-Beta Production --//    
    Obeta_poisson = g.amplitude*g.Obeta;
    half_width = .025; // duration of Glutamate Poisson Pulses || too small for glu input
    printf("half_width = %f\n",half_width);

    //-- Agonist-dependent IP3 Production  --//
    vbeta_bar=0.8; // .2 for AM .5 for FM
    KR=1.3;
    Kp=10.;
    Kpi=0.6;
    glu = 0.0;
    alpha = 0.;


    //-- Endogenous IP3 Production --//
    Odelta = .025; // .02 for AM .05 for FM | COMPGLIO pg130: .025
    kPLCdelta=0.5; //Ca affinity // Oschmann20009: .1 | COMPGLIO pg130: .5
    kdelta=1.; //IP3 inhibiting affinity   // Oschmann2009: 1.5 | COMPGLIO pg130: 1.


    
    //----------------------------------------//

}
//calls exact at time zero, aka sets the initial exact solution..
//----------------------------------------------------------------------//
void ModelOschmann::initialSolution(const Branches& branches)
{
    // another implementation would be to define a function, and call it from 
    // branching model
    Globals g;
    printf("Inside ModelOschmann::initialSolution\n");

    // Random initial conditions

    for (int bb=0; bb < branches.size(); bb++) {
        Branch& b = *branches[bb];
        g.model->setBranch(&b); // set current_branch in model

        b.left_boundary_t0.resize(n_equ);
        b.right_boundary_t0.resize(n_equ);

        RVEC& c  = b.sol[tk0][0]; // equ 0
        RVEC& cer= b.sol[tk0][1]; // equ 1
        RVEC& p  = b.sol[tk0][2]; // equ 1
        RVEC& h  = b.sol[tk0][3]; // equ 2

        for (int i=0; i < b.nx; i++) {
             //init_sol(b.x[i], c[i], p[i], h[i]);
             init_sol(b.x[i], c[i], cer[i], p[i], h[i]);

             printf("initial sol on %s : %f, %f, %f, %f\n",b.name.c_str() ,c[i],cer[i], p[i], h[i]);
        }

        // I should have a way to do this more generally to move it to basea
        // All branches have a right juncture, (which can be a real juncture or the right boundary point)
        RVEC& cj  = b.sol_junc[tk0][0]; // equ 0
        RVEC& cerj = b.sol_junc[tk0][1]; // equ 1
        RVEC& pj  = b.sol_junc[tk0][2]; // equ 1
        RVEC& hj  = b.sol_junc[tk0][3]; // equ 2

        //printf("x[nx-2]= %f\n", b.x[b.nx-2]);
        //printf("x[nx-1]= %f\n", b.x[b.nx-1]);
        //printf("x_junc= %f\n", b.x_junc);
    //
        //init_sol(b.x_junc, cj[0], pj[0], hj[0]);
        init_sol(b.x_junc, cj[0], cerj[0], pj[0], hj[0]);

        //printf("initial sol at juncture: %f, %f, %f, %f\n", cj[0], cej[0], pj[0], rj[0]);
    }
}
//----------------------------------------------------------------------
//----------------------------------------------------------------------
//void ModelOschmann::init_sol(REAL x, REAL& c, REAL& p, REAL& h) 
void ModelOschmann::init_sol(REAL x, REAL& c, REAL& cer, REAL& p, REAL& h) 
{
    // Initialize solution at a single point (x,c,ce)
    Branch& b = *current_branch;
    //printf("rho_A of %s = %f\n", b.name.c_str(),b.rho_A);
    Globals g;
    Normal normal;
    if(g.noise=="None"){
         c =.108; // equilibrium value after 50 s of no input. - randomly initialized at low value to stop initial spiking.
        cer =cai0*(b.rho_A); // equilibrium value after 50 s of no input. - randomly initialized at low value to stop initial spiking.
        h = .6; // equilibrium value after 50 s of no input.
       p = .05; // randomly initialized at low value to stop initial spiking.       
    }
    else{

    REAL scale_I = normal.r4_uniform_01 ( &g.seed );
        REAL scale_C = normal.r4_uniform_01 ( &g.seed );
        REAL scale_H = normal.r4_uniform_01 ( &g.seed );
        c = scale_C*.108;     // equilibrium value after 50 s of no input. - randomly initialized at low value to stop initial spiking.

    printf("Defining inital C_ER with rho_A = %f\n", b.rho_A );

        cer = cai0*(b.rho_A);     // c_cyt  + rho_A * c_ER = c_T c_cit~nM => c_ER = rho_A * c_T
    //          ^    .15 should be b.rho_A however branches are initialized afeter model so we must go back and change. 

        h = scale_H*.6;     // equilibrium value after 50 s of no input.
        p = scale_I*.1;     // randomly initialized at low value to stop initial spiking.
         //c =.108; // equilibrium value after 50 s of no input. - randomly initialized at low value to stop initial spiking.
        //cer =cai0*(b.rho_A); // equilibrium value after 50 s of no input. - randomly initialized at low value to stop initial spiking.
        //h = .6; // equilibrium value after 50 s of no input.
       //p = .05; // randomly initialized at low value to stop initial spiking.       

    }
    return;
}
//
//
//
//----------------------------------------------------------------------
REAL ModelOschmann::Hill(REAL x, REAL cst, REAL m) {
        return pow(x,m) / (pow(x,m) + pow(cst,m));
}


//--------------------------------------------------//
//---------- IP3 Production ------------------------//

REAL ModelOschmann::v_prod(REAL t, int *seed)
{
    Normal normal;
    scale_beta = normal.r4_uniform_01 ( seed );
    //scale_delta = normal.r4_uniform_01 ( seed );

    REAL res;// =     Obeta_poisson * scale_beta * Poisson_Pulses1(t,half_width) 
              //+ Odelta_poisson * scale_delta * Poisson_Pulses2(t,half_width);

    #if 0
    if (Poisson_Pulses1(t,half_width) > 1.){
        res =     Obeta_poisson * scale_beta  ; // DO NOT INCORPORATE STOCH PLC-DELTA
    }
    else{
        es =     Obeta_poisson * scale_beta * Poisson_Pulses1(t,half_width)  ; // DO NOT INCORPORATE STOCH PLC-DELTA
    }
    #endif
    res =     Obeta_poisson * scale_beta * Poisson_Pulses1(t,half_width)  ; // DO NOT INCORPORATE STOCH PLC-DELTA
    return res;
}
// IP3-Beta production with explicit glu val
REAL ModelOschmann::vglu(REAL t, REAL& glu, REAL c, int *seed)
{
    Normal normal;
    REAL scale_glu = normal.r4_uniform_01 ( seed );

    glu = 70.* scale_glu * Poisson_Pulses1(t,half_width);
    //printf("in vglu glu=%f\n", glu);
    
    REAL res = vbeta_bar * Hill( glu , KR * ( 1. + (Kp / KR) * Hill(c , Kpi, 1.) ), .7 );

    return res;
}
// Endogenous IP3 Production
REAL ModelOschmann::vdelta(REAL c, REAL p)
{
    REAL res =  Odelta * (1. - Hill(kdelta, p, 1.) )* Hill(c,kPLCdelta,2.) ; // as defined on pg 122 of Comp. Glio.
    //REAL res =  Odelta * ( Hill(c,kPLCdelta,1.)  + (1-alpha) * Hill(kPLCdelta,c,1.) ); // as defined in pg 164 of Comp. Glio. RESULTS IN DETERMINISTIC OSCILLATIONS without glu input!!
    return res;
}
//--------------------------------------------------//
//---------- IP3 Degradation -----------------------//
REAL ModelOschmann::v5p(REAL p)
{
    REAL res =  r5p_bar * p ;
    return res;
}
REAL ModelOschmann::v3k(REAL c, REAL p)
{
    REAL res =  O3k * Hill(c,kd,4.) *Hill( p,k3,1.) ;
    return res;
}
//--------------------------------------------------//
//------------- h Utility Functions-----------------//
REAL ModelOschmann::q2(REAL p)
{
    REAL res =  d2 * (p + d1) / (p + d3) ;
    return res;
}
REAL ModelOschmann::tau_h(REAL c, REAL p)
{
    REAL res = 1. / ( a2 * ( q2(p) + c  )  ) ;
    return res;
}
REAL ModelOschmann::hinf(REAL c, REAL p)
{
    REAL res = q2(p) / (q2(p) + c) ;
    return res;
}
//--------------------------------------------------//
//-------------- Calcium Flux Utility Functions-----//
REAL ModelOschmann::minf(REAL p)
{
    REAL res =  Hill(p,d1,1) ;
    return res;
}
REAL ModelOschmann::ninf(REAL c)
{
    REAL res =  Hill(c,d5,1) ;
    return res;
}

//--------------------------------------------------//
//----------------- Calcium Fluxes -----------------//
REAL ModelOschmann::Jpump(REAL c)
{
    REAL res =  ver * Hill(c,ker,2.) ;
    return res;
}
//REAL ModelOschmann::Jleak(REAL c)
REAL ModelOschmann::Jleak(REAL c, REAL cer)
{
    //REAL res =  rl * (cai0 - (1. + c1) * c) ;
    REAL res =  rl * (cer - c) ;
    return res;
}
//REAL ModelOschmann::Jchan(REAL c, REAL p, REAL h)
REAL ModelOschmann::Jchan(REAL c, REAL cer, REAL p, REAL h)
{
    //REAL res =  rc * pow(minf(p),3.) * pow(ninf(c),3.) * pow(h,3.) * (cai0 - (1. + c1) * c) ;
    REAL res =  rc * pow(minf(p),3.) * pow(ninf(c),3.) * pow(h,3.) * (cer - c) ;
    return res;
}

//****************************************************************************80
double *ModelOschmann::rk4vec_reaction( double t, int n, double u[] )

////****************************************************************************80
////
////  Purpose:
////
////    RK4VEC_TEST_F evaluates the right hand side of a vector ODE.
////
////  Licensing:
////
////    This code is distributed under the GNU LGPL license. 
////
////  Modified:
////
////    09 October 2013
////
////  Author:
////
////    John Burkardt
////
////  Parameters:
////
////    Input, double T, the current time.
////
////    Input, int N, the dimension of the system.
////
////    Input, double U[N], the current solution value.
////
////    Output, double RK4VEC_TEST_F[N], the value of the derivative, dU/dT.
////
{



    Branch& b = *current_branch;
    Globals g;
      double *uprime;

        uprime = new double[n];

    REAL c_ = u[0];
    REAL cer_ = u[1];
    REAL p_ = u[2];
    REAL h_ = u[3];
    
    //REAL p_ = u[1];
    //REAL h_ = u[2];



    REAL input;
    if(g.noise=="Poisson"){
        input =  vglu(t, glu, c_,&g.seed);
        b.input_history.push_back(b.src_var * v_prod(t,&g.seed));
    }
    else if(g.noise=="None"){
        glu = .6;
        //input = vglu(t, glu, c_,&g.seed);
        //input = v_prod(t,&g.seed);
        REAL a = .1;
        input =  a+ a * sin(t);
    }
    else{
        printf("Noise type not implemented. Input failure in model_oschmann.cpp\n\n");
        exit(0);
    }

    if(g.writing && b.name.c_str() == "soma"){
        printf("Jchan = %f\n", Jchan(c_,cer_,p_,h_));
        printf("Jleak = %f\n", Jleak(c_,cer_)) ;
        printf("Jleak = %f\n", Jpump(c_) );
    }

    REAL ER_flux =  ( Jchan(c_,cer_,p_,h_)  + Jleak(c_,cer_)  - Jpump(c_)); //
    //ER_flux =  (  Jleak(c_,cer_)  - Jpump(c_)); //
      //--------------------------//
    
    b.input_history.push_back( b.src_var * input );
        uprime[0] =ER_flux ; 
        uprime[1] = - (1./b.rho_A )  * ER_flux ; // Comp Glio pg 342
    uprime[2] = b.src_var * v_prod(t,&g.seed) + vdelta(c_,p_) - v3k(c_,p_) - v5p(p_) ;
    uprime[3]  = (hinf(c_,p_) - h_) / tau_h(c_,p_) ;


    b.input_history.push_back( b.src_var * input );
    

    return uprime;
}


//----------------------------------------------------------------------

//----------------------------------------------------------------------//
//--------------------------------------------------//
void ModelOschmann::reactionTerm(const RVVEC& u, RVVEC& reaction_term, REAL t)
{
    //printf("ModelOschmann::inside reactionTerm\n");
    const RVEC& c  = u[0];
    const RVEC& cer  = u[1];
    const RVEC& p  = u[2];
    const RVEC& h  = u[3];

    //const RVEC& p  = u[1];
    //const RVEC& h  = u[2];
    //
    //
    RVEC& react_c  = reaction_term[0];
    RVEC& react_cer  = reaction_term[1];
    RVEC& react_p  = reaction_term[2];
    RVEC& react_h  = reaction_term[3];

    //RVEC& react_p  = reaction_term[1];
    //RVEC& react_h  = reaction_term[2];

    Branch& b = *current_branch;
    int nx = c.size();

    //printf("\nfor %s in reactionTerm  t[tk0] = %f\n",current_branch->name.c_str(),t);

    #if 0
    //printf("nx= %d\n", nx);
    for (int i=0; i < nx; i++) {
        printf("c[%d]= %f\n", i, c[i]);
    }
    #endif

    for (int i=0; i < nx; i++) {
        REAL c_  = c[i]; 
        REAL cer_  = cer[i]; 
        REAL p_  = p[i];
        REAL h_  = h[i];


    Globals g;
    /***************************************** Flux Description *********************************************************************
    * - as defined in De Pitta 2009a - Glutamate Regulation of Calcium and IP3 Oscillatiing and Pulsating Dynamics in Astorcytes    * 
    * Notes:                                                            *
    *    -v5p is a linearized approximation of michaelis-menten dynamics of the degradation.                    * 
    *        *this can only be used if IP3 values are less than 10 micoMole                            * 
    *    -vdelta is Ca2+ dependent IP3 production                                        *
    *    -maximal vglu rate is influenced by Ca2+ (>10 microMole) which is not considered due to Ca2+ constraints        *
        ********************************************************************************************************************************/

    //printf("glu = %f\n", glu);
    //printf("vglu = %f\n", vglu(glu, c_,&g.seed));
    //printf("vdelta1or2 = %f\n", vdelta2(c_,p_));
    
    //react_p[i] = v_prod(t,&g.seed) + vdelta1(c_,p_)  - v3k(c_,p_) - v5p(p_) ;
    REAL input;
    if(g.noise=="Poisson"){
        //input =  v_prod(t,&g.seed);
        input =  vglu(t, glu, c_,&g.seed);
    }
    else if(g.noise=="None"){
        //input = vglu(t, glu, c_,&g.seed);
        //input = v_prod(t,&g.seed);
        REAL a = .1;
        input =  a+ a * sin(t);
    }
    else{
        printf("No input defined in model_depitta.cpp\n");
        exit(0);
    }    
    
    react_p[i] = b.src_var *input + vdelta(c_,p_)  - v3k(c_,p_) - v5p(p_) ;
    if(g.graph=="benchmark_singPt")
    {
        // if benchmarking ensure that input is given.. disregard src_var
        react_p[i] = input + vdelta(c_,p_)  - v3k(c_,p_) - v5p(p_) ;
    }
    //react_p[i] =  vdelta(c_,p_)   - v5p(p_) ;


    if( v_prod(t,&g.seed) > Obeta_poisson ){
        printf("Input to %s = %f\n", b.name.c_str(),  v_prod(t,&g.seed));
        printf("INPUT EXCEEDS MAX VALUE FOR BETA IP3 PRODUCTION");
        exit(0);
    }
    //b.input_history.push_back(input);

    //-------------------------------- ELECTRODIFFUSION NOTE: -------------------------------------//
    //---------------------------------------------------------------------------------------------//
    // NCX is built up by defining currents and then changing to flux densities by incorporating:
            //    F -  Faradays Constant
            //     svr - membrane surface area / tissu volume

    // Steady state definition of ions regulating NCX activity. 
    // dependent on glutamate.
    // pg 50 Oschman Dissertation
    // ***********CHANGE*****************************
    REAL a_Na = -6.72; // mM
    REAL b_Na = -8.24; // mM
    REAL c_Na = 23.14; // mM
    //this does not match up with dissetation figure 2.4 of reduced model section
    a_Na = 9.72; // mM
    b_Na = 10.24; // mM
    c_Na = 15.; // mM


    REAL a_V = -28.;    // mV
    REAL b_V = 11.54;  // mV
    REAL c_V = 48.;    // mV
    //this does not match up with dissetation figure 2.4 of reduced model section
    c_V = -85.;
    b_V = .01;
    a_V = 40.;
    // **********************************************



    REAL Na_o0 = 145.;     // mM
    REAL Na0 = 15.;     // mM
    REAL Ca_o0 = 1800.;     // micoM
    // g is in mM here    
    REAL V = a_V * Hill( glu, b_V,1.) + c_V;
    REAL Na = a_Na * Hill( glu, b_Na,1.) + c_Na;
    REAL Na_o = Na_o0 + Na0 - Na;
    REAL Ca_o =  Ca_o0;
    REAL X_N = 87.5;
    REAL X_C = 1380.;

    REAL kbar = .1;
    kbar = .25; // from Luo Rudy 1994, fits experiments, not used because it leads to Ca2+ uptake.... not a bad thing for us
    REAL V_T = 26.; // thermal voltage Comp Glio pg 340
    REAL v_a = V;
    REAL numer = (pow(Na,3.)/pow(Na_o,3.)) * exp( eta * (v_a/V_T) )  -  (c_/Ca_o) * exp( (eta - 1) * (v_a/V_T) ) ;
    REAL denom = (1 + kbar * exp( (eta-1.)*(v_a/V_T) ) ) ;
    REAL incx0 = .1;

    REAL INCX =  incx0 * Hill( Na_o, X_N,3.) * Hill( Ca_o, X_C,1. )  * ( numer / denom ); 
    #if 0
    printf("glu = %f\n",  glu);
    printf("Na = %f\n", Na);
    printf("V = %f\n", V);
    printf("INCX = %f\n", INCX);
    #endif

    //---------------------------------------------------------------------------------------------//
    //---------------------------------------------------------------------------------------------//
        //react_c[i] = Jchan(c_,p_,h_) + Jleak(c_) - Jpump(c_);

    //printf("Jchan = %f\n", Jchan(c_,cer_,p_,h_));
    //printf("Jleak = %f\n", Jleak(c_,cer_)) ;
    //printf("Jleak = %f\n", Jpump(c_) );

    REAL ER_flux =  ( Jchan(c_,cer_,p_,h_)  + Jleak(c_,cer_)  - Jpump(c_)); //
    //ER_flux =  (  Jleak(c_,cer_)  - Jpump(c_)); //
    //react_c[i] = ER_flux  + ( F / (svr *pow(r_er,.5)) ) * INCX ; // Comp Glio pg 342
    react_c[i] = ER_flux ; 

        react_cer[i] = - (1./b.rho_A )  * ER_flux ; // Comp Glio pg 342
    react_h[i]  = (hinf(c_,p_) - h_) / tau_h(c_,p_) ;

      //--------------------------//

    b.input_history.push_back( b.src_var * input );
    //b.input_history.push_back(( F / (svr *pow(r_er,.5)) ) * INCX );
        // Set reaction term to zero
        #if 0
        react_c[i]  = 0.;
        react_p[i]  = 0.;
        react_h[i]  = 0.;
        #endif
    
    } // end for loop through points - this is essentially useless in compartmental model since we only have 1 pt per branch val
} //end reactionTerm()

