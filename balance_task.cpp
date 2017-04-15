/*============================================================================
==============================================================================
                      
                              balance_task.c
 
==============================================================================
Remarks:

      sekeleton to create the sample task

============================================================================*/

// system headers
#include "SL_system_headers.h"

/* SL includes */
#include "SL.h"
#include "SL_user.h"
#include "SL_tasks.h"
#include "SL_task_servo.h"
#include "SL_kinematics.h"
#include "SL_dynamics.h"
#include "SL_collect_data.h"
#include "SL_shared_memory.h"
#include "SL_man.h"

// defines
// #define VERBOSE_LOG
#define LOG

// local variables
static double      start_time = 0.0;
static SL_DJstate  target[N_DOFS+1];
static SL_DJstate  saved_target[N_DOFS+1];
static SL_Cstate   cog_target;
static SL_Cstate   cog_traj;
static SL_Cstate   cog_ref;
static double      delta_t = 0.01;
static double      duration = 10.0;
static double      time_to_go;
static int         which_step;
static int         init_balance_foot = RIGHT_FOOT;
static int         balance_foot;
static int         iter;

// possible states of a state machine
enum Steps {
  ASSIGN_COG_TARGET,
  MOVE_TO_COG_TARGET,
  ASSIGN_JOINT_TARGET_LIFT_UP,
  MOVE_JOINT_TARGET_LIFT_UP,
  ASSIGN_JOINT_TARGET_LOWER_DOWN,
  MOVE_JOINT_TARGET_LOWER_DOWN,
  ASSIGN_RETURN_CENTER,
  MOVE_RETURN_CENTER
};

// variables for COG control
static iMatrix     stat;
static Matrix      Jccogp;
static Matrix      NJccog;
static Matrix      fc;

// global functions 
extern "C" void
add_balance_task( void );

// local functions
static int  init_balance_task(void);
static int  run_balance_task(void);
static int  change_balance_task(void);

static int 
min_jerk_next_step (double x,double xd, double xdd, double t, double td, double tdd,
		    double t_togo, double dt,
		    double *x_next, double *xd_next, double *xdd_next);


/*****************************************************************************
******************************************************************************
Function Name	: add_balance_task
Date		: Feb 1999
Remarks:

adds the task to the task menu

******************************************************************************
Paramters:  (i/o = input/output)

none

*****************************************************************************/
void
add_balance_task( void )
{
  int i, j;
  
  addTask("Balance Task", init_balance_task, 
	  run_balance_task, change_balance_task);

}    

/*****************************************************************************
******************************************************************************
  Function Name	: init_balance_task
  Date		: Dec. 1997

  Remarks:

  initialization for task

******************************************************************************
  Paramters:  (i/o = input/output)

       none

 *****************************************************************************/
static int 
init_balance_task(void)
{
  int j, i;
  int ans;
  static int firsttime = TRUE;
  
  if (firsttime){
    firsttime = FALSE;

    // allocate memory
    stat   = my_imatrix(1,N_ENDEFFS,1,2*N_CART);
    Jccogp = my_matrix(1,N_DOFS,1,N_CART);
    NJccog = my_matrix(1,N_DOFS,1,N_DOFS+2*N_CART);
    fc     = my_matrix(1,N_ENDEFFS,1,2*N_CART);

    // this is an indicator which Cartesian components of the endeffectors are constraints
    // i.e., both feet are on the ground and cannot move in position or orientation
    stat[RIGHT_FOOT][1] = TRUE;
    stat[RIGHT_FOOT][2] = TRUE;
    stat[RIGHT_FOOT][3] = TRUE;
    stat[RIGHT_FOOT][4] = TRUE;
    stat[RIGHT_FOOT][5] = TRUE;
    stat[RIGHT_FOOT][6] = TRUE;

    stat[LEFT_FOOT][1] = TRUE;
    stat[LEFT_FOOT][2] = TRUE;
    stat[LEFT_FOOT][3] = TRUE;
    stat[LEFT_FOOT][4] = TRUE;
    stat[LEFT_FOOT][5] = TRUE;
    stat[LEFT_FOOT][6] = TRUE;

  }

  // prepare going to the default posture
  bzero((char *)&(target[1]),N_DOFS*sizeof(target[1]));
  for (i=1; i<=N_DOFS; i++)
    target[i] = joint_default_state[i];

  // TODO: set better initial position before shifting cog?

  // go to the target using inverse dynamics (ID)
  if (!go_target_wait_ID(target)) 
    return FALSE;

  // ready to go
  ans = 999;
  while (ans == 999) {
    if (!get_int("Enter 1 to start or anthing else to abort ...",ans,&ans))
      return FALSE;
  }
  
  // only go when user really types the right thing
  if (ans != 1) 
    return FALSE;

  start_time = task_servo_time;
  printf("start time = %.3f, task_servo_time = %.3f\n", 
	 start_time, task_servo_time);

  // start data collection
  scd();

  // state machine starts at ASSIGN_COG_TARGET
  which_step = ASSIGN_COG_TARGET;
  balance_foot = init_balance_foot;

  iter = 0;

  return TRUE;
}

/*****************************************************************************
******************************************************************************
  Function Name	: run_balance_task
  Date		: Dec. 1997

  Remarks:

  run the task from the task servo: REAL TIME requirements!

******************************************************************************
  Paramters:  (i/o = input/output)

  none

 *****************************************************************************/
static int 
run_balance_task(void)
{
  int j, i, n;

  double task_time;
  double kp = 0.1;

  double duration_scale = 1.0;

  double forward_offset = 0.00;
  // double forward_offset = 0.02;

  // ******************************************
  // NOTE: all array indices start with 1 in SL
  // ******************************************

  task_time = task_servo_time - start_time;

  // the following code computes the contraint COG Jacobian 
  // Jccogp is an N_DOFS x N_CART matrix
  // NJccog is an N_DOFS x N_DOF+2*N_CART matrix -- most likely this is not needed

  compute_cog_kinematics(stat, TRUE, FALSE, TRUE, Jccogp, NJccog);

  // switch according to the current state of the state machine
  switch (which_step) {
    
  case ASSIGN_COG_TARGET:

#if defined VERBOSE_LOG or defined LOG
    printf("assign cog target iter %d\n", iter);
#endif

    // what is the target for the COG?
    bzero((void *)&cog_target,sizeof(cog_target));

    // where_cog for target when on right foot:
    // 0.054 ( 0.055)   y= 0.012 ( 0.014)   z=-0.119 (-0.117)
    cog_target.x[_X_] =  (RIGHT_FOOT == balance_foot) ? base_state.x[_X_] + 0.054 : base_state.x[_X_] + -0.054;
    // cog_target.x[_X_] =  (RIGHT_FOOT == balance_foot) ? 0.04 : -0.04;
    cog_target.x[_Y_] =  base_state.x[_Y_] + 0.012 + forward_offset;
    cog_target.x[_Z_] =  base_state.x[_Z_] - 0.119;

    // the structure cog_des has the current position of the COG computed from the
    // joint_des_state of the robot. cog_des should track cog_traj
    bzero((void *)&cog_traj,sizeof(cog_traj));
    for (i=1; i<=N_CART; ++i)
      cog_traj.x[i] = cog_des.x[i];

    // time to go
    time_to_go = duration;

    // switch to next step of state machine
    which_step = MOVE_TO_COG_TARGET;

#ifdef LOG
    printf("move cog target iter %d\n", iter);
#endif

    break;

  case MOVE_TO_COG_TARGET: // this is for inverse kinematics control

#ifdef VERBOSE_LOG
    printf("move cog target iter %d\n", iter);
#endif

    // plan the next step of cog with min jerk
    for (i=1; i<=N_CART; ++i) {
      min_jerk_next_step(cog_traj.x[i],
			 cog_traj.xd[i],
			 cog_traj.xdd[i],
			 cog_target.x[i],
			 cog_target.xd[i],
			 cog_target.xdd[i],
			 time_to_go,
			 delta_t,
			 &(cog_traj.x[i]),
			 &(cog_traj.xd[i]),
			 &(cog_traj.xdd[i]));
    }

    // inverse kinematics: we use a P controller to correct for tracking erros
    for (i=1; i<=N_CART; ++i)
      cog_ref.xd[i] = kp*(cog_traj.x[i] - cog_des.x[i]) + cog_traj.xd[i];


    // compute the joint_des_state[i].th and joint_des_state[i].thd  
    for (i=1; i<=N_DOFS; ++i) {

      // intialize to zero
      joint_des_state[i].thd  = 0;
      joint_des_state[i].thdd = 0;
      joint_des_state[i].uff  = 0;

      // joint_des_state[i].thd = ?;
      // joint_des_state[i].thd = ?;

    }

    for (int rowIdx = 1; rowIdx <= N_DOFS; rowIdx++)
    {
        for (int colIdx = 1; colIdx <= N_CART; colIdx++)
        {
            joint_des_state[rowIdx].thd += Jccogp[rowIdx][colIdx] * cog_ref.xd[colIdx];
        }
    }

    for (int rowIdx = 1; rowIdx <= N_DOFS; rowIdx++)
    {
        joint_des_state[rowIdx].th += joint_des_state[rowIdx].thd * delta_t;
    }

    // decrement time to go
    time_to_go -= delta_t;
    if (time_to_go <= 0) {
      which_step = ASSIGN_JOINT_TARGET_LIFT_UP;

      // save current state to return leg to ground
      for (i=1; i<=N_DOFS; ++i)
          saved_target[i] = joint_des_state[i];
    }

    break;

  case ASSIGN_JOINT_TARGET_LIFT_UP:

#if defined VERBOSE_LOG or defined LOG
    printf("assign joint target lift up iter %d\n", iter);
#endif

    // initialize the target structure from the joint_des_state
    for (i=1; i<=N_DOFS; ++i)
      target[i] = joint_des_state[i];

    if (RIGHT_FOOT == balance_foot)
    {
        // target[L_HAA].th -=  0.4;
        // target[L_AAA].th  =  0.5;

        // target[R_HAA].th -=  0.25;
        // target[R_AAA].th +=  0.05;

        target[L_HFE].th += 1.0;
        target[L_KFE].th += 0.5;

        // stat[LEFT_FOOT][1] = FALSE;
        // stat[LEFT_FOOT][2] = FALSE;
        // stat[LEFT_FOOT][3] = FALSE;
        // stat[LEFT_FOOT][4] = FALSE;
        // stat[LEFT_FOOT][5] = FALSE;
        // stat[LEFT_FOOT][6] = FALSE;

        stat[RIGHT_FOOT][1] = TRUE;
        stat[RIGHT_FOOT][2] = TRUE;
        stat[RIGHT_FOOT][3] = TRUE;
        stat[RIGHT_FOOT][4] = TRUE;
        stat[RIGHT_FOOT][5] = TRUE;
        stat[RIGHT_FOOT][6] = TRUE;
    }
    else
    {
        // target[R_HAA].th -=  0.4;
        // target[R_AAA].th  = -0.5;

        // target[L_HAA].th -=  0.25;
        // target[L_AAA].th -=  0.05;

        target[R_HFE].th += 1.0;
        target[R_KFE].th += 0.5;

        // stat[RIGHT_FOOT][1] = FALSE;
        // stat[RIGHT_FOOT][2] = FALSE;
        // stat[RIGHT_FOOT][3] = FALSE;
        // stat[RIGHT_FOOT][4] = FALSE;
        // stat[RIGHT_FOOT][5] = FALSE;
        // stat[RIGHT_FOOT][6] = FALSE;

        stat[LEFT_FOOT][1] = TRUE;
        stat[LEFT_FOOT][2] = TRUE;
        stat[LEFT_FOOT][3] = TRUE;
        stat[LEFT_FOOT][4] = TRUE;
        stat[LEFT_FOOT][5] = TRUE;
        stat[LEFT_FOOT][6] = TRUE;
    }

    duration_scale = 1.0;
    // duration_scale = 5.0;

    time_to_go = duration_scale * duration;  // this may be too fast -- maybe a slower movement is better

    which_step = MOVE_JOINT_TARGET_LIFT_UP;

#ifdef LOG
    printf("move joint target lift up iter %d\n", iter);
#endif

    break;

  case MOVE_JOINT_TARGET_LIFT_UP:

#ifdef VERBOSE_LOG
    printf("move joint target lift up iter %d\n", iter);
#endif

    // compute the update for the desired states
    for (i=1; i<=N_DOFS; ++i) {
      min_jerk_next_step(joint_des_state[i].th,
			 joint_des_state[i].thd,
			 joint_des_state[i].thdd,
			 target[i].th,
			 target[i].thd,
			 target[i].thdd,
			 time_to_go,
			 delta_t,
			 &(joint_des_state[i].th),
			 &(joint_des_state[i].thd),
			 &(joint_des_state[i].thdd));
    }

    // decrement time to go
    time_to_go -= delta_t;

    if (time_to_go <= 0)
    {
        // freeze();
        which_step = ASSIGN_JOINT_TARGET_LOWER_DOWN;
    }

    break;

  case ASSIGN_JOINT_TARGET_LOWER_DOWN:

#if defined VERBOSE_LOG or defined LOG
    printf("assign joint target lower down iter %d\n", iter);
#endif

    // initialize the target structure from the saved target with both legs
    // on the ground
    for (i=1; i<=N_DOFS; ++i)
        target[i] = saved_target[i];

    duration_scale = 1.0;
    // duration_scale = 5.0;

    time_to_go = duration_scale * duration;  // this may be too fast -- maybe a slower movement is better

    which_step = MOVE_JOINT_TARGET_LOWER_DOWN;

#ifdef LOG
    printf("move joint target lower down iter %d\n", iter);
#endif

    break;

  case MOVE_JOINT_TARGET_LOWER_DOWN:

#ifdef VERBOSE_LOG
    printf("move joint target lower down iter %d\n", iter);
#endif

    // compute the update for the desired states
    for (i=1; i<=N_DOFS; ++i) {
      min_jerk_next_step(joint_des_state[i].th,
			 joint_des_state[i].thd,
			 joint_des_state[i].thdd,
			 target[i].th,
			 target[i].thd,
			 target[i].thdd,
			 time_to_go,
			 delta_t,
			 &(joint_des_state[i].th),
			 &(joint_des_state[i].thd),
			 &(joint_des_state[i].thdd));
    }

    // decrement time to go
    time_to_go -= delta_t;

    if (time_to_go <= 0)
    {
      // freeze();
      which_step = ASSIGN_RETURN_CENTER;
      balance_foot = (RIGHT_FOOT == balance_foot) ? LEFT_FOOT : RIGHT_FOOT;

      stat[RIGHT_FOOT][1] = TRUE;
      stat[RIGHT_FOOT][2] = TRUE;
      stat[RIGHT_FOOT][3] = TRUE;
      stat[RIGHT_FOOT][4] = TRUE;
      stat[RIGHT_FOOT][5] = TRUE;
      stat[RIGHT_FOOT][6] = TRUE;

      stat[LEFT_FOOT][1] = TRUE;
      stat[LEFT_FOOT][2] = TRUE;
      stat[LEFT_FOOT][3] = TRUE;
      stat[LEFT_FOOT][4] = TRUE;
      stat[LEFT_FOOT][5] = TRUE;
      stat[LEFT_FOOT][6] = TRUE;
    }

    break;

  case ASSIGN_RETURN_CENTER:

#if defined VERBOSE_LOG or defined LOG
    printf("assign return center iter %d\n", iter);
#endif

    // prepare going to the default posture
    bzero((char *)&(target[1]),N_DOFS*sizeof(target[1]));
    for (i=1; i<=N_DOFS; i++)
        target[i] = joint_default_state[i];

    duration_scale = 1.0;
    // duration_scale = 5.0;

    time_to_go = duration_scale * duration;  // this may be too fast -- maybe a slower movement is better
    which_step = MOVE_RETURN_CENTER;

#ifdef LOG
    printf("move return center iter %d\n", iter);
#endif

    break;

  case MOVE_RETURN_CENTER:

#ifdef VERBOSE_LOG
    printf("move return center iter %d\n", iter);
#endif

    // compute the update for the desired states
    for (i=1; i<=N_DOFS; ++i) {
      min_jerk_next_step(joint_des_state[i].th,
			 joint_des_state[i].thd,
			 joint_des_state[i].thdd,
			 target[i].th,
			 target[i].thd,
			 target[i].thdd,
			 time_to_go,
			 delta_t,
			 &(joint_des_state[i].th),
			 &(joint_des_state[i].thd),
			 &(joint_des_state[i].thdd));
    }

    // decrement time to go
    time_to_go -= delta_t;

    if (time_to_go <= 0)
    {
        which_step = ASSIGN_COG_TARGET;

        // increment counter if we have gone back to initial balance foot
        if (init_balance_foot == balance_foot)
        {
            iter++;
        }
    }

    break;
  }

  // if ((MOVE_JOINT_TARGET_LIFT_UP != which_step) && (MOVE_JOINT_TARGET_LOWER_DOWN != which_step))
  // if ((MOVE_JOINT_TARGET_LIFT_UP > which_step) || (MOVE_RETURN_CENTER == which_step))
  if ((MOVE_JOINT_TARGET_LIFT_UP > which_step))
  {
      // this is a special inverse dynamics computation for a free standing robot
      inverseDynamicsFloat(delta_t, stat, TRUE, joint_des_state, NULL, NULL, fc);
  }

  return TRUE;
}

/*****************************************************************************
******************************************************************************
  Function Name	: change_balance_task
  Date		: Dec. 1997

  Remarks:

  changes the task parameters

******************************************************************************
  Paramters:  (i/o = input/output)

  none

 *****************************************************************************/
static int 
change_balance_task(void)
{
  int    ivar;
  double dvar;

  get_int("This is how to enter an integer variable",ivar,&ivar);
  get_double("This is how to enter a double variable",dvar,&dvar);

  return TRUE;

}


/*!*****************************************************************************
 *******************************************************************************
\note  min_jerk_next_step
\date  April 2014
   
\remarks 

Given the time to go, the current state is updated to the next state
using min jerk splines

 *******************************************************************************
 Function Parameters: [in]=input,[out]=output

 \param[in]          x,xd,xdd : the current state, vel, acceleration
 \param[in]          t,td,tdd : the target state, vel, acceleration
 \param[in]          t_togo   : time to go until target is reached
 \param[in]          dt       : time increment
 \param[in]          x_next,xd_next,xdd_next : the next state after dt

 ******************************************************************************/
static int 
min_jerk_next_step (double x,double xd, double xdd, double t, double td, double tdd,
		    double t_togo, double dt,
		    double *x_next, double *xd_next, double *xdd_next)

{
  double t1,t2,t3,t4,t5;
  double tau,tau1,tau2,tau3,tau4,tau5;
  int    i,j;

  // a safety check
  if (dt > t_togo || dt <= 0) {
    return FALSE;
  }

  t1 = dt;
  t2 = t1 * dt;
  t3 = t2 * dt;
  t4 = t3 * dt;
  t5 = t4 * dt;

  tau = tau1 = t_togo;
  tau2 = tau1 * tau;
  tau3 = tau2 * tau;
  tau4 = tau3 * tau;
  tau5 = tau4 * tau;

  // calculate the constants
  const double dist   = t - x;
  const double p1     = t;
  const double p0     = x;
  const double a1t2   = tdd;
  const double a0t2   = xdd;
  const double v1t1   = td;
  const double v0t1   = xd;
  
  const double c1 = 6.*dist/tau5 + (a1t2 - a0t2)/(2.*tau3) - 
    3.*(v0t1 + v1t1)/tau4;
  const double c2 = -15.*dist/tau4 + (3.*a0t2 - 2.*a1t2)/(2.*tau2) +
    (8.*v0t1 + 7.*v1t1)/tau3; 
  const double c3 = 10.*dist/tau3+ (a1t2 - 3.*a0t2)/(2.*tau) -
    (6.*v0t1 + 4.*v1t1)/tau2; 
  const double c4 = xdd/2.;
  const double c5 = xd;
  const double c6 = x;
  
  *x_next   = c1*t5 + c2*t4 + c3*t3 + c4*t2 + c5*t1 + c6;
  *xd_next  = 5.*c1*t4 + 4*c2*t3 + 3*c3*t2 + 2*c4*t1 + c5;
  *xdd_next = 20.*c1*t3 + 12.*c2*t2 + 6.*c3*t1 + 2.*c4;
  
  return TRUE;
}
