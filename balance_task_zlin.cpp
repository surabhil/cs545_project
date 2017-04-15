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

// local variables
static double      start_time = 0.0;
static SL_DJstate  target[N_DOFS+1];
static SL_DJstate  target_right[N_DOFS+1];
static SL_DJstate  target_left[N_DOFS+1];
static SL_DJstate  target_base[N_DOFS+1];
static SL_Cstate   cog_target;
static SL_Cstate   cog_traj;
static SL_Cstate   cog_ref;
static double      delta_t = 0.01;
static double      duration = 10.0;
static double      time_to_go;
static int         which_step;
static int         right_foot;

// possible states of a state machine
enum Steps {
  ASSIGN_COG_TARGET,
  MOVE_TO_COG_TARGET,
  ASSIGN_COG_TARGET_BASE,
  MOVE_TO_COG_TARGET_BASE,
  ASSIGN_JOINT_TARGET_LIFT_UP,
  MOVE_JOINT_TARGET_LIFT_UP,
  MOVE_JOINT_TARGET_LEFT_DOWN,
  ASSIGN_JOINT_TARGET_LEFT_DOWN
};

// variables for COG control
static iMatrix     stat;
static Matrix      Jccogp;
static Matrix      NJccog;
static Matrix      fc;
static Vector      delta_xd;
static Vector      delta_thd;

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
    delta_xd  = my_vector(1, N_CART);
    delta_thd = my_vector(1, N_DOFS);

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
  bzero((char *)&(target_right[1]),N_DOFS*sizeof(target_right[1]));
  bzero((char *)&(target_left[1]),N_DOFS*sizeof(target_left[1]));
  bzero((char *)&(target_base[1]),N_DOFS*sizeof(target_base[1]));
  for (i=1; i<=N_DOFS; i++) {
    target[i] = joint_default_state[i];
    target_base[i] = joint_default_state[i];
  }

  right_foot = 1;
  // go to the target using inverse dynamics (ID)
  if (!go_target_wait_ID(target)) 
    return FALSE;

  // ready to go
  freezeBase(FALSE);
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

  // ******************************************
  // NOTE: all array indices start with 1 in SL
  // ******************************************

  task_time = task_servo_time - start_time;

  // the following code computes the contraint COG Jacobian 
  // Jccogp is an N_DOFS x N_CART matrix
  // NJccog is an N_DOFS x N_DOF+2*N_CART matrix -- most likely this is not needed


  // switch according to the current state of the state machine
  switch (which_step) {
    
  case ASSIGN_COG_TARGET:

    // what is the target for the COG?
    bzero((void *)&cog_target,sizeof(cog_target));
#if 0
    cog_target.x[_X_] = cart_des_state[RIGHT_FOOT].x[_X_];
    cog_target.x[_Y_] = cart_des_state[RIGHT_FOOT].x[_Y_];
    cog_target.x[_Z_] = cart_des_state[RIGHT_FOOT].x[_Z_];
# else
    if (right_foot == 1) {
//    cog_target.x[_X_] = (right_foot == 1) ? 0.054 : -0.054;
//    cog_target.x[_Y_] = 0.012;
//    cog_target.x[_Z_] = -0.119;
    cog_target.x[_X_] = cart_des_state[RIGHT_FOOT].x[_X_];
    cog_target.x[_Y_] = cart_des_state[RIGHT_FOOT].x[_Y_];
    cog_target.x[_Z_] = cart_des_state[RIGHT_FOOT].x[_Z_];

    } else {
    cog_target.x[_X_] = cart_des_state[LEFT_FOOT].x[_X_];
    cog_target.x[_Y_] = cart_des_state[LEFT_FOOT].x[_Y_];
    cog_target.x[_Z_] = cart_des_state[LEFT_FOOT].x[_Z_];
   
    }
#endif
    //cog_target.x[_Z_] = 0;

    // the structure cog_des has the current position of the COG computed from the
    // joint_des_state of the robot. cog_des should track cog_traj
    bzero((void *)&cog_traj,sizeof(cog_traj));
    for (i=1; i<=N_CART; ++i)
      cog_traj.x[i] = cog_des.x[i];

    // time to go
    time_to_go = duration;

    // switch to next step of state machine
    which_step = MOVE_TO_COG_TARGET;

    break;

  case MOVE_TO_COG_TARGET: // this is for inverse kinematics control

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
    for (i=1; i<=N_CART; ++i) {
      cog_ref.xd[i] = kp*(cog_traj.x[i] - cog_des.x[i]) + cog_traj.xd[i];
      // create a separate velocity vector for use of the provided mat_vec_mult routine 
      delta_xd[i] = cog_ref.xd[i];
    } 

    // IK solution with constraint COG Jacobian
    compute_cog_kinematics(stat, TRUE, FALSE, TRUE, Jccogp, NJccog);
    mat_vec_mult(Jccogp, delta_xd, delta_thd);
    // compute the joint_des_state[i].th and joint_des_state[i].thd  
    for (i=1; i<=N_DOFS; ++i) {
      // assign theta and theta_d accordingly
      joint_des_state[i].thd  = delta_thd[i];
      joint_des_state[i].th  += delta_thd[i] * delta_t;
      joint_des_state[i].thdd = 0;
      joint_des_state[i].uff  = 0;
    }

    // this is a special inverse dynamics computation for a free standing robot
    inverseDynamicsFloat(delta_t, stat, TRUE, joint_des_state, NULL, NULL, fc);
    // decrement time to go
    time_to_go -= delta_t;
    if (time_to_go <= 0) {
      which_step = ASSIGN_JOINT_TARGET_LIFT_UP;
    }

    break;

  case ASSIGN_COG_TARGET_BASE:

    // what is the target for the COG?
    bzero((void *)&cog_target,sizeof(cog_target));
#if 0
    cog_target.x[_X_] = cart_des_state[RIGHT_FOOT].x[_X_];
    cog_target.x[_Y_] = cart_des_state[RIGHT_FOOT].x[_Y_];
    cog_target.x[_Z_] = cart_des_state[RIGHT_FOOT].x[_Z_];
# else
    cog_target.x[_X_] = 0.5 * (cart_des_state[LEFT_FOOT].x[_X_] + cart_des_state[RIGHT_FOOT].x[_X_]);
    cog_target.x[_Y_] = 0.012;
    cog_target.x[_Z_] = cart_des_state[RIGHT_FOOT].x[_Z_];
#endif
    //cog_target.x[_Z_] = 0;

    // the structure cog_des has the current position of the COG computed from the
    // joint_des_state of the robot. cog_des should track cog_traj
    bzero((void *)&cog_traj,sizeof(cog_traj));
    for (i=1; i<=N_CART; ++i)
      cog_traj.x[i] = cog_des.x[i];

    // time to go
    time_to_go = duration;

    // switch to next step of state machine
    which_step = MOVE_TO_COG_TARGET_BASE;

    break;

  case MOVE_TO_COG_TARGET_BASE: // this is for inverse kinematics control

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
    for (i=1; i<=N_CART; ++i) {
      cog_ref.xd[i] = kp*(cog_traj.x[i] - cog_des.x[i]) + cog_traj.xd[i];
      // create a separate velocity vector for use of the provided mat_vec_mult routine 
      delta_xd[i] = cog_ref.xd[i];
    } 

    // IK solution with constraint COG Jacobian
    compute_cog_kinematics(stat, TRUE, FALSE, TRUE, Jccogp, NJccog);
    mat_vec_mult(Jccogp, delta_xd, delta_thd);
    // compute the joint_des_state[i].th and joint_des_state[i].thd  
    for (i=1; i<=N_DOFS; ++i) {
      // assign theta and theta_d accordingly
      joint_des_state[i].thd  = delta_thd[i];
      joint_des_state[i].th  += delta_thd[i] * delta_t;
      joint_des_state[i].thdd = 0;
      joint_des_state[i].uff  = 0;
    }

    // this is a special inverse dynamics computation for a free standing robot
    inverseDynamicsFloat(delta_t, stat, TRUE, joint_des_state, NULL, NULL, fc);
    // decrement time to go
    time_to_go -= delta_t;
    if (time_to_go <= 0) {
      which_step = ASSIGN_COG_TARGET;
    }

    break;


  case ASSIGN_JOINT_TARGET_LIFT_UP:

    // initialize the target structure from the joint_des_state
    for (i=1; i<=N_DOFS; ++i) {
      target[i] = joint_des_state[i];
      target_right[i] = joint_des_state[i];
      target_left[i] = joint_des_state[i];
    }
#if 1
    if (right_foot == 1) {
      target[L_HAA].th -=  0.4;
      target[L_AAA].th  =  0.5;
      target[R_HAA].th -= 0.25;
      target[R_AAA].th += 0.05;
      right_foot = 0;
    } else {
      target[R_HAA].th -=  0.35;
      target[R_AAA].th  = -0.4;
      target[L_HAA].th -= 0.2;
      target[L_AAA].th -= 0.03;
      right_foot = 1;
    }

#else
    target[L_HFE].th =  0.68;
    target[L_HAA].th = -0.1;
    target[L_KFE].th =  0.95;
    target[L_AFE].th =  0.28;
    target[L_AAA].th =  0.20;
#endif
    time_to_go = duration;  // this may be too fast -- maybe a slower movement is better

    which_step = MOVE_JOINT_TARGET_LIFT_UP;

    break;

  case MOVE_JOINT_TARGET_LIFT_UP:

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

    //SL_InvDynNE(joint_state,joint_des_state,endeff,&base_state,&base_orient);
    // decrement time to go
    time_to_go -= delta_t;

    if (time_to_go <= 0) {
      which_step = ASSIGN_JOINT_TARGET_LEFT_DOWN;
    }

    break;

  case ASSIGN_JOINT_TARGET_LEFT_DOWN:

    // initialize the target structure from the joint_des_state
    for (i=1; i<=N_DOFS; ++i) {
      if (right_foot == 0) {
        target[i] = target_right[i];
      } else {
        target[i] = target_left[i];
      }
    }

    time_to_go = duration;  // this may be too fast -- maybe a slower movement is better

    which_step = MOVE_JOINT_TARGET_LEFT_DOWN;

    break;

  case MOVE_JOINT_TARGET_LEFT_DOWN:

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

    //SL_InvDynNE(joint_state,joint_des_state,endeff,&base_state,&base_orient);
    // decrement time to go
    time_to_go -= delta_t;

    if (time_to_go <= 0) {
      which_step = ASSIGN_COG_TARGET_BASE;
    }

    break;

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
