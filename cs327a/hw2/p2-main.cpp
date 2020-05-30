/*======================================================================================
 * hw2-p2-main.cpp
 *
 * Implement operational space motion tracking for a 7-DOF redundant kuka 
 * by completing the "FILL ME IN" section 
 *
 * Elena Galbally, Spring 2019
 *
 *======================================================================================*/

/* --------------------------------------------------------------------------------------
   Include Required Libraries and Files
-----------------------------------------------------------------------------------------*/

#include <iostream>
#include <string>
#include <thread>
#include <math.h>

#include "Sai2Model.h"
#include "Sai2Graphics.h"
#include "Sai2Simulation.h"

#include "timer/LoopTimer.h"

#include <GLFW/glfw3.h> //must be loaded after loading opengl/glew as part of graphicsinterface

using namespace std;

const string world_fname = "resources/hw2/world_2_iiwa.urdf";
const string robot_fname = "../resources/kuka_iiwa/kuka_iiwa.urdf";
const string robot_name = "Kuka-IIWA";
const string ee_link_name = "link6";

const string camera_name = "camera_front";
// const string camera_name = "camera_back_top";
// const string camera_name = "camera_top";

/* --------------------------------------------------------------------------------------
	Simulation Loop Initialization
-------------------------------------------------------------------------------------*/
// simulation loop
bool fSimulationRunning = false;
void control(Sai2Model::Sai2Model* robot, Simulation::Sai2Simulation* sim);
void simulation(Sai2Model::Sai2Model* robot, Simulation::Sai2Simulation* sim);

// initialize window manager
GLFWwindow* glfwInitialize();

// callback to print glfw errors
void glfwError(int error, const char* description);

// callback when a key is pressed
void keySelect(GLFWwindow* window, int key, int scancode, int action, int mods);


/* =======================================================================================
   MAIN LOOP
========================================================================================== */
int main (int argc, char** argv) {
	cout << "Loading URDF world model file: " << world_fname << endl;

	// load graphics scene
	auto graphics = new Sai2Graphics::Sai2Graphics(world_fname, false);

	// load robots
	auto robot = new Sai2Model::Sai2Model(robot_fname, false);

	// load simulation world
	auto sim = new Simulation::Sai2Simulation(world_fname, false);

	// set initial condition
	robot->_q << 125.9/180.0*M_PI,
				39.2/180.0*M_PI,
				-49.2/180.0*M_PI,
				70.0/180.0*M_PI,
				-62.4/180.0*M_PI,
				80.2/180.0*M_PI,
				187.2/180.0*M_PI;
	sim->setJointPositions(robot_name, robot->_q);
	robot->updateModel();
	
	// Eigen::Affine3d ee_trans;
	// robot->transform(ee_trans, ee_link_name);
	// cout << ee_trans.translation().transpose() << endl;
	// cout << ee_trans.rotation() << endl;

	// initialize GLFW window
	GLFWwindow* window = glfwInitialize();

    // set callbacks
	glfwSetKeyCallback(window, keySelect);

	// start the simulation thread first
	fSimulationRunning = true;
	thread sim_thread(simulation, robot, sim);

	// next start the control thread
	thread ctrl_thread(control, robot, sim);
	
    // while window is open:
    while (!glfwWindowShouldClose(window)) {
		// update kinematic models
		// robot->updateModel();

		// update graphics. this automatically waits for the correct amount of time
		int width, height;
		glfwGetFramebufferSize(window, &width, &height);
		graphics->updateGraphics(robot_name, robot);
		graphics->render(camera_name, width, height);
		glfwSwapBuffers(window);
		glFinish();

	    // poll for events
	    glfwPollEvents();
	}

	// stop simulation
	fSimulationRunning = false;
	sim_thread.join();
	ctrl_thread.join();

    // destroy context
    glfwDestroyWindow(window);

    // terminate
    glfwTerminate();

	return 0;
}

/* =======================================================================================
   CONTROL LOOP
========================================================================================== */
void control(Sai2Model::Sai2Model* robot, Simulation::Sai2Simulation* sim) {
	// create a timer
	LoopTimer timer;
	timer.initializeTimer();
	timer.setLoopFrequency(200); //200Hz timer
	double last_time = timer.elapsedTime(); //secs

	// cache variables
	bool fTimerDidSleep = true;
	Eigen::VectorXd tau = Eigen::VectorXd::Zero(robot->dof());

	// additional cache variables
	Eigen::VectorXd g(robot->dof()); //joint space gravity vector
	Eigen::MatrixXd J0(6, robot->dof()); //end effector basic Jacobian
	Eigen::MatrixXd L0(6, 6); //Lambda_0 at end effector
	Eigen::VectorXd p(6); //gravity vector at end-effector
	const Eigen::MatrixXd In = Eigen::MatrixXd::Identity(robot->dof(), robot->dof()); // n x n identity matrix
	Eigen::MatrixXd Jbar(robot->dof(), 6);	//A-weighted generalized inverse of J0
	Eigen::MatrixXd Nbar(robot->dof(), robot->dof()); //I - Jbar*J0, null space projection matrix for Jbar
	Eigen::Vector3d ee_pos; //end effector position
	Eigen::Matrix3d ee_rot_mat; //end effector rotation
	robot->rotation(ee_rot_mat, ee_link_name); // initialize
	Eigen::Quaterniond ee_rot_lambda(ee_rot_mat); // end effector rotation in quaternion form
	Eigen::VectorXd ee_error(6); //end effector operational space instantaneous error
	Eigen::VectorXd v0(6); //end effector velocity
	Eigen::VectorXd v0d(6); //end effector desired velocity
	Eigen::VectorXd dv0d(6); //end effector desired acceleration

	// gains
	double kpj = 50;
	double kvj = 20;
	double kpx = 50;
	double kvx = 20;

	while (fSimulationRunning) { //automatically set to false when simulation is quit
		fTimerDidSleep = timer.waitForNextLoop();

		// update time
		double curr_time = timer.elapsedTime();
		double loop_dt = curr_time - last_time;

		// read joint positions, velocities
		sim->getJointPositions(robot_name, robot->_q);
		sim->getJointVelocities(robot_name, robot->_dq);
		robot->updateModel();
		
		// -------------------------------------------
		// -------------------------------------------
		// FILL ME IN: set joint torques for simulation

    // define desired trajectories, velocities and accelerations
    VectorXd xd(7), d_xd(7), dd_xd(7); 
    double c_temp = cos(2*M_PI*curr_time/5);
    double s_temp = sin(2*M_PI*curr_time/5);
    xd << 0., 
          0.5 + 0.1*cos(2*M_PI*curr_time/5),
          0.65 - 0.05*cos(4*M_PI*curr_time/5), 
          1/sqrt(2)*sin(M_PI/4*c_temp), 
          1/sqrt(2)*cos(M_PI/4*c_temp), 
          1/sqrt(2)*sin(M_PI/4*c_temp), 
          1/sqrt(2)*cos(M_PI/4*c_temp);

    d_xd << 0.,
          -0.1*2*M_PI/5*sin(2*M_PI*curr_time/5),
          0.05*4*M_PI/5*sin(4*M_PI*curr_time/5),
          -M_PI*M_PI/(10*sqrt(2))*cos(M_PI/4*c_temp)*s_temp, 
          M_PI*M_PI/(10*sqrt(2))*sin(M_PI/4*c_temp)*s_temp, 
          -M_PI*M_PI/(10*sqrt(2))*cos(M_PI/4*c_temp)*s_temp, 
          M_PI*M_PI/(10*sqrt(2))*sin(M_PI/4*c_temp)*s_temp;

    dd_xd << 0.,
            -0.1*pow(2*M_PI/5,2)*cos(2*M_PI*curr_time/5),
            0.05*pow(4*M_PI/5,2)*cos(4*M_PI*curr_time/5),
            -(M_PI/4 * 2*M_PI/5)*s_temp*d_xd[4] 
                            + xd[4]*(-M_PI/4 * pow(2*M_PI/5,2) * c_temp),
            (M_PI/4 * 2*M_PI/5)*s_temp*d_xd[3] 
                            + xd[3]*(M_PI/4 * pow(2*M_PI/5,2) * c_temp),
            -(M_PI/4 * 2*M_PI/5)*s_temp*d_xd[4] 
                            + xd[4]*(-M_PI/4 * pow(2*M_PI/5,2) * c_temp),
            (M_PI/4 * 2*M_PI/5)*s_temp*d_xd[3] 
                            + xd[3]*(M_PI/4 * pow(2*M_PI/5,2) * c_temp);

    // get rotation matrix
    Matrix3d rot;
    robot->rotation(rot, ee_link_name);

    Eigen::Quaterniond q(rot);

    double l1 = q.w();
    double l2 = q.x();
    double l3 = q.y();
    double l4 = q.z();

    // get representation matrix
    MatrixXd Er_inv(3,4);
    Er_inv << -l2, l1, -l4, l3,
              -l3, l4, l1, -l2,
              -l4, -l3, l2, l1;
    Er_inv = 2.*Er_inv;

    // get ee linear position
    Eigen::Vector3d ee_x;
    robot->position(ee_x, ee_link_name);
    // get ee linear velocity
    Eigen::Vector3d ee_dx;
    robot->linearVelocity(ee_dx, ee_link_name);
    // get ee angular velocity
    Eigen::Vector3d ee_dxr;
    robot->angularVelocity(ee_dxr, ee_link_name);
    // get delta phi
    Eigen::Vector3d delta_phi;
    delta_phi = -Er_inv * xd.tail(4);
    // get desired angular acceleration
    Eigen::Vector3d wd;
    wd = Er_inv * d_xd.tail(4);
    // get desired angular acceleration
    Eigen::Vector3d dwd;
    dwd = 2 * Er_inv * dd_xd.tail(4);
    
    // calculate control force and moment
    VectorXd f_star(3), m_star(3), f(6);
    f_star = dd_xd.head(3) - kpx * (ee_x - xd.head(3)) - kvx * (ee_dx - d_xd.head(3));
    m_star = dwd - kpx * delta_phi - kvx * (ee_dxr - wd);
    f.head(3) = f_star; f.tail(3) = m_star;   // put together wrench

    // get Jacobian and inverse jacobian
    robot->J_0(J0, ee_link_name);
    robot->operationalSpaceMatrices(L0, Jbar, Nbar, J0);
    Eigen::MatrixXd A(robot->dof(), robot->dof());
    A = robot->_M;

    // get null space damping
    Eigen::VectorXd ddq(robot->dof());
    ddq = (- kvj * robot->_dq);

    // get p
    robot->gravityVector(g);
    p = L0*J0*A.inverse()*g;

    // calculate tau
    tau = J0.transpose()*(L0 * f + p) + Nbar * (A * ddq + g);

		sim->setJointTorques(robot_name, tau);
		// -------------------------------------------
		// -------------------------------------------

		// update last time
		last_time = curr_time;
	}
}


/* =======================================================================================
   SIMULATION SETUP
   -----------------------
   * Simulation loop
   * Window initialization
   * Window error
   * Mouse click commands
========================================================================================== */
void simulation(Sai2Model::Sai2Model* robot, Simulation::Sai2Simulation* sim) {
	fSimulationRunning = true;

	// create a timer
	LoopTimer timer;
	timer.initializeTimer();
	timer.setLoopFrequency(1000); //1000Hz timer
	double last_time = timer.elapsedTime(); //secs

	bool fTimerDidSleep = true;
	while (fSimulationRunning) {
		fTimerDidSleep = timer.waitForNextLoop();

		// integrate forward
		double curr_time = timer.elapsedTime();
		double loop_dt = curr_time - last_time; 
		sim->integrate(loop_dt);

		// if (!fTimerDidSleep) {
		// 	cout << "Warning: timer underflow! dt: " << loop_dt << "\n";
		// }

		// update last time
		last_time = curr_time;
	}
}

//------------------------------------------------------------------------------
GLFWwindow* glfwInitialize() {
	
	/*------- Set up visualization -------*/
    // set up error callback
    glfwSetErrorCallback(glfwError);

    // initialize GLFW
    glfwInit();

    // retrieve resolution of computer display and position window accordingly
    GLFWmonitor* primary = glfwGetPrimaryMonitor();
    const GLFWvidmode* mode = glfwGetVideoMode(primary);

    // information about computer screen and GLUT display window
	int screenW = mode->width;
    int screenH = mode->height;
    int windowW = 0.8 * screenH;
    int windowH = 0.5 * screenH;
    int windowPosY = (screenH - windowH) / 2;
    int windowPosX = windowPosY;

    // create window and make it current
    glfwWindowHint(GLFW_VISIBLE, 0);
    GLFWwindow* window = glfwCreateWindow(windowW, windowH, "SAI2.0 - CS327a HW2", NULL, NULL);
	glfwSetWindowPos(window, windowPosX, windowPosY);
	glfwShowWindow(window);
    glfwMakeContextCurrent(window);
	glfwSwapInterval(1);

	return window;
}

//------------------------------------------------------------------------------

void glfwError(int error, const char* description) {
	cerr << "GLFW Error: " << description << endl;
	exit(1);
}

//------------------------------------------------------------------------------

void keySelect(GLFWwindow* window, int key, int scancode, int action, int mods)
{
    // option ESC: exit
    if (key == GLFW_KEY_ESCAPE && action == GLFW_PRESS)
    {
        // exit application
         glfwSetWindowShouldClose(window, 1);
    }
}
