/*======================================================================================
 * hw1-main.cpp
 *
 * Implement an algorithm to have a 7-DOF kuka track a desired trajectory
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

#include <GLFW/glfw3.h> 

using namespace std;

const string world_fname = "resources/hw1/world.urdf";
const string robot_fname = "../resources/kuka_iiwa/kuka_iiwa.urdf";
const string robot_name = "Kuka-IIWA";
const string camera_name = "camera_front";
// const string camera_name = "camera_top";
const string ee_link_name = "link6";

/* --------------------------------------------------------------------------------------
	Simulation Loop Initialization
-------------------------------------------------------------------------------------*/

// simulation loop
bool fSimulationRunning = false;
void simulation(Sai2Model::Sai2Model* robot);

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

	// set initial condition
	robot->_q << 125.9/180.0*M_PI,
				39.2/180.0*M_PI,
				-49.2/180.0*M_PI,
				70.0/180.0*M_PI,
				-62.4/180.0*M_PI,
				80.2/180.0*M_PI,
				187.2/180.0*M_PI;
	robot->updateModel();

	// initialize GLFW window
	GLFWwindow* window = glfwInitialize();

    // set callbacks
	glfwSetKeyCallback(window, keySelect);

	// start the simulation
	thread sim_thread(simulation, robot);
	
    // while window is open:
    while (!glfwWindowShouldClose(window)) {

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

    // destroy context
    glfwDestroyWindow(window);

    // terminate
    glfwTerminate();

	return 0;
}

/* =======================================================================================
   SIMULATION SETUP
   -----------------------
   * Simulation loop
   * Window initialization
   * Window error
   * Mouse click commands
========================================================================================== */

void simulation(Sai2Model::Sai2Model* robot) {  //--- Simulation loop ---//
	fSimulationRunning = true;

	// create a timer
	LoopTimer timer;
	timer.initializeTimer();
	timer.setLoopFrequency(500); //500Hz timer
	double last_time = timer.elapsedTime(); //secs

	Eigen::MatrixXd Jbar;
	Eigen::MatrixXd Jv;
	Eigen::MatrixXd Jw;
	Eigen::MatrixXd J0;
	Eigen::VectorXd opspace_vel(6); // operational space velocity

	bool fTimerDidSleep = true;

	while (fSimulationRunning) {
		fTimerDidSleep = timer.waitForNextLoop();

		// integrate joint velocity to joint positions
		double curr_time = timer.elapsedTime();
		double loop_dt = curr_time - last_time; 
		robot->_q += robot->_dq*loop_dt;

		// update kinematic models
		robot->updateModel();

		/* --------------------------------------------------------------------------------------
			FILL ME IN: set new joint velocities
		-------------------------------------------------------------------------------------*/
    // set the desired cartesian position and velocities
    VectorXd xd(7), d_xd(7); 
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

    MatrixXd Ew_inv(3,4);
    double l_norm = pow(xd[3], 2) + pow(xd[4], 2) + pow(xd[5], 2) + pow(xd[6], 2);
    Ew_inv << -xd[4], xd[3], -xd[6], xd[5],
              -xd[5], xd[6], xd[3], -xd[4],
              -xd[6], -xd[5], xd[4], xd[3];
    Ew_inv = 2.*Ew_inv/l_norm;
    
    MatrixXd E_inv(6,7);
    E_inv.topLeftCorner(3, 3) = MatrixXd::Identity(3, 3);
    E_inv.bottomRightCorner(3, 4) = Ew_inv;
    // get cartesian velocities
    VectorXd d_zd = E_inv*d_xd;

    // get the task jacobian fist
    robot->J_0(J0, ee_link_name);
    MatrixXd Lambda, N;
    robot->operationalSpaceMatrices(Lambda, Jbar, N, (const MatrixXd) J0);
	  robot->_dq = Jbar*d_zd;

		// update last time
		last_time = curr_time;
	}
}

//------------------------------------------------------------------------------
GLFWwindow* glfwInitialize() { //--- Window initialization ---//
	
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
    GLFWwindow* window = glfwCreateWindow(windowW, windowH, "SAI2.0 - CS327a HW1", NULL, NULL);
	glfwSetWindowPos(window, windowPosX, windowPosY);
	glfwShowWindow(window);
    glfwMakeContextCurrent(window);
	glfwSwapInterval(1);

	return window;
}

//------------------------------------------------------------------------------

void glfwError(int error, const char* description) {  //--- Window error ---// 
	cerr << "GLFW Error: " << description << endl;
	exit(1);
}

//------------------------------------------------------------------------------

void keySelect(GLFWwindow* window, int key, int scancode, int action, int mods) //--- Mouse click commands ---//
{
    // option ESC: exit
    if (key == GLFW_KEY_ESCAPE && action == GLFW_PRESS)
    {
        // exit application
         glfwSetWindowShouldClose(window, 1);
    }
}
