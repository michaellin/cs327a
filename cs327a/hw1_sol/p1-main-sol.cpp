/*======================================================================================
 * hw1-main-sol.cpp
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

const string world_fname = "resources/hw1_sol/world.urdf";
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

		// (1) Desired velocities (linear and angular)
		//----------------------------------------------------
		double R = 0.1; //m
		double T = 5.0; //sec

		// desired linear velocities
		double dx = 0;
		double dy = -2.0*M_PI*R/T*sin(2.0*M_PI*curr_time/T);
		double dz = 2.0*M_PI*R/T*sin(4.0*M_PI*curr_time/T);
		
		//desired orientations
		Eigen::Quaterniond lambda ( 
			1/sqrt(2.0)*sin(M_PI/4.0*cos(2.0*M_PI*curr_time/T)),
			1/sqrt(2.0)*cos(M_PI/4.0*cos(2.0*M_PI*curr_time/T)),
			1/sqrt(2.0)*sin(M_PI/4.0*cos(2.0*M_PI*curr_time/T)),
			1/sqrt(2.0)*cos(M_PI/4.0*cos(2.0*M_PI*curr_time/T)));

		//desired angular velocities
		Eigen::Quaterniond dlambda ( 
			-pow(M_PI, 2)/(2*T*sqrt(2.0))*cos(M_PI/4.0*cos(2.0*M_PI*curr_time/T))*sin(2.0*M_PI*curr_time/T),
			pow(M_PI, 2)/(2*T*sqrt(2.0))*sin(M_PI/4.0*cos(2.0*M_PI*curr_time/T))*sin(2.0*M_PI*curr_time/T),
			-pow(M_PI, 2)/(2*T*sqrt(2.0))*cos(M_PI/4.0*cos(2.0*M_PI*curr_time/T))*sin(2.0*M_PI*curr_time/T),
			pow(M_PI, 2)/(2*T*sqrt(2.0))*sin(M_PI/4.0*cos(2.0*M_PI*curr_time/T))*sin(2.0*M_PI*curr_time/T));
		
		// (2) Desired operational space velocities
		// ----------------------------------------------------
		// Calculate omega. 
		// 		Option 1 --> using the quaternion form of the E matrix: omega = 2 * dlambda * lambda.conjugate()
		// 		Option 2 --> omega = Er_inv * dlambda
		
		Eigen::Vector3d omega = 2*(dlambda*lambda.conjugate()).vec(); 
		opspace_vel << dx, dy, dz, omega.x(), omega.y(), omega.z();
		
		// (3) Desired joint velocities
		// ----------------------------------------------------
		// Obtain the basic Jacobian at the end-effector, J_0 = [Jv' Jw']'
		robot->J_0(J0, ee_link_name, Eigen::Vector3d::Zero());

		// Calculate the inertia weighted pseudo-inverse of J0
		Jbar = robot->_M_inv * J0.transpose() * (J0 * robot->_M_inv * J0.transpose()).inverse();
		
		// Finally, calculate the desired joint velocities needed to achieve the desired opspace_vel
		robot->_dq = Jbar*opspace_vel;

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
    GLFWwindow* window = glfwCreateWindow(windowW, windowH, "SAI2.0 - CS327a HW1 - SOLUTION", NULL, NULL);
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
