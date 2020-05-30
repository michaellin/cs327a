/*======================================================================================
 * hw4 - p1-main-sol.cpp
 *
 * Implement cooperative manipulation
 * by completing the "FILL ME IN" section 
 *
 *
 * Elena Galbally, Spring 2019 
 * Wesley Guo, Spring 2020
 *======================================================================================*/

/* --------------------------------------------------------------------------------------
   Include Required Libraries and Files
-----------------------------------------------------------------------------------------*/
#include <iostream>
#include <string>
#include <csignal>
#include <thread>
#include <chrono>
#include <math.h>

#include "Sai2Model.h"
#include "Sai2Graphics.h"
#include "Sai2Simulation.h"
#include "redis/RedisClient.h"
#include "timer/LoopTimer.h"

#include <dynamics3d.h>
#include "timer/LoopTimer.h"
#include <GLFW/glfw3.h> 

#include "keys.h"

using namespace std;

const double time_slowdown_factor = 10;

constexpr const char *sim_title = "SAI2.0 - CS327a HW4 P1 Solution";
const string world_fname = "resources/hw4/world_1_puma.urdf";
const string robot_fname = "../resources/puma/puma_gripper.urdf";
const string robot1_name = "Puma1";
const string robot2_name = "Puma2";
const string object_name = "CoordObject";
const string object_fname = "resources/hw4/object.urdf";
const string object_link_name = "object";
const string camera_name = "camera_front";
const string ee_link_name = "end-effector";
const string gripper_joint_name = "gripper";

RedisClient redis_client;
/* ----------------------------------------------------------------------------------
	Simulation and Control Loop Setup
-------------------------------------------------------------------------------------*/
// state machine setup
enum ControlMode {
	CONTROL_GRASP_STABILIZE = 0,
	CONTROL_AUGMENTED_OBJECT
};

// simulation loop
bool fSimulationRunning = false;
void control(Sai2Model::Sai2Model* robot1, Sai2Model::Sai2Model* robot2, Sai2Model::Sai2Model* object_model, Simulation::Sai2Simulation* sim);
void simulation(Simulation::Sai2Simulation* sim);

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

    // set up redis callbacks
    redis_client.connect();
    redis_client.createReadCallback(0);
    redis_client.createWriteCallback(0);

	// load graphics scene
	auto graphics = new Sai2Graphics::Sai2Graphics(world_fname, false);

	// load robots
	auto robot1 = new Sai2Model::Sai2Model(robot_fname, false);
	auto robot2 = new Sai2Model::Sai2Model(robot_fname, false);

	// load object
	auto coobject = new Sai2Model::Sai2Model(object_fname, false);

	// load simulation world
	auto sim = new Simulation::Sai2Simulation(world_fname, false);
	// set co-efficient of restition to zero to avoid bounce
    // see issue: https://github.com/manips-sai/sai2-simulation/issues/1
    sim->setCollisionRestitution(0.0);
    sim->setCoeffFrictionStatic(0.5);
    sim->setCoeffFrictionDynamic(0.5);

    // set joint damping on grippers: 
    auto base_1 = sim->_world->getBaseNode(robot1_name);
    auto gripper_1 = base_1->getJoint(gripper_joint_name);
    gripper_1->setDamping(10.0);
    gripper_1->setJointLimits(-0.005, 0.068, 0.005);
    auto base_2 = sim->_world->getBaseNode(robot2_name);
    auto gripper_2 = base_2->getJoint(gripper_joint_name);
    gripper_2->setDamping(10.0);
    gripper_2->setJointLimits(-0.005, 0.068, 0.005);

    // set initial conditions
	robot1->_q << 90/180.0*M_PI,
				-22.5/180.0*M_PI,
				212/180.0*M_PI,
				90.0/180.0*M_PI,
				100/180.0*M_PI,
				180/180.0*M_PI,
				0.0;
	sim->setJointPositions(robot1_name, robot1->_q);
	robot1->updateModel();
	robot2->_q << 90/180.0*M_PI,
				202.5/180.0*M_PI,
				-28/180.0*M_PI,
				-90.0/180.0*M_PI,
				97/180.0*M_PI,
				180/180.0*M_PI,
				0.0;
	sim->setJointPositions(robot2_name, robot2->_q);
	robot2->updateModel();
	Eigen::Affine3d ee_trans;

	// set up error callback
    glfwSetErrorCallback(glfwError);

    // initialize GLFW
    glfwInit();

    // retrieve resolution of computer display and position window accordingly
    GLFWmonitor *primary = glfwGetPrimaryMonitor();
    const GLFWvidmode *mode = glfwGetVideoMode(primary);

    // information about computer screen and GLUT display window
    int screenW = mode->width;
    int screenH = mode->height;
    int windowW = 0.8 * screenH;
    int windowH = 0.5 * screenH;
    int windowPosY = (screenH - windowH) / 2;
    int windowPosX = windowPosY;

    // create window and make it current
    glfwWindowHint(GLFW_VISIBLE, 0);
    GLFWwindow *window = glfwCreateWindow(windowW, windowH, sim_title, NULL, NULL);
    glfwSetWindowPos(window, windowPosX, windowPosY);
    glfwShowWindow(window);
    glfwMakeContextCurrent(window);
    glfwSwapInterval(1);

    // set callbacks
    glfwSetKeyCallback(window, keySelect);

    // cache variables
    double last_cursorx, last_cursory;

	// start the simulation thread first
    fSimulationRunning = true;
	thread sim_thread(simulation, sim);

	// next start the control thread
	thread ctrl_thread(control, robot1, robot2, coobject, sim);
	
    // while window is open:
    while (!glfwWindowShouldClose(window)) {

		// update graphics. this automatically waits for the correct amount of time
		int width, height;
		glfwGetFramebufferSize(window, &width, &height);
		graphics->updateGraphics(robot1_name, robot1);
		graphics->updateGraphics(robot2_name, robot2);
		graphics->updateGraphics(object_name, coobject);
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

/* ----------------------------------------------------------------------------------
	Utility functions
-------------------------------------------------------------------------------------*/
// Calculate the cross product matrix
Eigen::Matrix3d getCrossProductMat(const Eigen::Vector3d& t) {
	Eigen::Matrix3d ret;
	ret <<  0, -t(2), t(1),
			t(2), 0, -t(0),
			-t(1), t(0), 0;
	return ret;
}

/* =======================================================================================
   CONTROL LOOP
========================================================================================== */
void control(Sai2Model::Sai2Model* robot1, Sai2Model::Sai2Model* robot2, Sai2Model::Sai2Model* object_model, Simulation::Sai2Simulation* sim) {
	
	// create a timer
	LoopTimer timer;
	timer.initializeTimer();
	timer.setLoopFrequency(1000); //1000Hz timer
	double last_time = timer.elapsedTime()/time_slowdown_factor; //secs

	// load robot global frame to robot base transformations: todo: move out to a function
	Eigen::Affine3d robot1_base_frame = sim->getRobotBaseTransform(robot1_name);
	Eigen::Affine3d robot2_base_frame = sim->getRobotBaseTransform(robot2_name);

	// cache variables
	bool fTimerDidSleep = true;
	bool fTorqueUseDefaults = false; // set true when torques are overriden for the first time
	Eigen::VectorXd tau1 = Eigen::VectorXd::Zero(robot1->dof());
	Eigen::VectorXd tau2 = Eigen::VectorXd::Zero(robot2->dof());

	Eigen::Affine3d object_com_frame;
	Eigen::Vector3d object_current_pos;
	Eigen::MatrixXd object_inertia(6,6);
	Eigen::MatrixXd object_j(6,6);
	Eigen::VectorXd object_p(6);

	Eigen::VectorXd robot1_g(robot1->dof());
	Eigen::VectorXd robot2_g(robot2->dof());

	// **  Other sugested variables and gains **
	Eigen::Vector3d object_com_in_robot1_ee_frame;
	Eigen::Vector3d object_com_in_robot2_ee_frame;
	Eigen::MatrixXd robot1_j0_objcom(6, robot1->dof());
	Eigen::MatrixXd robot2_j0_objcom(6, robot2->dof());
	Eigen::MatrixXd robot1_j0_objectcom_bar(robot1->dof(), 6);
	Eigen::MatrixXd robot2_j0_objectcom_bar(robot2->dof(), 6);
	Eigen::MatrixXd robot1_nbar_objectcom(robot1->dof(), robot1->dof());
	Eigen::MatrixXd robot2_nbar_objectcom(robot2->dof(), robot1->dof());
	Eigen::MatrixXd robot1_objcom_inertia(6,6);
	Eigen::MatrixXd robot2_objcom_inertia(6,6);

	Eigen::MatrixXd augmented_object_inertia(6,6);
	Eigen::VectorXd augmented_object_p(6);

	Eigen::MatrixXd G(2*6, 2*6);
	// Eigen::MatrixXd W(6, 2*6);

	// Eigen::Vector3d obj_des_pos;
	// Eigen::Vector3d obj_ori_error;
	// Eigen::VectorXd obj_task_err(6);
	Eigen::VectorXd force_des_vec(12);
	Eigen::VectorXd force_ee_vec(12);

	double kp = 30;
	double kv = 10;

	// ** Control Mode **
	// 		0 = grasp stabilizing controller
	//		1 = augmented object controller
	ControlMode control_mode = CONTROL_GRASP_STABILIZE; 

	while (fSimulationRunning) { //automatically set to false when simulation is quit
		fTimerDidSleep = timer.waitForNextLoop();

		// update time
		double curr_time = timer.elapsedTime()/time_slowdown_factor;
		double loop_dt = curr_time - last_time;

        // read joint positions, velocities
        sim->getJointPositions(robot1_name, robot1->_q);
        sim->getJointVelocities(robot1_name, robot1->_dq);
        robot1->updateModel();
        sim->getJointPositions(robot2_name, robot2->_q);
        sim->getJointVelocities(robot2_name, robot2->_dq);
        robot2->updateModel();

        // read object position
        sim->getJointPositions(object_name, object_model->_q);
        sim->getJointVelocities(object_name, object_model->_dq);
        object_model->updateModel();

        // update object dynamics
		// - find object COM frame in global frame
		object_model->transform(object_com_frame, object_link_name);
		object_current_pos << object_com_frame.translation();
		redis_client.addEigenToWriteCallback(0, CURRENT_OBJECT_POSITION_KEY, object_current_pos);
		// - obtain object inertia matrix in world frame
		object_model->J_0(object_j, object_link_name, Eigen::Vector3d::Zero());
		object_inertia = (object_j*object_model->_M_inv*object_j.transpose()).inverse();
		// - obtain object p
		object_p << 0, 0, 9.8*0.5, 0, 0, 0;

        // ---------------------------------------------------------------------------------
        /* ---------------------- FILL ME IN ----------------------------------------------- */
		
		// --------------------------------------------------------------------
		// (1) Manipulator Jacobians
		//---------------------------------------------------------------------
		// ** FILL ME IN **
    // get the translation of the object COM in the end-effector frame
    // robot 1 jacobian
    Eigen::Affine3d robot1_T_0_ee;
    robot1->transform(robot1_T_0_ee, ee_link_name);
    Eigen::Affine3d robot1_T_W_ee = robot1_base_frame * robot1_T_0_ee;
    Eigen::Affine3d robot1_T_ee_objcom = robot1_T_W_ee.inverse() * object_com_frame;
    object_com_in_robot1_ee_frame = robot1_T_ee_objcom.translation();
    robot1->J_0(robot1_j0_objcom, ee_link_name, object_com_in_robot1_ee_frame);


    //std::cout << "object com " << object_com_frame.translation() << std::endl;
    //std::cout << "robot 1 relative rot" << robot1_T_ee_objcom.linear() << std::endl;
    //std::cout << "robot1 ee pos " << robot1_T_W_ee.translation() << std::endl;

    // robot 2 jacobian
    Eigen::Affine3d robot2_T_0_ee;
    robot2->transform(robot2_T_0_ee, ee_link_name);
    Eigen::Affine3d robot2_T_W_ee = robot2_base_frame * robot2_T_0_ee;
    Eigen::Affine3d robot2_T_ee_objcom = robot2_T_W_ee.inverse() * object_com_frame;
    object_com_in_robot2_ee_frame = robot2_T_ee_objcom.translation();
    robot2->J_0(robot2_j0_objcom, ee_link_name, object_com_in_robot2_ee_frame);

    //std::cout << "robot2 ee pos " << robot2_T_W_ee.translation() << std::endl;
    //std::cout << "robot 2 relative rot" << robot2_T_ee_objcom.linear() << std::endl;

		// --------------------------------------------------------------------
		// (2) Augmented Object Model
		//---------------------------------------------------------------------
		// ** FILL ME IN **
    robot1->operationalSpaceMatrices(robot1_objcom_inertia, 
                                    robot1_j0_objectcom_bar,
                                    robot1_nbar_objectcom,
                                    robot1_j0_objcom);
    robot2->operationalSpaceMatrices(robot2_objcom_inertia, 
                                    robot2_j0_objectcom_bar,
                                    robot2_nbar_objectcom,
                                    robot2_j0_objcom);
    augmented_object_inertia = object_inertia
                              + robot1_objcom_inertia
                              + robot2_objcom_inertia;

    //std::cout << "auginertia " << augmented_object_inertia << std::endl;
    //std::cout << "objinertia " << object_inertia<< std::endl;
    //std::cout << "rob1inertia " << robot1_objcom_inertia<< std::endl;
    //std::cout << "rob2inertia " << robot2_objcom_inertia<< std::endl;

    // get the gravity vector for each robot
    Eigen::VectorXd robot1_g(robot1->dof());
    Eigen::VectorXd robot1_p(6);
    Eigen::VectorXd robot2_g(robot2->dof());
    Eigen::VectorXd robot2_p(6);
    robot1->gravityVector(robot1_g);
    robot1_p = robot1_objcom_inertia*robot1_j0_objcom*robot1->_M_inv*robot1_g;
    robot2->gravityVector(robot2_g);
    robot2_p = robot2_objcom_inertia*robot2_j0_objcom*robot2->_M_inv*robot2_g;
    augmented_object_p = object_p
                          + robot1_p
                          + robot2_p;

		// --------------------------------------------------------------------
		// (3) Grasp Matrix
		//---------------------------------------------------------------------
		// ** FILL ME IN **
    // external forces matrices
    Eigen::MatrixXd Wf(6, 6);
    Eigen::MatrixXd Wm(6, 6);
    Eigen::Matrix3d robot1_R = robot1_T_W_ee.linear();
    Eigen::Matrix3d robot2_R = robot2_T_W_ee.linear();
    //robot1_R << 0., -1., 0.,
     //           -1., 0., 0.,
      //           0., 0., -1.;

    //robot2_R << 0., 1., 0.,
     //            1., 0., 0.,
      //           0., 0., -1.;
    Eigen::Vector3d robot1_rhat;
    //robot1_rhat << 0., -1., 0.;
    robot1_rhat  = robot1_T_W_ee.translation() - object_com_frame.translation();
    //robot1_rhat  = robot1_T_ee_objcom.inverse().translation();
    //Eigen::Vector3d robot1_rhat = (object_com_frame.inverse() * robot1_T_W_ee).translation();
    //std::cout << "robot1 rhat " << robot1_rhat << std::endl;

    Eigen::Vector3d robot2_rhat;
    //robot2_rhat << 0., 1., 0.;
    robot2_rhat  = robot2_T_W_ee.translation() - object_com_frame.translation();
    //robot2_rhat  = robot2_T_ee_objcom.inverse().translation();
    //Eigen::Vector3d robot2_rhat = (object_com_frame.inverse() * robot2_T_W_ee).translation();
    //std::cout << "robot2 rhat " << robot2_rhat << std::endl;

    Wf.topLeftCorner(3, 3) = robot1_R;
    Wf.topRightCorner(3, 3) = robot2_R;
    Wf.bottomLeftCorner(3, 3) = getCrossProductMat(-robot1_rhat);
    Wf.bottomRightCorner(3, 3) = getCrossProductMat(-robot2_rhat);

    Wm.topLeftCorner(3, 3) = Matrix3d::Zero();
    Wm.topRightCorner(3, 3) = Matrix3d::Zero();
    Wm.bottomLeftCorner(3, 3) = robot1_R;
    Wm.bottomRightCorner(3, 3) = robot2_R;

    // internal forces matrices
    Eigen::Matrix3d W_objcom_R = object_com_frame.linear();
    //Vector3d e12 = W_objcom_R*(robot2_T_W_ee.translation() - robot1_T_W_ee.translation());
    Vector3d e12 = (robot2_T_W_ee.translation() - robot1_T_W_ee.translation());
    e12 = e12/e12.norm();
    Eigen::MatrixXd Ebar_row(1, 12);
    //Ebar_row << 0.5*e12(0), 0.5*e12(1), 0.5*e12(2), -0.5*e12(0), -0.5*e12(1), -0.5*e12(2),
    //            0., 0., 0., 0., 0., 0.;
    Ebar_row << 0.5, 0, 0, 0.5, 0, 0,
                0., 0., 0., 0., 0., 0.;

    G.topLeftCorner(6, 6) = Wf;
    G.topRightCorner(6, 6) = Wm;
    G.row(6) = Ebar_row;
    Vector3d e1z, e2z;
    //e1z << 0., 0., -1.;
    e1z << robot1_T_W_ee(0,2), 
          robot1_T_W_ee(1,2), 
          robot1_T_W_ee(2,2);
    //e2z << 0., 0., -1.;
    e2z << robot2_T_W_ee(0,2), 
          robot2_T_W_ee(1,2), 
          robot2_T_W_ee(2,2);
    Vector3d e1y = e1z.cross(e12);
    Vector3d e2y = e2z.cross(-e12);
    //G.bottomRightCorner(5,6) << 
    //        0.5*e12(0), 0.5*e12(1), 0.5*e12(2), -0.5*e12(0), -0.5*e12(1), -0.5*e12(2),
    //        e1y(0), e1y(1), e1y(2), 0., 0., 0.,
    //        e1z(0), e1z(1), e1z(2), 0., 0., 0.,
    //        0., 0., 0., e2y(0), e2y(1), e2y(2),
    //        0., 0., 0., e2z(0), e2z(1), e2z(2);

    G.bottomRightCorner(5,6) << 
            0.5*e12(0), 0.5*e12(1), 0.5*e12(2), -0.5*e12(0), -0.5*e12(1), -0.5*e12(2),
            e1y(0), e1y(1), e1y(2), 0., 0., 0.,
            e1z(0), e1z(1), e1z(2), 0., 0., 0.,
            0., 0., 0., e2y(0), e2y(1), e2y(2),
            0., 0., 0., e2z(0), e2z(1), e2z(2);

    //vector<Vector3d> contact_locations;
    //contact_locations.push_back(robot1_T_W_ee);
    //contact_locations.push_back(robot2_T_W_ee);
    //vector<Sai2Model::ContactType> contact_types;
    //contact_types.push_back(Sai2Model::SurfaceContact);
    //contact_types.push_back(Sai2Model::SurfaceContact);
    //object_model->graspMatrix(G, Ginv, object_com_frame.linear(), object_com_frame.translation(),
    //            contact_locations, contact_types)

    //std::cout << "Ginv:\n " << G.inverse() << std::endl;
    std::cout << "G:\n " << G << std::endl;
		// --------------------------------------------------------------------
		// (4) Force Control 
		//---------------------------------------------------------------------
		if (control_mode == CONTROL_AUGMENTED_OBJECT) {
		
      // ** FILL ME IN ** Compute tau1 and tau2 	

      // first define the desired position and orientation
      VectorXd xd(7);
      xd << 0.0,
            //0.15 * sin(2*M_PI/3 * curr_time),
            0.0,
            0.4,
            1.0,
            0.0,
            0.0,
            0.0;

      VectorXd d_xd = VectorXd::Zero(7);
      //d_xd(1) = 0.15 * (2*M_PI/3) * cos(2*M_PI/3 * curr_time);
      VectorXd dd_xd = VectorXd::Zero(7);
      //dd_xd(1) = -0.15 * pow(2*M_PI/3,2) * sin(2*M_PI/3 * curr_time);

      // get object rotation matrix
      Matrix3d rot;
      object_model->rotation(rot, object_link_name);

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

      // get object linear position
      Eigen::Vector3d obj_x;
      object_model->position(obj_x, object_link_name);
      // get object linear velocity
      Eigen::Vector3d obj_dx;
      object_model->linearVelocity(obj_dx, object_link_name);
      // get delta phi
      Eigen::Vector3d delta_phi;
      delta_phi = -Er_inv * xd.tail(4);

      VectorXd f_star(3), m_star(3), f_motion(6);
      VectorXd f_ext = VectorXd::Zero(7);
      //f_star = dd_xd.head(3) - kp * (obj_x - xd.head(3)) - kv * (obj_dx - d_xd.head(3));
      f_star = - kp * (obj_x - xd.head(3));
      //f_star = VectorXd::Zero(3);
      m_star = - kp * delta_phi;
      //m_star = VectorXd::Zero(3);
      f_motion.head(3) = f_star; f_motion.tail(3) = m_star;   // put together wrench

      f_ext = augmented_object_inertia * f_motion + object_p;
      //std::cout << "fmotion" << f_motion << std::endl;
      
      force_des_vec << f_ext, -15, 0, 0, 0, 0, 0;

      std::cout << "fdes " << force_des_vec << std::endl;

      // compute the end effector forces with grasp matrix
      force_ee_vec = G.inverse()*force_des_vec;

      VectorXd f1_ee(6), f2_ee(6);
      f1_ee << force_ee_vec(0), force_ee_vec(1), force_ee_vec(2),
               force_ee_vec(6), force_ee_vec(7), force_ee_vec(8);
      f2_ee << force_ee_vec(3), force_ee_vec(4), force_ee_vec(5),
               force_ee_vec(9), force_ee_vec(10), force_ee_vec(11);
      std::cout << "f1ee " << f1_ee << std::endl;
      std::cout << "f2ee " << f2_ee << std::endl;
      //std::cout << "fee " << force_ee_vec << std::endl;
      //std::cout << "G " << G << std::endl;
      //std::cout << "force_des " << G*force_ee_vec << std::endl;

      // get robot end effector jacobians
      MatrixXd robot1_j0_ee(6, robot1->dof());
      MatrixXd robot2_j0_ee(6, robot2->dof());
			robot1->J_0(robot1_j0_ee, ee_link_name);
			robot2->J_0(robot2_j0_ee, ee_link_name);
      tau1 = robot1_j0_ee.transpose()*(f1_ee) + robot1_g;
      tau2 = robot2_j0_ee.transpose()*(f2_ee) + robot2_g;

		} else if (control_mode == CONTROL_GRASP_STABILIZE) { // initial grasp stabilization
			
			Eigen::MatrixXd robot1_j0_ee(6, robot1->dof());
			Eigen::MatrixXd robot2_j0_ee(6, robot2->dof());
			robot1->J_0(robot1_j0_ee, ee_link_name, Eigen::Vector3d::Zero());
			robot2->J_0(robot2_j0_ee, ee_link_name, Eigen::Vector3d::Zero());

			// Joint Torques
			tau1 = robot1_j0_ee.transpose()*(object_p/2) + robot1->_M*(-10.0*robot1->_dq) + robot1_g;
			tau2 = robot2_j0_ee.transpose()*(object_p/2) + robot2->_M*(-10.0*robot2->_dq) + robot2_g;
			
			// Grasp stabilization
			static uint grasp1Counter = 0;
			static uint grasp2Counter = 0;
			if (robot1->_dq[6] < 0.1) {
				grasp1Counter += 1;
			} else {
				grasp1Counter = 0;
			}
			if (robot2->_dq[6] < 0.1) {
				grasp2Counter += 1;
			} else {
				grasp2Counter = 0;
			}
			if (grasp1Counter > 40 && grasp2Counter > 40) {
				cout << " ** Switch Control Mode to Augmented Object Model ** " << endl;
				control_mode = CONTROL_AUGMENTED_OBJECT;
			}
		}

		/* ----------------------------------------------------------------------------------
			Safety torques 
		-------------------------------------------------------------------------------------*/ 
		
		// Set constant gripper forces
		tau1[6] = 15;
		tau2[6] = 15;

        // Default values if torques are exceeded:
        bool fTorqueOverride = false; // to avoid robot blow-ups
        const double tau1_max = 200;
        const double tau2_max = 200;
        if (!fTorqueUseDefaults) {
        	if (tau1.cwiseAbs().maxCoeff() > tau1_max || tau2.cwiseAbs().maxCoeff() > tau2_max) {
	        	fTorqueOverride = true;
	        	cerr << "Torque overriden. User asked torques beyond safe limit: \n";
	        	cerr << "Robot 1: " << tau1.transpose() << "\n";
	        	cerr << "Robot 2: " << tau2.transpose() << "\n";
	        	fTorqueUseDefaults = true;
	        }
	        // Also default values if object is dropped
	        const double object_thickness = 0.05;
	        bool fRobot1DroppedObject = robot1->_q[6] > object_thickness/2;
	        bool fRobot2DroppedObject = robot2->_q[6] > object_thickness/2;
	        if (fRobot1DroppedObject || fRobot2DroppedObject) {
	        	cerr << "Torque overriden. Robot 1 or 2 dropped object. \n";
	        	fTorqueUseDefaults = true;
	        }
        }
        else {
        	robot1->gravityVector(tau1);
          tau1 = tau1 + robot1->_M*(-10.0*robot1->_dq);
          tau1 = (tau1.array() >= tau1_max).select(tau1_max*Eigen::VectorXd::Ones(robot1->dof()), 
                                                                                  tau1);
          tau1 = (tau1.array() <= -tau1_max).select(-tau1_max*Eigen::VectorXd::Ones(robot1->dof()), 
                                                                                    tau1);
          robot2->gravityVector(tau2);
          tau2 = tau2 + robot2->_M*(-10.0*robot2->_dq);
          tau2 = (tau2.array() >= tau2_max).select(tau2_max*Eigen::VectorXd::Ones(robot2->dof()), 
                                                                                  tau2);
          tau2 = (tau2.array() <= -tau2_max).select(-tau2_max*Eigen::VectorXd::Ones(robot2->dof()), 
                                                                                    tau2);
        }

        /* ----------------------------------------------------------------------------------
			Send robot torques from controller
		-------------------------------------------------------------------------------------*/ 
		sim->setJointTorques(robot1_name, tau1);
		sim->setJointTorques(robot2_name, tau2);
		
		// write object location key out to redis
		redis_client.executeWriteCallback(0);

		// update last time
		last_time = curr_time;
	}
}

/* =======================================================================================
   SIMULATION SETUP
   -----------------------
   * Simulation loop
   * Select trajectory
   * Window initialization
   * Window error
   * Mouse click commands
========================================================================================== */
void simulation(Simulation::Sai2Simulation* sim) {
	fSimulationRunning = true;

	// create a timer
	LoopTimer timer;
	timer.initializeTimer();
	timer.setLoopFrequency(1000); //10000Hz timer

	double last_time = timer.elapsedTime()/time_slowdown_factor; //secs
	bool fTimerDidSleep = true;
	while (fSimulationRunning) {
		fTimerDidSleep = timer.waitForNextLoop();
		// integrate forward
		double curr_time = timer.elapsedTime()/time_slowdown_factor;
		double loop_dt = curr_time - last_time; 
		// sim->integrate(loop_dt);
		sim->integrate(loop_dt);
		// update last time
		last_time = curr_time;
	}
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
