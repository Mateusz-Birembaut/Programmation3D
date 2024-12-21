#pragma once

#include <GLFW/glfw3.h>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/quaternion.hpp>
#include <glm/gtx/quaternion.hpp>
#include <glm/gtx/euler_angles.hpp>
#include <glm/gtx/norm.hpp>

class Camera
{
public: 

	void init();
	void resetCamera();

    void update(float _deltaTime, GLFWwindow *_window);
    void onModeChanged(int m_newMode);
    void updateInterface(float _deltaTime);
    void updateFreeInput(float _deltaTime, GLFWwindow* _window);

    void computeFinalView();
	void centerMouse();

	void startResetAnimation();

    static void glfw_cursor_enter_callback(GLFWwindow *window, int entered);	
	static void keyCallback(GLFWwindow* window, int key, int scancode, int action, int mods);
	static void mouseCallback(GLFWwindow* window, double xpos, double ypos);

	void trackActor(const glm::vec3& _position, const glm::quat& _rotation);

    glm::quat getRotation() const {return m_rotation;}
	glm::mat4 getViewMatrix() const {return m_viewMatrix;}
	glm::mat4 getProjectionMatrix() const {return m_projectionMatrix;}

	int getMode() const {return m_mode;}
	void setMode(int _mode) {m_mode = _mode;}

	bool getIsTrackingTarget() const {return m_isTrackingTarget;}

	void setWindow(GLFWwindow* window) { m_window = window; }
	GLFWwindow* getWindow() const { return m_window; }

private:

	//Camera parameters 
	float		m_fovDegree{ 45.0f };
	glm::vec3	m_position{ glm::vec3(0.f, 40.f, 0.f) };
	glm::vec3	m_eulerAngle{ glm::vec3(0.f, 0.f, 0.f) };
	glm::quat	m_rotation{};

	float m_distanceBehind { 5.0f};

	float m_translationSpeed { 1.f }; 
    float m_rotationSpeed { 1.f };

	int m_mode { 0 };

	//Animation
	bool m_animating { false };
	glm::vec3 m_initialPosition;
    glm::vec3 m_targetPosition;
	glm::quat m_initialRotation;
	glm::quat m_targetRotation;
    float m_animationStartTime;
	float m_animationLength { 1000.f };

	GLFWwindow* m_window;

	float m_lastX, m_lastY;
    bool m_firstMouse;

	//Interface option
	bool m_showImguiDemo{ false };

	bool m_isTrackingTarget{ false };

	//View
	glm::mat4 m_viewMatrix;
	glm::mat4 m_projectionMatrix;
};