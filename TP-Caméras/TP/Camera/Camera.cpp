#include <TP/Camera/Camera.hpp>
#include <TP/Camera/Camera_Helper.hpp>
#include <glm/gtx/quaternion.hpp>

#include <iostream>
#include <algorithm>

// Include GLM
#include <imgui/imgui.h>
#include <imgui/imgui_impl_glfw.h>
#include <imgui/imgui_impl_opengl3.h>
#include "Camera.hpp"


void Camera::init()
{
	resetCamera();
}

void Camera::updateInterface(float _deltaTime)
{
	// ImGUI window creation
	if (ImGui::Begin("Interface"))
	{
        ImGui::Separator();
		ImGui::Checkbox("Track Target##target", &m_isTrackingTarget);

		ImGui::Separator();
 		ImGui::Text("Camera Position : ");
        ImGui::InputFloat3("Position", &m_position[0]);
        ImGui::InputFloat("Distance avec la target", &m_distanceBehind);


		ImGui::Separator();
        ImGui::Text("Camera Orientation : ");
		ImGui::InputFloat("Pitch", &m_eulerAngle.x);
        ImGui::InputFloat("Yaw", &m_eulerAngle.y);  
        Camera_Helper::clipAngle180(m_eulerAngle.y);
		m_rotation = glm::quat(glm::vec3(glm::radians(m_eulerAngle.x), glm::radians(m_eulerAngle.y), glm::radians(m_eulerAngle.z)));
	

        ImGui::Separator();
        ImGui::Text("Camera Speed : ");
        ImGui::InputFloat("Translation Speed", &m_translationSpeed);
        ImGui::InputFloat("Rotation Speed", &m_rotationSpeed);

		ImGui::Separator();
		if (ImGui::Button("Reset Camera")){
			startResetAnimation();
		}

		ImGui::Separator();
        ImGui::Text("Movement Mode : ");
        const char* modes[] = { "Mode 0", "Mode 1" };
        int currentMode = m_mode;
        if (ImGui::Combo("Mode", &currentMode, modes, IM_ARRAYSIZE(modes)))
        {
            onModeChanged(currentMode);
        }

        ImGui::Separator();
        ImGui::Text("Animation : ");
        ImGui::InputFloat("Animation Length", &m_animationLength);
        

		ImGui::Separator();
		ImGui::Checkbox("show Imgui Demo", &m_showImguiDemo);
	}
	ImGui::End();

	if (m_showImguiDemo)
	{
		ImGui::Text("hey");
		ImGui::ShowDemoWindow();
	}

}

void Camera::updateFreeInput(float _deltaTime, GLFWwindow* _window)
{

    if(m_isTrackingTarget){
        return;
    }

	if (m_mode == 0){
		if (glfwGetKey(_window, GLFW_KEY_W) == GLFW_PRESS) {
		m_position += m_rotation * glm::vec3(0.f, 0.f, 1.f) * m_translationSpeed;
		}
		if (glfwGetKey(_window, GLFW_KEY_A) == GLFW_PRESS) {
			m_position += m_rotation * glm::vec3(1.f, 0.f, 0.f) * m_translationSpeed;
		}
		if (glfwGetKey(_window, GLFW_KEY_S) == GLFW_PRESS) {
			m_position += m_rotation * glm::vec3(0.f, 0.f, -1.f) * m_translationSpeed;
		}
		if (glfwGetKey(_window, GLFW_KEY_D) == GLFW_PRESS) {
			m_position += m_rotation * glm::vec3(-1.f, 0.f, 0.f) * m_translationSpeed;
		}
	}
	else if (m_mode == 1){

        glm::vec3 forward = m_rotation * glm::vec3(0.0f, 0.0f, 1.0f);
        forward.y = 0.0f; 
        forward = glm::normalize(forward);
        
        glm::vec3 right = glm::cross(forward, glm::vec3(0.0f, 1.0f, 0.0f));

        if (glfwGetKey(_window, GLFW_KEY_W) == GLFW_PRESS) {
            m_position += forward * m_translationSpeed;
        }
        if (glfwGetKey(_window, GLFW_KEY_S) == GLFW_PRESS) {
            m_position -= forward * m_translationSpeed;
        }
        if (glfwGetKey(_window, GLFW_KEY_D) == GLFW_PRESS) {
            m_position += right * m_translationSpeed;
        }
        if (glfwGetKey(_window, GLFW_KEY_A) == GLFW_PRESS) {
            m_position -= right * m_translationSpeed;
        }
        
        if (glfwGetKey(_window, GLFW_KEY_Q) == GLFW_PRESS) {
            m_position.y -= m_translationSpeed;
        }
        if (glfwGetKey(_window, GLFW_KEY_E) == GLFW_PRESS) {
            m_position.y += m_translationSpeed;
        }

        if (glfwGetKey(_window, GLFW_KEY_RIGHT) == GLFW_PRESS) {
            m_eulerAngle.y -= m_rotationSpeed;
            Camera_Helper::clipAngle180(m_eulerAngle.y);
        }
        if (glfwGetKey(_window, GLFW_KEY_LEFT) == GLFW_PRESS) {
            m_eulerAngle.y += m_rotationSpeed;
            Camera_Helper::clipAngle180(m_eulerAngle.y);
        }

        if (glfwGetKey(_window, GLFW_KEY_UP) == GLFW_PRESS) {
            m_eulerAngle.x -= m_rotationSpeed;
            m_eulerAngle.x = std::clamp(m_eulerAngle.x, -90.0f, 90.0f);
        }
        if (glfwGetKey(_window, GLFW_KEY_DOWN) == GLFW_PRESS) {
            m_eulerAngle.x += m_rotationSpeed;
            m_eulerAngle.x = std::clamp(m_eulerAngle.x, -90.0f, 90.0f);
        }
        

	}
	

}

void Camera::keyCallback(GLFWwindow* _window, int _key, int _scancode, int _action, int _mods)
{
    Camera* camera = static_cast<Camera*>(glfwGetWindowUserPointer(_window));

    if (_key == GLFW_KEY_Z && _action == GLFW_PRESS)
    {
		int newMode = camera->getMode() + 1;
        if (newMode > 1)
        {
            newMode = 0;
        }
		camera->setMode(newMode);
		camera->onModeChanged(newMode);
    }

}

void Camera::mouseCallback(GLFWwindow* _window, double _xpos, double _ypos)
{
    Camera* camera = static_cast<Camera*>(glfwGetWindowUserPointer(_window));

    if (camera->getMode() != 0) {
        return; 
    }

    if (camera->m_isTrackingTarget)
    {
        return;
    }
    

    if (camera->m_firstMouse){
        camera->m_lastX = _xpos;
        camera->m_lastY = _ypos;
        camera->m_firstMouse = false;
    }

    float xoffset = camera->m_lastX - _xpos;
    float yoffset = _ypos - camera->m_lastY; 
    camera->m_lastX = _xpos;
    camera->m_lastY = _ypos;

    xoffset *= camera->m_rotationSpeed;
    yoffset *= camera->m_rotationSpeed;

    camera->m_eulerAngle.y += xoffset;
    camera->m_eulerAngle.x += yoffset;

    camera->m_eulerAngle.x = std::clamp(camera->m_eulerAngle.x, -90.0f, 90.0f);

    camera->m_rotation = glm::quat(glm::vec3(glm::radians(camera->m_eulerAngle.x), glm::radians(camera->m_eulerAngle.y), glm::radians(camera->m_eulerAngle.z)));
}



void Camera::glfw_cursor_enter_callback(GLFWwindow* _window, int _entered)
{
    ImGuiIO& io = ImGui::GetIO();
    if (!_entered){
        io.MousePos = ImVec2(-FLT_MAX, -FLT_MAX);
    }
}

void Camera::update(float _deltaTime, GLFWwindow* _window)
{
    if (m_animating) {
        glfwSetCursorPosCallback(m_window, NULL);
        float currentTime = static_cast<float>(glfwGetTime()) * 1000.0f; 
        float elapsedTime = currentTime - m_animationStartTime;

        if (elapsedTime < m_animationLength) {
            float t = std::clamp(elapsedTime / m_animationLength, 0.0f, 1.0f);
            t = Camera_Helper::interpolate(t, InterpolationType::EXP); 

            m_position = (m_targetPosition - m_initialPosition) * t + m_initialPosition;


            m_rotation = glm::slerp(m_initialRotation, m_targetRotation, t);
            m_rotation = glm::normalize(m_rotation);
            m_eulerAngle = glm::degrees(glm::eulerAngles(m_rotation)); 

        } else {
            m_position = m_targetPosition;
            m_rotation = m_targetRotation;
            m_rotation = glm::normalize(m_rotation);
            m_eulerAngle = glm::degrees(glm::eulerAngles(m_rotation)); 
            Camera_Helper::clipAngle180(m_eulerAngle.y);
            m_eulerAngle.x = std::clamp(m_eulerAngle.x, -90.0f, 90.0f);
            m_animationStartTime = 0;
            m_animating = false;

            centerMouse();

            glfwSetCursorPosCallback(m_window, Camera::mouseCallback);

        }
    }else {
        updateFreeInput(_deltaTime, _window);
    }

    updateInterface(_deltaTime);


	Camera_Helper::computeFinalView(m_projectionMatrix, m_viewMatrix, m_position, m_rotation, m_fovDegree);
}

void Camera::onModeChanged(int _newMode)
{
    setMode(_newMode);
    centerMouse();
    if (_newMode == 0)
    {
        glfwSetInputMode(m_window, GLFW_CURSOR, GLFW_CURSOR_DISABLED);
    }
    else
    {
        glfwSetInputMode(m_window, GLFW_CURSOR, GLFW_CURSOR_NORMAL);
    }
}

void Camera::resetCamera(){
	m_fovDegree = 45.0f;
	m_showImguiDemo = false;
	m_translationSpeed = 0.1f;
	m_rotationSpeed = 0.1f;
	m_mode = 0;

	glfwSetKeyCallback(m_window, Camera::keyCallback);
    glfwSetCursorPosCallback(m_window, Camera::mouseCallback);
    glfwSetWindowUserPointer(m_window, this);
    glfwSetCursorEnterCallback(m_window, Camera::glfw_cursor_enter_callback);

    glfwSetInputMode(m_window, GLFW_CURSOR, GLFW_CURSOR_DISABLED);

}

void Camera::startResetAnimation() {
    resetCamera();
    m_animationStartTime = static_cast<float>(glfwGetTime()) * 1000.0f; // Convertir en millisecondes
    m_initialPosition = m_position;
    m_targetPosition = glm::vec3(0.f, 40.f, 0.f);
    m_distanceBehind = 5.0f;
    m_initialRotation = m_rotation;
    m_targetRotation = glm::quat{};
    m_animating = true;
}


void Camera::centerMouse(){
    int width, height;
    glfwGetWindowSize(m_window, &width, &height);
    glfwSetCursorPos(m_window,  width/ 2, height / 2);
    m_lastX = width / 2;
    m_lastY = height / 2;
}

void Camera::trackActor(const glm::vec3& _position, const glm::quat& _rotation){
    glm::vec3 front = _rotation * glm::vec3(0.0f, 0.0f, 1.0f);
    m_position = _position - front * m_distanceBehind; 
    m_rotation = _rotation; 
    m_eulerAngle = glm::degrees(glm::eulerAngles(m_rotation));
}

