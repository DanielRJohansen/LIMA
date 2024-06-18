#include "DisplayV2.h"

#ifndef __linux__

const float deg2rad = 2.f * PI / 360.f;
const float rad2deg = 1.f / deg2rad;
// Function to set up the projection
void setupProjection(int screenWidth, int screenHeight) {
     // Set up a perspective projection
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    double aspectRatio = 1.0;  // Assuming a window size of 640x480
    double fovY = 45.0;  // Field of view angle in the y direction (degrees)
    double nearPlane = 0.1;
    double farPlane = 10.0;
    double fH = tan(fovY / 360.0 * 3.14) * nearPlane;
    double fW = fH * aspectRatio;
    glFrustum(-fW, fW, -fH, fH, nearPlane, farPlane);
}

void Display::updateCamera(float pitch, float yaw, float delta_distance) {

    // Reset previous rotations
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    camera_distance += delta_distance;
    glTranslatef(0.0f, .0f, camera_distance);  // Move the scene away from the camera    
    glRotatef(-90.0f, 1.0f, 0.0f, 0.0f);  // Rotate around the x-axis
    

    // Apply rotation for pitch and yaw
    // Rotate around the x-axis for pitch
    glRotatef(pitch*rad2deg, 1.0f, 0.0f, 0.0f);
    // Rotate around the y-axis for yaw
    glRotatef(yaw*rad2deg, 0.0f, 0.0f, 1.0f);
    
    camera_pitch = pitch;
    camera_yaw = yaw;

    camera_normal = Float3{
    sin(yaw) * cos(pitch),
    cos(yaw) * cos(pitch),
    -sin(pitch)
    };// .norm();
}

Display::Display(EnvMode envmode) :
    logger(LimaLogger::LogMode::compact, envmode, "display") 
{    
    int success = initGLFW();
    logger.finishSection("Display initialized");

    setupProjection(display_width, display_height);
    updateCamera(0.f, 0.f);
    //updateCamera(PI/2.F, 0.f);    // Tilt the camera slightly up

    glfwSetWindowUserPointer(window, this);

    auto keyCallback = [](GLFWwindow* window, int key, int scancode, int action, int mods) {
        if (action == GLFW_PRESS) {
            // Retrieve the Display instance from the window user pointer
            Display* display = static_cast<Display*>(glfwGetWindowUserPointer(window));
            if (display) {
                const float delta = 3.1415 / 4.f;
                switch (key) {
                case GLFW_KEY_UP:
                    display->updateCamera(display->camera_pitch + delta, display->camera_yaw);
                    break;
                case GLFW_KEY_DOWN:
                    display->updateCamera(display->camera_pitch - delta, display->camera_yaw);
                    break;
                case GLFW_KEY_LEFT:
                    display->updateCamera(display->camera_pitch, display->camera_yaw + delta);
                    break;
                case GLFW_KEY_RIGHT:
                    display->updateCamera(display->camera_pitch, display->camera_yaw - delta);
                    break;
                case GLFW_KEY_PAGE_UP:
                    display->updateCamera(display->camera_pitch, display->camera_yaw, 0.5f);
                    break;
                case GLFW_KEY_PAGE_DOWN:
					display->updateCamera(display->camera_pitch, display->camera_yaw, -0.5f);
					break;
                }                              
            }
        }
        };

    glfwSetKeyCallback(window, keyCallback);
}

Display::~Display() {
    glfwTerminate();
}



void findPerpendicularVector(const Float3& v, Float3& perp1, Float3& perp2) {
    if (v.x != 0.f || v.y != 0.f) {
        perp1 = Float3{ -v.y, v.x, 0.f }.norm();
    }
    else {
        perp1 = Float3{ 0.f, -v.z, v.y }.norm();
    }
    perp2 = v.cross(perp1);
}

void Display::drawFilledCircle(const RenderAtom& ball) {
    float light = 0.5;
    Int3 shaded_color = ball.color * light;
    glColor3ub((uint8_t)shaded_color.x, (uint8_t)shaded_color.y, (uint8_t)shaded_color.z);

    glBegin(GL_TRIANGLE_FAN);
    glVertex3f(ball.pos.x, ball.pos.y, ball.pos.z);

    Float3 perp1, perp2;
    findPerpendicularVector(camera_normal, perp1, perp2);

    for (int i = 0; i <= triangleAmount; i++) {
        light = (sin(i * 2 * PI / triangleAmount) + 1.f) / 2.f;
        shaded_color = ball.color * light;

        glColor3ub((uint8_t)shaded_color.x, (uint8_t)shaded_color.y, (uint8_t)shaded_color.z);

        const float angle = i * twicePi / triangleAmount;

        // Calculate point on circle in the plane defined by perp1 and perp2
        const Float3 point = (perp1 * cos(angle) + perp2 * sin(angle)) * ball.radius;

        // Apply position offset
        const Float3 finalPoint = ball.pos + point;

        glVertex3f(finalPoint.x, finalPoint.y, finalPoint.z);
    }
    glEnd();
}


void Display::drawBalls(const std::vector<RenderAtom>& renderAtoms, int n_balls) {
    for (const auto& atom : renderAtoms) {
        if (!atom.Disabled()) {
            drawFilledCircle(atom);
        }
    }
}

void drawBoxOutline() {
    glBegin(GL_LINES);
    const float min = -0.5f;
    const float max = 0.5f;

    glLineWidth(20.0f);  // Set the line width

    // Bottom
    glColor3ub(50, 50, 200);

    glVertex3f(min, min, min);
    glVertex3f(max, min, min);

    glVertex3f(max, min, min);
    glVertex3f(max, max, min);

    glVertex3f(max, max, min);
    glVertex3f(min, max, min);

    glVertex3f(min, max, min);
    glVertex3f(min, min, min);

    glColor3ub(100, 100, 100);

    // Top
    glVertex3f(min, min, max);
    glVertex3f(max, min, max);

    glVertex3f(max, min, max);
    glVertex3f(max, max, max);

    glVertex3f(max, max, max);
    glVertex3f(min, max, max);

    glVertex3f(min, max, max);
    glVertex3f(min, min, max);

    // Sides
    glVertex3f(min, min, min);
    glVertex3f(min, min, max);

    glVertex3f(max, min, min);
    glVertex3f(max, min, max);

    glVertex3f(max, max, min);
    glVertex3f(max, max, max);

    glVertex3f(min, max, min);
    glVertex3f(min, max, max);
    glEnd();
}


void Display::render(const Float3* positions, const std::vector<Compound>& compounds, const BoxParams& boxparams, int64_t step, float temperature, ColoringMethod coloringMethod) {
    auto start = std::chrono::high_resolution_clock::now();

    const std::vector<RenderAtom>& renderAtoms = rasterizer.render(positions, compounds, boxparams, step, camera_normal, coloringMethod);
    //glClearColor(0x1c / 255.f, 0x24 / 255.f, 0x3f / 255.f, 1);
    //glClearColor(0.2f, 0.2f, 0.2f, 1.f);
    glClear(GL_COLOR_BUFFER_BIT);

    drawBoxOutline();
    //drawFilledCircle(balls[0]);
    drawBalls(renderAtoms, boxparams.total_particles_upperbound);

    // Swap front and back buffers
    glfwSwapBuffers(window);

    std::string window_text = std::format("{}        Step: {}    Temperature: {:.1f}[k]", window_title, step, temperature);
    glfwSetWindowTitle(window, window_text.c_str());

    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
    //printf("\tRender time: %4d ys  ", (int) duration.count());
}

bool Display::checkWindowStatus() {
    glfwPollEvents();
    if (glfwWindowShouldClose(window)) {
        glfwTerminate();

        return false;
    }

    return true;
}

bool Display::initGLFW() {
    
    logger.print("Initializing display...\n");

    // Initialize the library
    if (!glfwInit()) {
        throw std::exception("\nGLFW failed to initialize");
    }

    // Create a windowed mode window and its OpenGL context
    logger.print("Loading window --->");
    window = glfwCreateWindow(display_width, display_height, window_title.c_str(), NULL, NULL);
    if (!window)
    {
        glfwTerminate();
        return 0;
    }
    glfwSetWindowPos(window, screensize[0] - display_width - 550, 50);
    logger.print("done\n");

    // Make the window's context current
    glfwMakeContextCurrent(window);
    return 1;
}

#endif
