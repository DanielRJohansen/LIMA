#include "DisplayV2.h"
//#include "glu"

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

void Display::updateCamera(float pitch, float yaw) {

    // Reset previous rotations
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    glTranslatef(0.0f, .0f, -3.5f);  // Move the scene away from the camera    
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

Display::Display() :
    logger(LimaLogger::LogMode::compact, EnvMode::Full, "display") 
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

void Display::drawFilledCircle(const RenderBall& ball) {
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


void Display::drawBalls(const std::vector<RenderBall>& balls, int n_balls) {
    for (const auto& ball : balls) {
        if (!ball.disable) {
            drawFilledCircle(ball);
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
    return;



    //glLineWidth(20.0f);  // Set the line width
    //glColor3f(1.0f, 0.0f, 0.0f);  // Set the color to red

    //glBegin(GL_LINES);

    //float boxWidth = 800;  // 80% of 1000 pixels
    //float boxHeight = 800; // Assuming square box for simplicity
    //float boxDepth = 800;  // Depth can be the same or different

    //Float3 minCorner = { -boxWidth / 2, -boxHeight / 2, -boxDepth / 2 };
    //Float3 maxCorner = { boxWidth / 2, boxHeight / 2, boxDepth / 2 };

    //Float3 points[8] = {
    //    {minCorner.at(0), minCorner.at(1), minCorner.at(2)},  // 0: min-min-min
    //    {maxCorner.at(0), minCorner.at(1), minCorner.at(2)},  // 1: max-min-min
    //    {maxCorner.at(0), minCorner.at(1), maxCorner.at(2)},  // 2: max-min-max
    //    {minCorner.at(0), minCorner.at(1), maxCorner.at(2)},  // 3: min-min-max
    //    {minCorner.at(0), maxCorner.at(1), minCorner.at(2)},  // 4: min-max-min
    //    {maxCorner.at(0), maxCorner.at(1), minCorner.at(2)},  // 5: max-max-min
    //    {maxCorner.at(0), maxCorner.at(1), maxCorner.at(2)},  // 6: max-max-max
    //    {minCorner.at(0), maxCorner.at(1), maxCorner.at(2)}   // 7: min-max-max
    //};

    //int edges[12][2] = {
    //    {0, 1}, {1, 2}, {2, 3}, {3, 0},  // Bottom edges
    //    {4, 5}, {5, 6}, {6, 7}, {7, 4},  // Top edges
    //    {0, 4}, {1, 5}, {2, 6}, {3, 7}   // Side edges
    //};

    //for (int i = 0; i < 12; ++i) {
    //    glVertex3f(points[edges[i][0]].at(0), points[edges[i][0]].at(1), points[edges[i][0]].at(2));
    //    glVertex3f(points[edges[i][1]].at(0), points[edges[i][1]].at(1), points[edges[i][1]].at(2));
    //}

    //glEnd();
}


void Display::render(const Float3* positions, const std::vector<Compound>& compounds, const BoxParams& boxparams, int64_t step, float temperature) {
    auto start = std::chrono::high_resolution_clock::now();

    auto balls = rasterizer.render(positions, compounds, boxparams, step, camera_normal);
    glClear(GL_COLOR_BUFFER_BIT);

    drawBoxOutline();
    //drawFilledCircle(balls[0]);
    drawBalls(balls, boxparams.total_particles_upperbound);

    // Swap front and back buffers
    glfwSwapBuffers(window);

    std::string window_text = std::format("{}        Step: {}    Temperature: {:.1f}[k]", window_title, step, temperature);
    glfwSetWindowTitle(window, window_text.c_str());

    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
    //printf("\tRender time: %4d ys  ", (int) duration.count());
}



//void keyCallback(GLFWwindow* window, int key, int scancode, int action, int mods) {
//    const float delta = PI / 4.f;
//
//    if (action == GLFW_PRESS) {
//        switch (key) {
//        case GLFW_KEY_UP:
//            updateCamera(camera_pitch - delta, camera_yaw);
//            break;
//        case GLFW_KEY_DOWN:
//            updateCamera(camera_pitch + delta, camera_yaw);
//            break;
//        case GLFW_KEY_LEFT:
//            updateCamera(camera_pitch, camera_yaw - delta);
//            break;
//        case GLFW_KEY_RIGHT:
//            updateCamera(camera_pitch, camera_yaw + delta);
//            break;
//        }
//    }
//}

bool Display::checkWindowStatus() {
    glfwPollEvents();
    if (glfwWindowShouldClose(window)) {
        glfwTerminate();

        return false;
    }


    //// Listen for arrow keys
    //const float delta = PI / 4.f;
    //if (glfwGetKey(window, GLFW_KEY_UP) == GLFW_PRESS) {
    //    updateCamera(camera_pitch - delta, camera_yaw);
    //}
    //if (glfwGetKey(window, GLFW_KEY_DOWN) == GLFW_PRESS) {
    //    updateCamera(camera_pitch + delta, camera_yaw);
    //}
    //if (glfwGetKey(window, GLFW_KEY_LEFT) == GLFW_PRESS) {
    //    updateCamera(camera_pitch, camera_yaw - delta);
    //}
    //if (glfwGetKey(window, GLFW_KEY_RIGHT) == GLFW_PRESS) {
    //    updateCamera(camera_pitch, camera_yaw + delta);
    //}





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

