#include "LIMA_MD/include/DisplayV2.h"


Display::Display() :
    logger(LimaLogger::LogMode::compact, EnvMode::Full, "display") 
{    
#ifdef ENABLE_DISPLAY
    int success = initGLFW();
#endif
    logger.finishSection("Display initialized");
}

Display::~Display() {
#ifdef ENABLE_DISPLAY
    glfwTerminate();
#endif
}

void Display::drawFilledCircle(const RenderBall& ball) {
    float light = 0.5;
    Int3 shaded_color = ball.color * light;
    glColor3ub((uint8_t)shaded_color.x, (uint8_t)shaded_color.y, (uint8_t)shaded_color.z);

    glBegin(GL_TRIANGLE_FAN);
    glVertex3f(ball.pos.x, ball.pos.z, ball.pos.y + 1.f); // center of circle

    for (int i = 0; i <= triangleAmount; i++) {
        light = (sin(i * 2 * PI / triangleAmount) + 1.f) / 2.f;
        shaded_color = ball.color * light;

        glColor3ub((uint8_t)shaded_color.x, (uint8_t)shaded_color.y, (uint8_t)shaded_color.z);
        glVertex3f(
            ball.pos.x + (ball.radius * cos(i * twicePi / triangleAmount)),
            ball.pos.z + (ball.radius * sin(i * twicePi / triangleAmount)),
            -ball.pos.y
        );
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








void Display::render(Simulation* simulation) {
#ifdef ENABLE_DISPLAY
    auto start = std::chrono::high_resolution_clock::now();

    auto balls = rasterizer.render(simulation);
    glClear(GL_COLOR_BUFFER_BIT);

    drawBalls(balls, simulation->total_particles_upperbound);

    /* Swap front and back buffers */
    glfwSwapBuffers(window);

    std::string window_text = std::format("{}        Step: {}    Temperature: {:.1f}[k]", window_title, simulation->getStep(), simulation->temperature);
    glfwSetWindowTitle(window, window_text.c_str());

    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
    //printf("\tRender time: %4d ys  ", (int) duration.count());
#endif
}

bool Display::checkWindowStatus() {
#ifdef ENABLE_DISPLAY
    glfwPollEvents();
    if (glfwWindowShouldClose(window)) {
        glfwTerminate();

        return false;
    }
#endif
    return true;
}

bool Display::initGLFW() {
    
    logger.print("Initializing display...\n");

    /* Initialize the library */
    if (!glfwInit()) {
        printf("\nGLFW failed to initialize\n");
        exit(0);
    }

    /* Create a windowed mode window and its OpenGL context */
    logger.print("Loading window --->");
    window = glfwCreateWindow(display_width, display_height, window_title.c_str(), NULL, NULL);
    if (!window)
    {
        glfwTerminate();
        return 0;
    }
    glfwSetWindowPos(window, screensize[0] - display_width - 550, 50);
    logger.print("done\n");

    /* Make the window's context current */
    glfwMakeContextCurrent(window);
    return 1;
}