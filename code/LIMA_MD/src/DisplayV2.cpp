#include "DisplayV2.h"


DisplayV2::DisplayV2() {
#ifdef ENABLE_DISPLAY
    int success = initGLFW();
#endif
    printf("Display initialized\n");
}


void DisplayV2::drawFilledCircle(const RenderBall& ball) {
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
            ball.pos.y
        );
    }
    glEnd();
}


void DisplayV2::drawBalls(RenderBall* balls, int n_balls) {
    for (int i = n_balls-1; i >= 0; i--) {
        RenderBall& ball = balls[i];
        if (!ball.disable) { drawFilledCircle(ball); }
    }
}

void DisplayV2::render(Simulation* simulation) {
#ifdef ENABLE_DISPLAY
    auto start = std::chrono::high_resolution_clock::now();

    RenderBall* balls = rasterizer.render(simulation);
    glClear(GL_COLOR_BUFFER_BIT);


    drawBalls(balls, simulation->total_particles_upperbound);

    /* Swap front and back buffers */
    glfwSwapBuffers(window);

    delete[] balls;

    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
    //printf("\tRender time: %4d ys  ", (int) duration.count());
#endif
}

bool DisplayV2::checkWindowStatus() {
#ifdef ENABLE_DISPLAY
    glfwPollEvents();
    if (glfwWindowShouldClose(window)) {
        glfwTerminate();

        return false;
    }
#endif
    return true;
}

bool DisplayV2::initGLFW() {
    printf("Initializing display...  ");

    /* Initialize the library */
    if (!glfwInit()) {
        printf("GLFW failed to initialize\n");
        exit(0);
    }

    /* Create a windowed mode window and its OpenGL context */
    printf("Loading window --->");
    window = glfwCreateWindow(display_width, display_height, "LIMA - Molecular Dynamics Engine", NULL, NULL);
    if (!window)
    {
        glfwTerminate();
        return 0;
    }
    glfwSetWindowPos(window, screensize[0] - display_width, 0);
    printf("done\n");

    /* Make the window's context current */
    glfwMakeContextCurrent(window);
    return 1;
}