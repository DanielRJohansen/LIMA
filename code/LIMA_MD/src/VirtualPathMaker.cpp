#include "VirtualPathMaker.h"

Float3* VirtualPathMaker::makeVirtualPath(Float3* particle_positions, Float3* forces, int n_steps)
{
    Float3* virtual_path = new Float3[n_steps];


    Float3** all_paths = generateAllPaths(particle_positions, forces, n_steps);

    int best_path_index = findIndexOfShortestPath(all_paths, n_steps);
    for (int i = 0; i < n_steps; i++) {
        virtual_path[i] = all_paths[best_path_index][i] + particle_positions[i];
    }

    for (int i = 0; i < 8; i++)
        delete[] all_paths[i];
    delete[] all_paths;

    return virtual_path;
}

Float3** VirtualPathMaker::generateAllPaths(Float3* particle_positions, Float3* forces, int n_steps)
{
    Float3** all_paths = new Float3 * [8];
    for (int i = 0; i < 8; i++) {
        all_paths[i] = new Float3[n_steps];
    }

    Float3* initial_positions = generateAllPositions(particle_positions[0], forces[0]);
    for (int i = 0; i < 8; i++) {
        all_paths[i][0] = initial_positions[i];
    }
    delete[] initial_positions;


    for (int step = 1; step < n_steps; step++) {
        Float3* possible_positions = generateAllPositions(particle_positions[step], forces[step]);

        for (int j = 0; j < 8; j++) {
            int pos_index = findIndexOfShortestStep(all_paths[j][step - 1], possible_positions);
            all_paths[j][step] = possible_positions[pos_index];
        }
        delete[] possible_positions;
    }




    return all_paths;
}

int VirtualPathMaker::findIndexOfShortestStep(Float3 pos_prev, Float3* possible_positions)
{
    double shortest_step = 999999;
    int best_index = 0;

    for (int i = 0; i < 8; i++) {
        double dist = (pos_prev - possible_positions[i]).len();
        if (dist < shortest_step) {
            shortest_step = dist;
            best_index = i;
        }
    }
    if (shortest_step > 999)
        printf("FUCKIIIIING LONG STEP\n");
    return best_index;
}

int VirtualPathMaker::findIndexOfShortestPath(Float3** all_paths, int n_steps)
{
    double shortest_path = 1e+10;
    int best_index = 0;

    for (int i = 0; i < 8; i++) {
        double dist = lengthOfPath(all_paths[i], n_steps);
        //printf("Dist: %f\n", dist);
        if (dist < shortest_path) {
            shortest_path = dist;
            best_index = i;
        }
    }
    printf("Shortest path %f\tindex: %d", shortest_path, best_index);
    return best_index;
}



Float3* VirtualPathMaker::generateAllPositions(Float3 particle_position, Float3 force)
{
    Float3* possible_positions = new Float3[8];

    float buffer[6];

    for (int dim = 0; dim < 3; dim++) {
        float pos_repellent = static_cast<float>(getRepulsivePosition(force.at(dim)));
        float pos_attractive = static_cast<float>(getAtttractivePosition(force.at(dim)));
        buffer[dim * 2 + 0] = pos_repellent;
        buffer[dim * 2 + 1] = pos_attractive;
    }

    int cnt = 0;
    for (int x = 0; x < 2; x++) {
        for (int y = 0; y < 2; y++) {
            for (int z = 0; z < 2; z++) {
                possible_positions[cnt++] = Float3(buffer[x], buffer[y + 2], buffer[z + 4]);
            }
        }
    }

    //for (int i = 0; i < 8; i++)
      //  possible_positions[i].print();

    return possible_positions;
}

double VirtualPathMaker::getAtttractivePosition(double force)
{
    double lower = pow(2, 1.f / 6.f) * sigma;
    double upper = 1.2445 * sigma;
    double m = (lower + upper) / 2.f;
    if (abs(calcForce(upper)) < (abs(force))) { // If we cannot attract enough, we make the position very unlike to be chosen, by making it very expensive to choose this path.
        return force > 0 ? 9999 : -9999;
    }

    double pos_right_of_particle = binarySearch(lower, m, upper, abs(force));
    double pos = force < 0 ? pos_right_of_particle * -1 : pos_right_of_particle;


    return pos;
}

double VirtualPathMaker::getRepulsivePosition(double force)
{
    //printf("force: %f\n", force);
    double lower = pow(2, 1.f / 6.f) * sigma;
    double upper = 0.1;
    double m = (lower + upper) / 2.f;
    //printf("m: %f\n", m);
    double pos_right_of_particle = binarySearch(lower, m, upper, abs(force));
    double pos = force < 0 ? pos_right_of_particle : pos_right_of_particle * -1;

    return pos;
}



double VirtualPathMaker::binarySearch(double lower, double middle, double upper, double force)
{
    double dif = abs(calcForce(middle)) - force;
    //printf("m %f    dif: %f\n", middle, dif);

    if (abs(dif) < threshold)
        return middle;
    else if (dif > 0)
        return binarySearch(lower, (lower + middle) / 2.f, middle, force);
    else
        return binarySearch(middle, (middle + upper) / 2.f, upper, force);
}

double VirtualPathMaker::calcForce(double dist)
{
    double fraction = sigma / dist;		//nm/nm, so unitless

    double f2 = fraction * fraction;
    double f6 = f2 * f2 * f2;


    double force = -24.f * (epsilon / (dist * dist)) * f6 * (1.f - 2.f * f6);
    return force;
}

double VirtualPathMaker::lengthOfPath(Float3* path, int n_steps)
{
    double length = 0;
    for (int i = 1; i < n_steps; i++) {
        length += (path[i] - path[i - 1]).len();
    }
    return length;
}
