#pragma once

#include "LimaTypes.cuh"

#include <cstdio>
#include <vector>
#include <complex>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <cuda_runtime.h>
#include <cufft.h>




namespace PMEtest {

    const int gridpointsPerNm = 20;
    constexpr float gridpointsPerNm_f = static_cast<float>(gridpointsPerNm);
    constexpr float modifiedCoulombConstant = COULOMBCONSTANT / NANO / AVOGADROSNUMBER * KILO * KILO;	// [J/mol*nm / (kilo C/mol)^2]


    int GetGridIndexRealspace(const Int3& gridIndex, int gridpointsPerDim) {
        return gridIndex.x + gridIndex.y * gridpointsPerDim + gridIndex.z * gridpointsPerDim * gridpointsPerDim;
    }
    int GetGridIndexReciprocalspace(const Int3& gridIndex, int gridpointsPerDim, int nGridpointsHalfdim) {
        return gridIndex.x + gridIndex.y * nGridpointsHalfdim + gridIndex.z * nGridpointsHalfdim * gridpointsPerDim;
    }
    NodeIndex Get3dIndexRealspace(int index1d, int gridpointsPerDim) {
        int z = index1d / (gridpointsPerDim * gridpointsPerDim);
        index1d -= z * gridpointsPerDim * gridpointsPerDim;
        int y = index1d / gridpointsPerDim;
        index1d -= y * gridpointsPerDim;
        int x = index1d;
        return NodeIndex{ x, y, z };
    }
    NodeIndex Get3dIndexReciprocalspace(int index1d, int gridpointsPerDim, int nGridpointsHalfdim) {
        int z = index1d / (nGridpointsHalfdim * gridpointsPerDim);
        index1d -= z * nGridpointsHalfdim * gridpointsPerDim;
        int y = index1d / nGridpointsHalfdim;
        index1d -= y * nGridpointsHalfdim;
        int x = index1d;
        return NodeIndex{ x, y, z };
    }

    std::array<float, 4> CalcBspline(float f) {
        return {
            (1.f - f) * (1.f - f) * (1.f - f) / 6.f,
            (4.f - 6.f * f * f + 3.f * f * f * f) / 6.f,
            (1.f + 3.f * f + 3.f * f * f - 3.f * f * f * f) / 6.f,
            (f * f * f) / 6.f
        };
    }


    void DistributeChargesToGridKernel(const std::vector<Float3>& positions, const std::vector<float>& charges, float* chargeGrid, int gridpointsPerDim) {

        constexpr float delta = 1.f / gridpointsPerNm_f;
        constexpr float invCellVolume = 1.f / (delta * delta * delta);

        for (int pid = 0; pid < positions.size(); pid++) {
            const float charge = charges[pid];	        // [kC/mol]
            const Float3 absPos = positions[pid];       // [nm]

            // Map absolute position to fractional grid coordinates
            Float3 gridPos = absPos * gridpointsPerNm_f;
            int ix = static_cast<int>(floorf(gridPos.x));
            int iy = static_cast<int>(floorf(gridPos.y));
            int iz = static_cast<int>(floorf(gridPos.z));

            float fx = gridPos.x - static_cast<float>(ix);
            float fy = gridPos.y - static_cast<float>(iy);
            float fz = gridPos.z - static_cast<float>(iz);

            // Cubic B-spline weights along x
            float wx[4];
            wx[0] = (1.f - fx) * (1.f - fx) * (1.f - fx) / 6.f;
            wx[1] = (4.f - 6.f * fx * fx + 3.f * fx * fx * fx) / 6.f;
            wx[2] = (1.f + 3.f * fx + 3.f * fx * fx - 3.f * fx * fx * fx) / 6.f;
            wx[3] = (fx * fx * fx) / 6.f;

            // Cubic B-spline weights along y
            float wy[4];
            wy[0] = (1.f - fy) * (1.f - fy) * (1.f - fy) / 6.f;
            wy[1] = (4.f - 6.f * fy * fy + 3.f * fy * fy * fy) / 6.f;
            wy[2] = (1.f + 3.f * fy + 3.f * fy * fy - 3.f * fy * fy * fy) / 6.f;
            wy[3] = (fy * fy * fy) / 6.f;

            // Cubic B-spline weights along z
            float wz[4];
            wz[0] = (1.f - fz) * (1.f - fz) * (1.f - fz) / 6.f;
            wz[1] = (4.f - 6.f * fz * fz + 3.f * fz * fz * fz) / 6.f;
            wz[2] = (1.f + 3.f * fz + 3.f * fz * fz - 3.f * fz * fz * fz) / 6.f;
            wz[3] = (fz * fz * fz) / 6.f;


            // Distribute the charge to the surrounding 4x4x4 cells
            for (int dx = 0; dx < 4; dx++) {
                int X = -1 + ix + dx;
                float wxCur = wx[dx];
                for (int dy = 0; dy < 4; dy++) {
                    int Y = -1 + iy + dy;
                    float wxyCur = wxCur * wy[dy];
                    for (int dz = 0; dz < 4; dz++) {
                        int Z = -1 + iz + dz;
                        const NodeIndex index3d = NodeIndex{ X,Y,Z };
                        const int index1D = GetGridIndexRealspace(index3d, gridpointsPerDim);
                        chargeGrid[index1D] += charge * wxyCur * wz[dz] * invCellVolume; // TODO This must be applyBC
                    }
                }
            }
        }
    }





    static void InterpolateForcesAndPotentialKernel(
        const std::vector<Float3>& positions, const std::vector<float>& charges,
        const float* potentialGrid,
        int gridpointsPerDim,
        std::vector<ForceEnergy>& forceEnergies,
        float selfenergyCorrection			// [J/mol]
    )
    {
        for (int pid = 0; pid < positions.size(); pid++) {
            const float charge = charges[pid];	        // [kC/mol]
            const Float3 absPos = positions[pid];       // [nm]

            const Float3 gridPos = absPos * gridpointsPerNm_f;
            int ix = static_cast<int>(floorf(gridPos.x));
            int iy = static_cast<int>(floorf(gridPos.y));
            int iz = static_cast<int>(floorf(gridPos.z));

            float fx = gridPos.x - static_cast<float>(ix);
            float fy = gridPos.y - static_cast<float>(iy);
            float fz = gridPos.z - static_cast<float>(iz);

            float wx[4];
            wx[0] = (1.f - fx) * (1.f - fx) * (1.f - fx) / 6.f;
            wx[1] = (4.f - 6.f * fx * fx + 3.f * fx * fx * fx) / 6.f;
            wx[2] = (1.f + 3.f * fx + 3.f * fx * fx - 3.f * fx * fx * fx) / 6.f;
            wx[3] = (fx * fx * fx) / 6.f;

            float wy[4];
            wy[0] = (1.f - fy) * (1.f - fy) * (1.f - fy) / 6.f;
            wy[1] = (4.f - 6.f * fy * fy + 3.f * fy * fy * fy) / 6.f;
            wy[2] = (1.f + 3.f * fy + 3.f * fy * fy - 3.f * fy * fy * fy) / 6.f;
            wy[3] = (fy * fy * fy) / 6.f;

            float wz[4];
            wz[0] = (1.f - fz) * (1.f - fz) * (1.f - fz) / 6.f;
            wz[1] = (4.f - 6.f * fz * fz + 3.f * fz * fz * fz) / 6.f;
            wz[2] = (1.f + 3.f * fz + 3.f * fz * fz - 3.f * fz * fz * fz) / 6.f;
            wz[3] = (fz * fz * fz) / 6.f;

            Float3 force{};			// [J/mol/nm]
            float potential{};		// [J/mol]

            for (int dx = 0; dx < 4; dx++) {
                int X = ix - 1 + dx;
                float wxCur = wx[dx];
                for (int dy = 0; dy < 4; dy++) {
                    int Y = iy - 1 + dy;
                    float wxyCur = wxCur * wy[dy];
                    for (int dz = 0; dz < 4; dz++) {
                        int Z = iz - 1 + dz;
                        float wxyzCur = wxyCur * wz[dz];

                        const NodeIndex node = NodeIndex{ X, Y, Z };
                        const int gridIndex = GetGridIndexRealspace(node, gridpointsPerDim);

                        float phi = potentialGrid[gridIndex];

                        NodeIndex plusX = NodeIndex{ node.x + 1, node.y,     node.z };
                        NodeIndex minusX = NodeIndex{ node.x - 1, node.y,     node.z };
                        NodeIndex plusY = NodeIndex{ node.x,     node.y + 1, node.z };
                        NodeIndex minusY = NodeIndex{ node.x,     node.y - 1, node.z };
                        NodeIndex plusZ = NodeIndex{ node.x,     node.y,     node.z + 1 };
                        NodeIndex minusZ = NodeIndex{ node.x,     node.y,     node.z - 1 };

                        float phi_plusX = potentialGrid[GetGridIndexRealspace(plusX, gridpointsPerDim)];
                        float phi_minusX = potentialGrid[GetGridIndexRealspace(minusX, gridpointsPerDim)];
                        float phi_plusY = potentialGrid[GetGridIndexRealspace(plusY, gridpointsPerDim)];
                        float phi_minusY = potentialGrid[GetGridIndexRealspace(minusY, gridpointsPerDim)];
                        float phi_plusZ = potentialGrid[GetGridIndexRealspace(plusZ, gridpointsPerDim)];
                        float phi_minusZ = potentialGrid[GetGridIndexRealspace(minusZ, gridpointsPerDim)];

                        float E_x = -(phi_plusX - phi_minusX) * (gridpointsPerNm_f / 2.0f);
                        float E_y = -(phi_plusY - phi_minusY) * (gridpointsPerNm_f / 2.0f);
                        float E_z = -(phi_plusZ - phi_minusZ) * (gridpointsPerNm_f / 2.0f);

                        force += Float3{ E_x, E_y, E_z } *wxyzCur;
                        potential += phi * wxyzCur;
                    }
                }
            }

            // Now add self charge to calculations
            force *= charge;
            potential *= charge;

            // Ewald self-energy correction
            potential += selfenergyCorrection;

            potential *= 0.5f; // Potential is halved because we computing for both this and the other particle's

            forceEnergies[pid] = ForceEnergy{ force, potential };
        }
    }



    static void PrecomputeGreensFunctionKernel(float* d_greensFunction, int gridpointsPerDim,
        double boxLen,		// [nm]
        double ewaldKappa	// [nm^-1]
    ) {

        const int halfNodes = gridpointsPerDim / 2;
        int nGridpointsHalfdim = gridpointsPerDim / 2 + 1;

        for (int z = 0; z < gridpointsPerDim; z++) {
            for (int y = 0; y < gridpointsPerDim; y++) {
                for (int x = 0; x < nGridpointsHalfdim; x++) {

                    NodeIndex freqIndex{ x,y,z };

                    // Remap frequencies to negative for indices > N/2
                    int kxIndex = freqIndex.x;
                    int kyShiftedIndex = (freqIndex.y <= halfNodes) ? freqIndex.y : freqIndex.y - gridpointsPerDim;
                    int kzShiftedIndex = (freqIndex.z <= halfNodes) ? freqIndex.z : freqIndex.z - gridpointsPerDim;

                    // Ewald kappa fixed
                    double volume = boxLen * boxLen * boxLen;				// [nm^3]
                    double delta = boxLen / (double)gridpointsPerDim;		// [nm]

                    // Physical wavevectors
                    double kx = (2.0 * PI * (double)kxIndex) / boxLen;
                    double ky = (2.0 * PI * (double)kyShiftedIndex) / boxLen;
                    double kz = (2.0 * PI * (double)kzShiftedIndex) / boxLen;

                    double kSquared = kx * kx + ky * ky + kz * kz;

                    double currentGreensValue = 0.0f;

                    // Compute B-spline structure factor (4th order)
                    double kHalfX = kx * (delta * 0.5);
                    double kHalfY = ky * (delta * 0.5);
                    double kHalfZ = kz * (delta * 0.5);

                    const double epsilon = 1e-24;

                    auto splineFactor = [epsilon](double kh) {
                        if (fabs(kh) < epsilon) return 1.0;
                        double ratio = sin(kh) / kh;
                        return pow(ratio, 4);
                        };

                    double Sx = splineFactor(kHalfX);
                    double Sy = splineFactor(kHalfY);
                    double Sz = splineFactor(kHalfZ);

                    double splineCorrection = (Sx * Sy * Sz);
                    splineCorrection = splineCorrection * splineCorrection; // squared for forward+back interpolation

                    if (kSquared > epsilon) {
                        currentGreensValue = (4.0 * PI / (kSquared))
                            * exp(-kSquared / (4.0 * ewaldKappa * ewaldKappa))
                            * splineCorrection
                            * modifiedCoulombConstant
                            ;
                    }

                    const int index1D = GetGridIndexReciprocalspace(freqIndex, gridpointsPerDim, nGridpointsHalfdim);
                    d_greensFunction[index1D] = static_cast<float>(currentGreensValue);
                }
            }
        }
    }

    // Apply Green's function in reciprocal space (in-place)
    static void applyGreensFunction(
        cufftComplex* reciprocalGrad,
        const std::vector<float>& greens,
        int gridpointsPerDim)
    {
        const int halfNodes = gridpointsPerDim / 2;
        int nGridpointsHalfdim = gridpointsPerDim / 2 + 1;

        for (int i = 0; i < gridpointsPerDim * gridpointsPerDim * nGridpointsHalfdim; i++) {
            reciprocalGrad[i].x *= greens[i];
            reciprocalGrad[i].y *= greens[i];
        }
    }


    static float CalcEnergyCorrection(std::vector<float> charges, float ewaldKappa) {
        double chargeSquaredSum = 0;
        for (auto charge : charges)
                chargeSquaredSum += charge * charge;              
        return static_cast<float>(-ewaldKappa / sqrt(PI) * chargeSquaredSum * PhysicsUtils::modifiedCoulombConstant);
    }

    static Float3 ComputeDipoleMoment(const std::vector<Float3>& positions, const std::vector<float>& charges) {
        Float3 dipoleMoment = { 0.0f, 0.0f, 0.0f };
        for (size_t i = 0; i < charges.size(); ++i) {
            dipoleMoment += positions[i] * charges[i];
        }
        return dipoleMoment;
    }
    static float ComputeDipoleEnergyCorrection(const Float3& dipoleMoment, float boxVolume) {
        const float twoPiOverThree = 2.0f * PI / 3.0f;  // Constant factor
        return twoPiOverThree * dipoleMoment.lenSquared() / boxVolume;
    }
    static void ApplyDipoleForceCorrection(
        const Float3& dipoleMoment,
        const std::vector<float>& charges,
        std::vector<ForceEnergy>& forceEnergies,
        float boxVolume
    ) {
        const float dipoleForceFactor = (4.0f * PI) / (3.0f * boxVolume);

        for (size_t i = 0; i < charges.size(); ++i) {
            // Correction to the force due to the dipole
            float charge = charges[i];
            forceEnergies[i].force -= dipoleMoment * dipoleForceFactor * charge;
			forceEnergies[i].potE += ComputeDipoleEnergyCorrection(dipoleMoment, boxVolume) * 0.5f;
        }
    }

    // Full PME routine
    void computePME(
        const std::vector<Float3>& positions,
        const std::vector<float>& charges,
        float                      boxLenNm,
        float                      ewaldKappa,
        std::vector<ForceEnergy>& forceEnergy
        )
    {
		const int gridpointsPerDim = gridpointsPerNm * boxLenNm;
        size_t gridsizeRealspace = gridpointsPerDim * gridpointsPerDim * gridpointsPerDim;
		size_t gridsizeReciprocalspace = gridpointsPerDim * gridpointsPerDim * (gridpointsPerDim / 2 + 1);

        // Allocate host memory for charge grid
        std::vector<float> realspaceGrid_Host(gridsizeRealspace);

        // Precompute Green's function
        std::vector<float> greens(gridsizeReciprocalspace);
        PrecomputeGreensFunctionKernel(greens.data(), gridpointsPerDim, boxLenNm, ewaldKappa);

        // Spread charges
        DistributeChargesToGridKernel(positions, charges, realspaceGrid_Host.data(), gridpointsPerDim);

        // Allocate device memory
        cufftComplex* reciprocalGrid_Device;
        cudaMalloc(&reciprocalGrid_Device, sizeof(cufftComplex) * gridsizeReciprocalspace);
        float* realspaceGrid_Device;
		cudaMalloc(&realspaceGrid_Device, sizeof(float) * gridsizeRealspace);


        // Copy real charge grid to device as cufftComplex       
        cudaMemcpy(realspaceGrid_Device, realspaceGrid_Host.data(), sizeof(float) * gridsizeRealspace, cudaMemcpyHostToDevice);

        // Create plan and do forward FFT
        cufftHandle planForward;
        cufftPlan3d(&planForward, gridpointsPerDim, gridpointsPerDim, gridpointsPerDim, CUFFT_R2C);
        cufftExecR2C(planForward, realspaceGrid_Device, reciprocalGrid_Device);
        
        // Copy device data to host, apply Green's, copy back
		std::vector<cufftComplex> reciprocalGrid_Host(gridsizeReciprocalspace);
        cudaMemcpy(reciprocalGrid_Host.data(), reciprocalGrid_Device, sizeof(cufftComplex) * gridsizeReciprocalspace, cudaMemcpyDeviceToHost);
        applyGreensFunction(reciprocalGrid_Host.data(), greens, gridpointsPerDim);
        cudaMemcpy(reciprocalGrid_Device, reciprocalGrid_Host.data(), sizeof(cufftComplex) * gridsizeReciprocalspace, cudaMemcpyHostToDevice);

        // Inverse FFT
        cufftHandle planInverse;
        cufftPlan3d(&planInverse, gridpointsPerDim, gridpointsPerDim, gridpointsPerDim, CUFFT_C2R);
        cufftExecC2R(planInverse, reciprocalGrid_Device, realspaceGrid_Device);

        // Copy back to host real space
        cudaMemcpy(realspaceGrid_Host.data(), realspaceGrid_Device, sizeof(float) * gridsizeRealspace, cudaMemcpyDeviceToHost);

        for (auto& val : realspaceGrid_Host) 
			val /= (float)gridsizeRealspace;

        // Compute forces and energy
		InterpolateForcesAndPotentialKernel(positions, charges, realspaceGrid_Host.data(), gridpointsPerDim, forceEnergy, CalcEnergyCorrection(charges, ewaldKappa));

		//const Float3 dipoleMoment = ComputeDipoleMoment(positions, charges);
        //ApplyDipoleForceCorrection(dipoleMoment, charges, forceEnergy, boxLenNm * boxLenNm * boxLenNm);


    }
}