#include "TestUtils.h"
#include <random>
#include <format>

namespace InternalTypes {


    LimaUnittestResult TestCompression(EnvMode envmode) {
        constexpr int n = 1000;
        constexpr float maxDist = 0.13;
        constexpr float maxErrorOxygen = 1e-6f;
        constexpr float maxErrorHydrogen = 0.0021;

        std::mt19937 rng(1337); // fixed seed
        std::uniform_real_distribution<float> oDist(-2.68f, 2.68f);
        std::uniform_real_distribution<float> hDist(-maxDist, maxDist);

        for (int i = 0; i < n; ++i) {
            Float3 o(oDist(rng), oDist(rng), oDist(rng));

            Float3 h1 = o + Float3(hDist(rng), hDist(rng), hDist(rng));
            while ((h1 - o).len() > maxDist)
                h1 = o + Float3(hDist(rng), hDist(rng), hDist(rng));

            Float3 h2 = o + Float3(hDist(rng), hDist(rng), hDist(rng));
            while ((h2 - o).len() > maxDist)
                h2 = o + Float3(hDist(rng), hDist(rng), hDist(rng));

            CompressedSolvent packed(o, h1, h2);
            Float3 oUnpacked, h1Unpacked, h2Unpacked;
            packed.UnPack(oUnpacked, h1Unpacked, h2Unpacked);


			float oErr = (oUnpacked - o).len();
			float h1Err = (h1Unpacked - h1).len();
			float h2Err = (h2Unpacked - h2).len();

			if (oErr > maxErrorOxygen || h1Err > maxErrorHydrogen || h2Err > maxErrorHydrogen)
            {				
				return TestUtils::LimaUnittestResult(
					false,
					std::format("Failed to unpack: O dist {}, H1 dist {}, H2 dist {}",
						oErr, h1Err, h2Err
					),
					envmode == Full
				);
			}

            /*ASSERT((oUnpacked - o).len() < maxErrorOxygen, "Failed");
            ASSERT((h1Unpacked - h1).len() < maxErrorHydrogen, "Failed");
            ASSERT((h2Unpacked - h2).len() < maxErrorHydrogen, "Failed");*/
            ASSERT((oUnpacked - o).len() < maxErrorOxygen, std::format("O failed: dist {}", (oUnpacked - o).len()));
			ASSERT((h1Unpacked - h1).len() < maxErrorHydrogen, std::format("H1 failed: dist {}", (h1Unpacked - h1).len()));
			ASSERT((h2Unpacked - h2).len() < maxErrorHydrogen, std::format("H2 failed: dist {}", (h2Unpacked - h2).len()));



        }

        return LimaUnittestResult(true, "", envmode == Full);
    }
}