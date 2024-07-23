from MakeLipids import MakeLipids
import struct
import matplotlib.pyplot as plt
import numpy as np

def read_histogram_data(filename):
    with open(filename, 'rb') as f:
        # Read the number of steps
        numSteps = struct.unpack('i', f.read(4))[0]

        all_bins = []
        all_counts = []

        # Read bins size
        f.seek(0, 2)  # Move to the end of the file to get file size
        file_size = f.tell()
        f.seek(4)  # Move back to the first histogram data after the steps size
        remaining_size = file_size - 4  # Subtract the size of the steps_size integer

        for _ in range(numSteps):
            # Calculate number of bins based on remaining size
            num_bins = remaining_size // (8 + 4) // numSteps

            # Read bins
            bins = struct.unpack('q' * num_bins, f.read(8 * num_bins))
            all_bins.append(bins)

            # Read counts
            counts = struct.unpack('i' * num_bins, f.read(4 * num_bins))
            all_counts.append(counts)

    return all_bins, all_counts


def plot_histogram(all_bins, all_counts):
    nPlots = len(all_counts)
    fig, axes = plt.subplots(nPlots, 1)

    if nPlots == 1:
        axes = [axes]

    for ax, bins, counts in zip(axes, all_bins, all_counts):
        print(counts)
        x = np.arange(len(bins))
        ax.bar(x, counts, width=0.8, align='center', alpha=0.7)
        ax.set_yscale('log')  # Set the y-axis to a logarithmic scale
        ax.set_xticks(x)
        ax.set_xticklabels(bins, rotation=45)
        ax.set_xlabel('Bins')
        ax.set_ylabel('Frequency')

    plt.tight_layout()
    plt.show()


def plotStuff():
        kinE = [3.61565, 10.0414, 19.6751, 32.5108, 48.5405, 67.7545, 90.141, 115.686, 144.375, 176.189, 211.11,
                249.116, 290.184, 334.291, 381.409, 431.511, 484.567, 540.545, 599.415, 661.14, 725.686, 793.016,
                863.091, 935.872, 1011.32, 1089.39, 1170.04, 1253.23, 1338.91, 1427.03, 1517.56, 1610.44, 1705.62,
                1803.06, 1902.7, 2004.5, 2108.41, 2214.37, 2322.33, 2432.25, 2544.07, 2657.73, 2773.2, 2890.4, 3009.31,
                3129.85, 3251.98, 3375.64, 3500.79, 3627.38, 3755.34, 3884.64, 4015.21, 4147.02, 4280, 4414.12, 4549.31,
                4685.54, 4822.75, 4960.9, 5099.94, 5239.83, 5380.51, 5521.94, 5664.08, 5806.89, 5950.32, 6094.33,
                6238.88, 6383.92, 6529.42, 6675.34, 6821.64, 6968.28, 7115.22, 7262.44, 7409.88, 7557.53, 7705.34,
                7853.28, 8001.33, 8149.44, 8297.58, 8445.74, 8593.87, 8741.96, 8889.96, 9037.87, 9185.64, 9333.26,
                9480.7, 9627.93, 9774.94, 9921.7, 10068.2, 10214.4, 10360.3, 10505.8, 10651, 10795.9, 10940.3, 11084.4,
                11228, 11371.2, 11513.9, 11656.2, 11798, 11939.3, 12080.1, 12220.4, 12360.2, 12499.4, 12638.1, 12776.2,
                12913.7, 13050.7, 13187, 13322.8, 13458, 13592.5, 13726.5, 13859.8, 13992.4, 14124.5, 14255.9, 14386.6,
                14516.7, 14646.1, 14774.9, 14902.9, 15030.4, 15157.1, 15283.2, 15408.6, 15533.3, 15657.3, 15780.6,
                15903.3, 16025.2, 16146.5, 16267, 16386.9, 16506.1, 16624.6, 16742.4, 16859.5, 16975.9, 17091.6,
                17206.6, 17320.9, 17434.6, 17547.5, 17659.8, 17771.3, 17882.2, 17992.4, 18101.9, 18210.7, 18318.9,
                18426.4, 18533.2, 18639.3, 18744.7, 18849.5, 18953.6, 19057.1, 19159.9, 19262, 19363.5, 19464.3,
                19564.5, 19664, 19762.9, 19861.1, 19958.7, 20055.7, 20152.1, 20247.8, 20342.9, 20437.3, 20531.2,
                20624.4, 20717, 20809, 20900.5, 20991.3, 21081.5, 21171.1, 21260.1, 21348.6, 21436.4, 21523.7, 21610.4,
                21696.6, 21782.2, 21867.2, 21951.6, 22035.5, 22118.8, 22201.6, 22283.9, 22365.6, 22446.7, 22527.3,
                22607.4, 22687, 22766.1, 22844.6, 22922.6, 23000.1, 23077.1, 23153.6, 23229.5, 23305, 23380, 23454.5,
                23528.5, 23602, 23675.1, 23747.6, 23819.7, 23891.4, 23962.5, 24033.2, 24103.4, 24173.2, 24242.5,
                24311.4, 24379.9, 24447.8, 24515.4, 24582.5, 24649.2, 24715.4, 24781.3, 24846.7, 24911.7, 24976.2,
                25040.4, 25104.1, 25167.5, 25230.4, 25292.9, 25355.1, 25416.8, 25478.2, 25539.1, 25599.7, 25659.9,
                25719.7, 25779.2, 25838.2, 25896.9, 25955.3, 26013.2, 26070.8, 26128.1, 26185, 26241.5, 26297.7,
                26353.5, 26409, 26464.1, 26518.9, 26573.4, 26627.5, 26681.3, 26734.8, 26788, 26840.8, 26893.3, 26945.4,
                26997.3, 27048.9, 27100.1, 27151, 27201.6, 27251.9, 27301.9, 27351.7, 27401.1, 27450.2, 27499, 27547.5,
                27595.8, 27643.7, 27691.4, 27738.8, 27785.9, 27832.7, 27879.3, 27925.6, 27971.6, 28017.3, 28062.8,
                28108, 28152.9, 28197.6, 28242, 28286.2, 28330.1, 28373.8, 28417.2, 28460.4, 28503.3, 28545.9, 28588.3,
                28630.5, 28672.5, 28714.2, 28755.6, 28796.9, 28837.9, 28878.6, 28919.2, 28959.5, 28999.5, 29039.4,
                29079, 29118.4, 29157.6, 29196.6, 29235.4, 29273.9, 29312.2, 29350.4, 29388.3, 29426, 29463.5, 29500.8,
                29537.8, 29574.7, 29611.4, 29647.9, 29684.2, 29720.3, 29756.2, 29791.9, 29827.4, 29862.7, 29897.8,
                29932.8, 29967.5, 30002.1, 30036.5, 30070.7, 30104.7, 30138.6, 30172.2, 30205.7, 30239, 30272.2,
                30305.2, 30337.9, 30370.6, 30403, 30435.3, 30467.4, 30499.4, 30531.1, 30562.8, 30594.2, 30625.5,
                30656.7, 30687.6, 30718.5, 30749.1, 30779.6, 30810, 30840.2, 30870.2, 30900.1, 30929.8, 30959.4,
                30988.9, 31018.2, 31047.3, 31076.3, 31105.2, 31133.9, 31162.5, 31190.9, 31219.2, 31247.4, 31275.4,
                31303.3, 31331, 31358.6, 31386.1, 31413.4, 31440.6, 31467.7, 31494.6, 31521.4, 31548.1, 31574.7,
                31601.1, 31627.4, 31653.6, 31679.6, 31705.5, 31731.3, 31757, 31782.6, 31808, 31833.3, 31858.5, 31883.6,
                31908.5, 31933.4, 31958.1, 31982.7, 32007.2, 32031.6, 32055.8, 32060.8, 32046.6, 32032.4, 32018.1,
                32003.7, 31989.3, 31974.8, 31960.3, 31945.7, 31931.1, 31916.4, 31901.7, 31886.9, 31872.1, 31857.2,
                31842.2, 31827.2, 31812.2, 31797, 31781.9, 31766.6, 31751.3, 31736, 31720.6, 31705.1, 31689.6, 31674,
                31658.4, 31642.7, 31626.9, 31611.1, 31595.2, 31579.3, 31563.3, 31547.2, 31531.1, 31514.9, 31498.6,
                31482.3, 31465.9, 31449.5, 31433, 31416.4, 31399.8, 31383.1, 31366.4, 31349.5, 31332.6, 31315.7,
                31298.7, 31281.6, 31264.4, 31247.2, 31229.9, 31212.5, 31195.1, 31177.6, 31160, 31142.4, 31124.7,
                31106.9, 31089, 31071.1, 31053.1, 31035.1, 31016.9, 30998.7, 30980.4, 30962, 30943.6, 30925.1, 30906.5,
                30887.8, 30869.1, 30850.3, 30831.4, 30812.4, 30793.4, 30774.2, 30755, 30735.8, 30716.4, 30696.9,
                30677.4, 30657.8, 30638.1, 30618.4, 30598.5, 30578.6, 30558.5, 30538.4, 30518.2, 30498, 30477.6,
                30457.2, 30436.6, 30416, 30395.3, 30374.5, 30353.6, 30332.6, 30311.6, 30290.4, 30269.2, 30247.8,
                30226.4, 30204.9, 30183.3, 30161.6, 30139.8, 30117.9, 30095.9, 30073.8, 30051.6, 30029.3, 30006.9,
                29984.4, 29961.9, 29939.2, 29916.4, 29893.5, 29870.5, 29847.4, 29824.3, 29801, 29777.6, 29754.1,
                29730.5, 29706.8, 29682.9, 29659, 29635, 29610.8, 29586.6, 29562.2, 29537.7, 29513.2, 29488.5, 29463.7,
                29438.7, 29413.7, 29388.5, 29363.3, 29337.9, 29312.4, 29286.8, 29261, 29235.2, 29209.2, 29183.1,
                29156.9, 29130.5, 29104.1, 29077.5, 29050.8, 29023.9, 28996.9, 28969.8, 28942.6, 28915.3, 28887.8,
                28860.2, 28832.4, 28804.5, 28776.5, 28748.4, 28720.1, 28691.7, 28663.1, 28634.4, 28605.6, 28576.6,
                28547.5, 28518.2, 28488.8, 28459.3, 28429.6, 28399.8, 28369.8, 28339.7, 28309.4, 28279, 28248.4,
                28217.7, 28186.8, 28155.8, 28124.6, 28093.2, 28061.7, 28030.1, 27998.3, 27966.3, 27934.2, 27901.9,
                27869.4, 27836.8, 27804, 27771, 27737.9, 27704.6, 27671.2, 27637.5, 27603.7, 27569.7, 27535.6, 27501.3,
                27466.7, 27432.1, 27397.2, 27362.1, 27326.9, 27291.5, 27255.9, 27220.1, 27184.2, 27148, 27111.6,
                27075.1, 27038.4, 27001.4, 26964.3, 26927, 26889.5, 26851.8, 26813.9, 26775.7, 26737.4, 26698.9,
                26660.2, 26621.2, 26582.1, 26542.7, 26503.1, 26463.3, 26423.3, 26383.1, 26342.7, 26302, 26261.1, 26220,
                26178.7, 26137.1, 26095.3, 26053.3, 26011.1, 25968.6, 25925.9, 25882.9, 25839.7, 25796.3, 25752.7,
                25708.7, 25664.6, 25620.2, 25575.5, 25530.6, 25485.5, 25440.1, 25394.4, 25348.5, 25302.3, 25255.8,
                25209.1, 25162.2, 25114.9, 25067.4, 25019.6, 24971.6, 24923.3, 24874.7, 24825.8, 24776.6, 24727.2,
                24677.4, 24627.4, 24577.1, 24526.5, 24475.6, 24424.4, 24373, 24321.2, 24269.1, 24216.7, 24164, 24111,
                24057.7, 24004.1, 23950.1, 23895.9, 23841.3, 23786.4, 23731.2, 23675.6, 23619.7, 23563.5, 23507,
                23450.1, 23392.9, 23335.3, 23277.4, 23219.2, 23160.6, 23101.7, 23042.4, 22982.7, 22922.7, 22862.4,
                22801.6, 22740.5, 22679.1, 22617.2, 22555, 22492.5, 22429.5, 22366.2, 22302.5, 22238.3, 22173.9, 22109,
                22043.7, 21978, 21911.9, 21845.5, 21778.6, 21711.3, 21643.6, 21575.5, 21507, 21438.1, 21368.7, 21298.9,
                21228.7, 21158.1, 21087, 21015.5, 20943.6, 20871.2, 20798.4, 20725.1, 20651.4, 20577.3, 20502.7,
                20427.6, 20352.1, 20276.1, 20199.7, 20122.7, 20045.4, 19967.5, 19889.2, 19810.4, 19731.1, 19651.3,
                19571.1, 19490.3, 19409.1, 19327.4, 19245.2, 19162.5, 19079.3, 18995.5, 18911.3, 18826.6, 18741.4,
                18655.6, 18569.4, 18482.6, 18395.3, 18307.5, 18219.1, 18130.3, 18040.9, 17951, 17860.5, 17769.5, 17678,
                17586, 17493.4, 17400.3, 17306.6, 17212.4, 17117.6, 17022.3, 16926.5, 16830.1, 16733.1, 16635.6,
                16537.6, 16439, 16339.8, 16240.1, 16139.8, 16039, 15937.6, 15835.7, 15733.2, 15630.2, 15526.6, 15422.4,
                15317.7, 15212.4, 15106.6, 15000.2, 14893.3, 14785.8, 14677.8, 14569.2, 14460.1, 14350.5, 14240.2,
                14129.5, 14018.2, 13906.4, 13794, 13681.2, 13567.7, 13453.8, 13339.4, 13224.4, 13108.9, 12992.9,
                12876.4, 12759.4, 12642, 12524, 12405.5, 12286.6, 12167.2, 12047.4, 11927.1, 11806.3, 11685.1, 11563.5,
                11441.4, 11319, 11196.1, 11072.8, 10949.2, 10825.2, 10700.8, 10576, 10450.9, 10325.5, 10199.7, 10073.7,
                9947.32, 9820.69, 9693.79, 9566.64, 9439.24, 9311.63, 9183.81, 9055.79, 8927.61, 8799.26, 8670.78,
                8542.17, 8413.46, 8284.66, 8155.81, 8026.91, 7897.98, 7769.05, 7640.15, 7511.28, 7382.48, 7253.77,
                7125.17, 6996.7, 6868.4, 6740.29, 6612.38, 6484.72, 6357.32, 6230.22, 6103.44, 5977.01, 5850.96,
                5725.33, 5600.13, 5475.4, 5351.18, 5227.49, 5104.36, 4981.84, 4859.94, 4738.71, 4618.19, 4498.39,
                4379.36, 4261.14, 4143.76, 4027.25, 3911.65, 3797, 3683.34, 3570.7, 3459.11, 3348.62, 3239.27, 3131.09,
                3024.12, 2918.4, 2813.97, 2710.86, 2609.11, 2508.77, 2409.87, 2312.44, 2216.54, 2122.19, 2029.43,
                1938.3, 1848.84, 1761.09, 1675.08, 1590.84, 1508.43, 1427.86, 1349.18, 1272.43, 1197.63, 1124.81,
                1054.03, 985.291, 918.644, 854.116, 791.735, 731.533, 673.538, 617.778, 564.28, 513.07, 464.173,
                417.613, 373.413, 331.596, 292.182, 255.191, 220.642, 188.553, 158.939, 131.816, 107.199, 85.0985,
                65.5273, 48.4952, 34.0111, 22.0825, 12.7155, 5.91508, 1.68478, 0.0268042, 0.942012, 4.42993, 10.4887,
                19.1153, 30.305, 44.0522, 60.3496, 79.1888, 100.56, 124.452, 150.853, 179.749, 211.126, 244.966,
                281.254, 319.971, 361.096, 404.61, 450.491, 498.716, 549.261, 602.102, 657.212, 714.565, 774.132,
                835.887, 899.798, 965.837, 1033.97, 1104.17]
        potE = [41677.5, 41671.1, 41661.4, 41648.6, 41632.5, 41613.3, 41590.9, 41565.4, 41536.7, 41504.9, 41470, 41432,
                41390.9, 41346.8, 41299.7, 41249.6, 41196.5, 41140.5, 41081.7, 41019.9, 40955.4, 40888, 40818, 40745.2,
                40669.7, 40591.7, 40511, 40427.8, 40342.1, 40254, 40163.5, 40070.6, 39975.4, 39878, 39778.3, 39676.5,
                39572.6, 39466.7, 39358.7, 39248.8, 39136.9, 39023.3, 38907.8, 38790.6, 38671.7, 38551.1, 38429,
                38305.3, 38180.2, 38053.6, 37925.6, 37796.3, 37665.8, 37533.9, 37401, 37266.8, 37131.6, 36995.4,
                36858.2, 36720, 36581, 36441.1, 36300.4, 36159, 36016.8, 35874, 35730.6, 35586.6, 35442.1, 35297,
                35151.5, 35005.6, 34859.3, 34712.6, 34565.7, 34418.5, 34271, 34123.4, 33975.6, 33827.6, 33679.6,
                33531.5, 33383.3, 33235.2, 33087, 32938.9, 32790.9, 32643, 32495.2, 32347.6, 32200.2, 32053, 31905.9,
                31759.2, 31612.7, 31466.5, 31320.6, 31175, 31029.8, 30885, 30740.5, 30596.5, 30452.9, 30309.7, 30166.9,
                30024.6, 29882.8, 29741.5, 29600.7, 29460.4, 29320.7, 29181.5, 29042.8, 28904.7, 28767.1, 28630.2,
                28493.8, 28358, 28222.9, 28088.3, 27954.4, 27821.1, 27688.4, 27556.4, 27425, 27294.3, 27164.2, 27034.8,
                26906, 26777.9, 26650.5, 26523.8, 26397.7, 26272.3, 26147.6, 26023.6, 25900.3, 25777.6, 25655.7,
                25534.4, 25413.8, 25294, 25174.8, 25056.3, 24938.5, 24821.4, 24705, 24589.3, 24474.3, 24359.9, 24246.3,
                24133.4, 24021.1, 23909.6, 23798.7, 23688.5, 23579, 23470.1, 23362, 23254.5, 23147.7, 23041.6, 22936.2,
                22831.4, 22727.3, 22623.8, 22521, 22418.9, 22317.4, 22216.6, 22116.4, 22016.9, 21918, 21819.8, 21722.2,
                21625.2, 21528.8, 21433.1, 21338, 21243.6, 21149.7, 21056.5, 20963.9, 20871.9, 20780.4, 20689.6,
                20599.4, 20509.8, 20420.8, 20332.3, 20244.4, 20157.2, 20070.5, 19984.3, 19898.7, 19813.7, 19729.3,
                19645.4, 19562.1, 19479.3, 19397, 19315.3, 19234.2, 19153.6, 19073.5, 18993.9, 18914.8, 18836.3,
                18758.3, 18680.8, 18603.8, 18527.4, 18451.4, 18375.9, 18300.9, 18226.4, 18152.4, 18078.9, 18005.8,
                17933.3, 17861.2, 17789.6, 17718.4, 17647.7, 17577.5, 17507.7, 17438.4, 17369.5, 17301.1, 17233.1,
                17165.5, 17098.4, 17031.7, 16965.5, 16899.6, 16834.2, 16769.2, 16704.7, 16640.5, 16576.8, 16513.4,
                16450.5, 16388, 16325.8, 16264.1, 16202.7, 16141.8, 16081.2, 16021, 15961.2, 15901.7, 15842.7, 15784,
                15725.7, 15667.7, 15610.1, 15552.8, 15496, 15439.4, 15383.2, 15327.4, 15271.9, 15216.8, 15162, 15107.5,
                15053.4, 14999.6, 14946.1, 14893, 14840.1, 14787.7, 14735.5, 14683.6, 14632.1, 14580.8, 14529.9,
                14479.3, 14429, 14379, 14329.3, 14279.9, 14230.7, 14181.9, 14133.4, 14085.1, 14037.2, 13989.5, 13942.1,
                13895, 13848.2, 13801.6, 13755.4, 13709.3, 13663.6, 13618.1, 13572.9, 13528, 13483.3, 13438.9, 13394.7,
                13350.8, 13307.1, 13263.7, 13220.6, 13177.7, 13135, 13092.6, 13050.4, 13008.5, 12966.8, 12925.3,
                12884.1, 12843.1, 12802.3, 12761.8, 12721.5, 12681.4, 12641.5, 12601.9, 12562.5, 12523.3, 12484.3,
                12445.6, 12407, 12368.7, 12330.6, 12292.7, 12254.9, 12217.5, 12180.2, 12143.1, 12106.2, 12069.5, 12033,
                11996.7, 11960.6, 11924.7, 11889, 11853.5, 11818.2, 11783.1, 11748.1, 11713.4, 11678.8, 11644.4,
                11610.2, 11576.2, 11542.4, 11508.7, 11475.2, 11441.9, 11408.7, 11375.8, 11343, 11310.4, 11277.9,
                11245.6, 11213.5, 11181.6, 11149.8, 11118.2, 11086.7, 11055.4, 11024.3, 10993.3, 10962.5, 10931.8,
                10901.3, 10871, 10840.8, 10810.7, 10780.8, 10751.1, 10721.5, 10692.1, 10662.8, 10633.6, 10604.6,
                10575.7, 10547, 10518.4, 10490, 10461.7, 10433.6, 10405.6, 10377.7, 10349.9, 10322.3, 10294.9, 10267.5,
                10240.3, 10213.2, 10186.3, 10159.5, 10132.8, 10106.3, 10079.8, 10053.5, 10027.4, 10001.3, 9975.41,
                9949.61, 9923.94, 9898.39, 9872.96, 9847.64, 9822.45, 9797.38, 9772.42, 9747.58, 9722.86, 9698.26,
                9673.76, 9649.39, 9625.12, 7365.77, 7379.97, 7394.23, 7408.53, 7422.89, 7437.3, 7451.77, 7466.29,
                7480.86, 7495.48, 7510.16, 7524.9, 7539.68, 7554.52, 7569.42, 7584.38, 7599.38, 7614.45, 7629.57,
                7644.75, 7659.98, 7675.27, 7690.62, 7706.03, 7721.49, 7737.01, 7752.59, 7768.23, 7783.94, 7799.69,
                7815.51, 7831.39, 7847.33, 7863.34, 7879.4, 7895.52, 7911.71, 7927.96, 7944.27, 7960.64, 7977.08,
                7993.58, 8010.15, 8026.78, 8043.47, 8060.23, 8077.06, 8093.95, 8110.9, 8127.93, 8145.02, 8162.18,
                8179.4, 8196.7, 8214.06, 8231.49, 8248.99, 8266.56, 8284.2, 8301.91, 8319.69, 8337.54, 8355.47, 8373.46,
                8391.53, 8409.67, 8427.89, 8446.18, 8464.54, 8482.98, 8501.49, 8520.07, 8538.74, 8557.48, 8576.29,
                8595.19, 8614.16, 8633.21, 8652.33, 8671.54, 8690.83, 8710.19, 8729.63, 8749.16, 8768.77, 8788.46,
                8808.23, 8828.09, 8848.02, 8868.04, 8888.15, 8908.34, 8928.62, 8948.98, 8969.43, 8989.96, 9010.58,
                9031.29, 9052.09, 9072.97, 9093.95, 9115.02, 9136.17, 9157.42, 9178.76, 9200.19, 9221.71, 9243.33,
                9265.04, 9286.84, 9308.74, 9330.74, 9352.83, 9375.02, 9397.3, 9419.68, 9442.17, 9464.75, 9487.43,
                9510.21, 9533.09, 9556.07, 9579.16, 9602.34, 9625.63, 9649.03, 9672.53, 9696.13, 9719.85, 9743.66,
                9767.59, 9791.62, 9815.77, 9840.02, 9864.38, 9888.85, 9913.43, 9938.13, 9962.94, 9987.86, 10012.9,
                10038.1, 10063.3, 10088.7, 10114.2, 10139.8, 10165.6, 10191.4, 10217.4, 10243.5, 10269.7, 10296.1,
                10322.5, 10349.1, 10375.8, 10402.7, 10429.7, 10456.8, 10484, 10511.3, 10538.8, 10566.4, 10594.2,
                10622.1, 10650.1, 10678.2, 10706.5, 10734.9, 10763.5, 10792.2, 10821, 10850, 10879.1, 10908.4, 10937.8,
                10967.3, 10997, 11026.8, 11056.8, 11086.9, 11117.2, 11147.6, 11178.2, 11208.9, 11239.8, 11270.8, 11302,
                11333.4, 11364.9, 11396.5, 11428.3, 11460.3, 11492.4, 11524.7, 11557.2, 11589.8, 11622.6, 11655.6,
                11688.7, 11722, 11755.4, 11789.1, 11822.9, 11856.9, 11891, 11925.4, 11959.9, 11994.5, 12029.4, 12064.5,
                12099.7, 12135.1, 12170.7, 12206.5, 12242.5, 12278.6, 12315, 12351.5, 12388.2, 12425.2, 12462.3,
                12499.6, 12537.1, 12574.8, 12612.7, 12650.9, 12689.2, 12727.7, 12766.5, 12805.4, 12844.5, 12883.9,
                12923.5, 12963.3, 13003.3, 13043.5, 13083.9, 13124.6, 13165.5, 13206.6, 13247.9, 13289.5, 13331.3,
                13373.3, 13415.5, 13458, 13500.7, 13543.7, 13586.8, 13630.3, 13673.9, 13717.9, 13762, 13806.4, 13851.1,
                13896, 13941.1, 13986.5, 14032.2, 14078.1, 14124.3, 14170.8, 14217.5, 14264.4, 14311.7, 14359.2, 14407,
                14455, 14503.3, 14551.9, 14600.8, 14650, 14699.4, 14749.2, 14799.2, 14849.5, 14900.1, 14951, 15002.2,
                15053.6, 15105.4, 15157.5, 15209.9, 15262.6, 15315.6, 15368.9, 15422.5, 15476.5, 15530.7, 15585.3,
                15640.2, 15695.4, 15751, 15806.9, 15863.1, 15919.6, 15976.5, 16033.7, 16091.2, 16149.1, 16207.4, 16266,
                16324.9, 16384.2, 16443.9, 16503.9, 16564.2, 16625, 16686.1, 16747.5, 16809.3, 16871.5, 16934.1,
                16997.1, 17060.4, 17124.1, 17188.2, 17252.7, 17317.6, 17382.9, 17448.6, 17514.6, 17581.1, 17648,
                17715.3, 17783, 17851.1, 17919.6, 17988.5, 18057.9, 18127.7, 18197.9, 18268.5, 18339.6, 18411.1, 18483,
                18555.4, 18628.2, 18701.4, 18775.2, 18849.3, 18923.9, 18999, 19074.5, 19150.5, 19226.9, 19303.8,
                19381.2, 19459.1, 19537.4, 19616.2, 19695.5, 19775.3, 19855.5, 19936.2, 20017.5, 20099.2, 20181.4,
                20264.1, 20347.3, 20431, 20515.2, 20600, 20685.2, 20770.9, 20857.2, 20944, 21031.3, 21119.1, 21207.4,
                21296.3, 21385.7, 21475.6, 21566, 21657, 21748.5, 21840.6, 21933.2, 22026.3, 22120, 22214.2, 22308.9,
                22404.3, 22500.1, 22596.5, 22693.5, 22791, 22889, 22987.6, 23086.8, 23186.5, 23286.7, 23387.6, 23488.9,
                23590.9, 23693.4, 23796.4, 23900, 24004.2, 24108.9, 24214.1, 24319.9, 24426.3, 24533.2, 24640.7,
                24748.7, 24857.3, 24966.4, 25076.1, 25186.3, 25297.1, 25408.3, 25520.2, 25632.5, 25745.4, 25858.8,
                25972.7, 26087.2, 26202.2, 26317.6, 26433.6, 26550.1, 26667.1, 26784.6, 26902.6, 27021, 27139.9,
                27259.3, 27379.2, 27499.5, 27620.2, 27741.4, 27863.1, 27985.1, 28107.6, 28230.4, 28353.7, 28477.4,
                28601.4, 28725.8, 28850.6, 28975.7, 29101.1, 29226.8, 29352.9, 29479.2, 29605.9, 29732.8, 29859.9,
                29987.3, 30114.9, 30242.8, 30370.8, 30499, 30627.3, 30755.8, 30884.4, 31013.1, 31141.9, 31270.8,
                31399.7, 31528.6, 31657.5, 31786.4, 31915.3, 32044.1, 32172.8, 32301.4, 32429.9, 32558.2, 32686.3,
                32814.2, 32941.9, 33069.3, 33196.4, 33323.1, 33449.6, 33575.6, 33701.3, 33826.5, 33951.2, 34075.4,
                34199.1, 34322.2, 34444.8, 34566.7, 34687.9, 34808.4, 34928.2, 35047.2, 35165.5, 35282.9, 35399.4,
                35515, 35629.6, 35743.3, 35855.9, 35967.5, 36078, 36187.4, 36295.6, 36402.5, 36508.2, 36612.7, 36715.8,
                36817.5, 36917.9, 37016.8, 37114.2, 37210.1, 37304.5, 37397.2, 37488.4, 37577.8, 37665.6, 37751.6,
                37835.8, 37918.2, 37998.8, 38077.5, 38154.3, 38229.1, 38301.9, 38372.7, 38441.4, 38508.1, 38572.6,
                38635, 38695.2, 38753.2, 38808.9, 38862.4, 38913.6, 38962.5, 39009.1, 39053.3, 39095.1, 39134.5,
                39171.5, 39206.1, 39238.2, 39267.8, 39294.9, 39319.5, 39341.6, 39361.2, 39378.2, 39392.7, 39404.6,
                39414, 39420.8, 39425, 39426.7, 39425.8, 39422.3, 39416.2, 39407.6, 39396.4, 39382.7, 39366.4, 39347.5,
                39326.2, 39302.3, 39275.9, 39247, 39215.6, 39181.8, 39145.5, 39106.8, 39065.6, 39022.1, 38976.2, 38928,
                38877.5, 38824.6, 38769.5, 38712.1, 38652.6, 38590.8, 38526.9, 38460.9, 38392.7, 38322.5]
        totE = [41681.1, 41681.1, 41681.1, 41681.1, 41681.1, 41681.1, 41681.1, 41681.1, 41681.1, 41681.1, 41681.1,
                41681.1, 41681.1, 41681.1, 41681.1, 41681.1, 41681.1, 41681.1, 41681.1, 41681.1, 41681.1, 41681.1,
                41681.1, 41681.1, 41681.1, 41681.1, 41681, 41681, 41681, 41681, 41681, 41681, 41681, 41681, 41681,
                41681, 41681, 41681, 41681, 41681, 41681, 41681, 41681, 41681, 41681, 41681, 41681, 41681, 41681, 41681,
                41681, 41681, 41681, 41681, 41681, 41681, 41681, 41681, 41680.9, 41680.9, 41680.9, 41680.9, 41680.9,
                41680.9, 41680.9, 41680.9, 41680.9, 41680.9, 41680.9, 41680.9, 41680.9, 41680.9, 41680.9, 41680.9,
                41680.9, 41680.9, 41680.9, 41680.9, 41680.9, 41680.9, 41680.9, 41680.9, 41680.9, 41680.9, 41680.9,
                41680.9, 41680.9, 41680.9, 41680.9, 41680.9, 41680.9, 41680.9, 41680.9, 41680.9, 41680.9, 41680.9,
                41680.9, 41680.9, 41680.9, 41680.9, 41680.9, 41680.9, 41680.9, 41680.9, 41680.9, 41680.9, 41680.9,
                41680.9, 41680.9, 41680.9, 41680.9, 41680.9, 41680.9, 41680.9, 41680.9, 41680.9, 41680.9, 41680.9,
                41680.9, 41680.9, 41680.9, 41680.9, 41680.9, 41680.9, 41680.9, 41680.9, 41680.9, 41680.9, 41680.9,
                41680.9, 41680.9, 41680.9, 41680.9, 41680.9, 41680.9, 41680.9, 41680.9, 41680.9, 41680.9, 41680.9,
                41680.9, 41680.9, 41680.9, 41680.9, 41680.9, 41680.9, 41680.9, 41680.9, 41680.9, 41680.9, 41680.9,
                41680.9, 41680.9, 41680.9, 41680.9, 41680.9, 41680.9, 41680.9, 41680.9, 41680.9, 41680.9, 41680.9,
                41680.9, 41680.9, 41680.9, 41680.9, 41680.9, 41680.9, 41680.9, 41680.9, 41680.9, 41680.9, 41680.9,
                41680.9, 41680.9, 41680.9, 41680.9, 41680.9, 41680.9, 41680.9, 41680.9, 41680.9, 41680.9, 41680.9,
                41680.9, 41680.9, 41680.9, 41680.9, 41680.9, 41680.9, 41680.9, 41680.9, 41680.9, 41680.9, 41680.9,
                41680.9, 41680.9, 41680.9, 41680.9, 41680.9, 41680.9, 41680.9, 41680.9, 41680.9, 41680.9, 41680.9,
                41680.9, 41680.9, 41680.9, 41680.9, 41680.9, 41680.9, 41680.9, 41680.9, 41680.9, 41680.9, 41680.9,
                41680.9, 41680.9, 41680.9, 41680.9, 41680.9, 41680.9, 41680.9, 41680.9, 41680.9, 41680.9, 41680.9,
                41680.9, 41680.9, 41680.9, 41680.9, 41680.9, 41680.9, 41680.9, 41680.9, 41680.9, 41680.9, 41680.9,
                41680.9, 41680.9, 41680.9, 41680.9, 41680.9, 41680.9, 41680.9, 41680.9, 41680.9, 41680.9, 41680.9,
                41680.9, 41680.9, 41680.9, 41680.9, 41680.9, 41680.9, 41680.9, 41680.9, 41680.9, 41680.9, 41680.9,
                41680.9, 41680.9, 41680.9, 41680.9, 41680.9, 41680.9, 41680.9, 41680.9, 41680.9, 41680.9, 41680.9,
                41680.9, 41680.9, 41680.9, 41680.9, 41680.9, 41680.9, 41680.9, 41680.9, 41680.9, 41680.9, 41680.9,
                41680.9, 41680.9, 41680.9, 41680.9, 41680.9, 41680.9, 41680.9, 41680.9, 41680.9, 41680.9, 41680.9,
                41680.9, 41680.9, 41680.9, 41680.9, 41680.9, 41680.9, 41680.9, 41680.9, 41680.9, 41680.9, 41680.9,
                41680.9, 41680.9, 41680.9, 41680.9, 41680.9, 41680.9, 41680.9, 41680.9, 41680.9, 41680.9, 41680.9,
                41680.9, 41680.9, 41680.9, 41680.9, 41680.9, 41680.9, 41680.9, 41680.9, 41680.9, 41680.9, 41680.9,
                41680.9, 41680.9, 41680.9, 41680.9, 41680.9, 41680.9, 41680.9, 41680.9, 41680.9, 41680.9, 41680.9,
                41680.9, 41680.9, 41680.9, 41680.9, 41680.9, 41680.9, 41680.9, 41680.9, 41680.9, 41680.9, 41680.9,
                41680.9, 41680.9, 41680.9, 41680.9, 41680.9, 41680.9, 41680.9, 41680.9, 41680.9, 41680.9, 41680.9,
                41680.9, 41680.9, 41680.9, 41680.9, 41680.9, 41680.9, 41680.9, 41680.9, 41680.9, 41680.9, 41680.9,
                41680.9, 41680.9, 41680.9, 41680.9, 41680.9, 41680.9, 41680.9, 41680.9, 41680.9, 41680.9, 41680.9,
                41680.9, 41680.9, 41680.9, 41680.9, 41680.9, 41680.9, 41680.9, 41680.9, 41680.9, 41680.9, 41680.9,
                41680.9, 41680.9, 41680.9, 41680.9, 41680.9, 41680.9, 41680.9, 41680.9, 41680.9, 41680.9, 41681,
                41680.9, 41680.9, 41680.9, 41680.9, 41680.9, 41680.9, 41680.9, 41680.9, 41680.9, 41680.9, 41680.9,
                39426.6, 39426.6, 39426.6, 39426.6, 39426.6, 39426.6, 39426.6, 39426.6, 39426.6, 39426.6, 39426.6,
                39426.6, 39426.6, 39426.6, 39426.6, 39426.6, 39426.6, 39426.6, 39426.6, 39426.6, 39426.6, 39426.6,
                39426.6, 39426.6, 39426.6, 39426.6, 39426.6, 39426.6, 39426.6, 39426.6, 39426.6, 39426.6, 39426.6,
                39426.6, 39426.6, 39426.6, 39426.6, 39426.6, 39426.6, 39426.6, 39426.6, 39426.6, 39426.6, 39426.6,
                39426.6, 39426.6, 39426.6, 39426.6, 39426.6, 39426.6, 39426.6, 39426.6, 39426.6, 39426.6, 39426.6,
                39426.6, 39426.6, 39426.6, 39426.6, 39426.6, 39426.6, 39426.6, 39426.6, 39426.6, 39426.6, 39426.6,
                39426.6, 39426.6, 39426.6, 39426.6, 39426.6, 39426.6, 39426.6, 39426.6, 39426.6, 39426.6, 39426.6,
                39426.6, 39426.6, 39426.6, 39426.6, 39426.6, 39426.6, 39426.6, 39426.6, 39426.6, 39426.6, 39426.6,
                39426.6, 39426.6, 39426.6, 39426.6, 39426.6, 39426.6, 39426.6, 39426.6, 39426.6, 39426.6, 39426.6,
                39426.6, 39426.6, 39426.6, 39426.6, 39426.6, 39426.6, 39426.6, 39426.6, 39426.6, 39426.6, 39426.6,
                39426.6, 39426.6, 39426.6, 39426.6, 39426.6, 39426.6, 39426.6, 39426.6, 39426.6, 39426.6, 39426.6,
                39426.6, 39426.6, 39426.6, 39426.6, 39426.6, 39426.6, 39426.6, 39426.6, 39426.6, 39426.6, 39426.6,
                39426.6, 39426.6, 39426.6, 39426.6, 39426.6, 39426.6, 39426.6, 39426.6, 39426.6, 39426.6, 39426.6,
                39426.6, 39426.6, 39426.6, 39426.6, 39426.6, 39426.6, 39426.6, 39426.6, 39426.6, 39426.6, 39426.6,
                39426.6, 39426.6, 39426.6, 39426.6, 39426.6, 39426.6, 39426.6, 39426.6, 39426.6, 39426.6, 39426.6,
                39426.6, 39426.6, 39426.6, 39426.6, 39426.6, 39426.6, 39426.6, 39426.6, 39426.6, 39426.6, 39426.6,
                39426.6, 39426.6, 39426.6, 39426.6, 39426.6, 39426.6, 39426.6, 39426.6, 39426.6, 39426.6, 39426.6,
                39426.6, 39426.6, 39426.6, 39426.6, 39426.6, 39426.6, 39426.6, 39426.6, 39426.6, 39426.6, 39426.6,
                39426.6, 39426.6, 39426.6, 39426.6, 39426.6, 39426.6, 39426.6, 39426.6, 39426.6, 39426.6, 39426.6,
                39426.6, 39426.6, 39426.6, 39426.6, 39426.6, 39426.6, 39426.6, 39426.6, 39426.6, 39426.6, 39426.6,
                39426.6, 39426.6, 39426.6, 39426.6, 39426.6, 39426.6, 39426.6, 39426.6, 39426.6, 39426.6, 39426.6,
                39426.6, 39426.6, 39426.6, 39426.6, 39426.6, 39426.6, 39426.6, 39426.6, 39426.6, 39426.6, 39426.6,
                39426.6, 39426.6, 39426.6, 39426.6, 39426.6, 39426.6, 39426.6, 39426.6, 39426.6, 39426.6, 39426.6,
                39426.6, 39426.6, 39426.6, 39426.6, 39426.6, 39426.6, 39426.6, 39426.6, 39426.6, 39426.6, 39426.6,
                39426.6, 39426.6, 39426.6, 39426.6, 39426.6, 39426.6, 39426.6, 39426.6, 39426.6, 39426.6, 39426.6,
                39426.6, 39426.6, 39426.6, 39426.6, 39426.6, 39426.6, 39426.6, 39426.6, 39426.6, 39426.6, 39426.6,
                39426.6, 39426.6, 39426.6, 39426.6, 39426.6, 39426.6, 39426.6, 39426.6, 39426.6, 39426.6, 39426.6,
                39426.6, 39426.6, 39426.6, 39426.6, 39426.6, 39426.6, 39426.6, 39426.6, 39426.6, 39426.6, 39426.6,
                39426.6, 39426.6, 39426.6, 39426.6, 39426.6, 39426.6, 39426.6, 39426.6, 39426.6, 39426.6, 39426.6,
                39426.6, 39426.6, 39426.6, 39426.6, 39426.6, 39426.6, 39426.6, 39426.6, 39426.6, 39426.6, 39426.6,
                39426.6, 39426.6, 39426.6, 39426.6, 39426.6, 39426.6, 39426.6, 39426.6, 39426.6, 39426.6, 39426.6,
                39426.6, 39426.6, 39426.6, 39426.6, 39426.6, 39426.6, 39426.6, 39426.6, 39426.6, 39426.6, 39426.6,
                39426.6, 39426.6, 39426.6, 39426.6, 39426.6, 39426.6, 39426.6, 39426.6, 39426.6, 39426.6, 39426.6,
                39426.6, 39426.6, 39426.6, 39426.6, 39426.6, 39426.6, 39426.6, 39426.6, 39426.6, 39426.6, 39426.6,
                39426.6, 39426.6, 39426.6, 39426.6, 39426.6, 39426.6, 39426.6, 39426.6, 39426.6, 39426.6, 39426.6,
                39426.6, 39426.6, 39426.6, 39426.6, 39426.6, 39426.6, 39426.6, 39426.6, 39426.6, 39426.6, 39426.6,
                39426.6, 39426.6, 39426.6, 39426.6, 39426.6, 39426.6, 39426.6, 39426.6, 39426.6, 39426.6, 39426.6,
                39426.6, 39426.6, 39426.6, 39426.6, 39426.6, 39426.6, 39426.6, 39426.6, 39426.6, 39426.6, 39426.6,
                39426.6, 39426.6, 39426.6, 39426.6, 39426.6, 39426.6, 39426.6, 39426.6, 39426.6, 39426.6, 39426.6,
                39426.6, 39426.6, 39426.6, 39426.6, 39426.6, 39426.6, 39426.6, 39426.6, 39426.6, 39426.6, 39426.6,
                39426.6, 39426.6, 39426.6, 39426.6, 39426.6, 39426.6, 39426.6, 39426.6, 39426.6, 39426.6, 39426.6,
                39426.6, 39426.6, 39426.6, 39426.6, 39426.6, 39426.6, 39426.6, 39426.6, 39426.6, 39426.6, 39426.6,
                39426.6, 39426.6, 39426.6, 39426.6, 39426.6, 39426.6, 39426.6, 39426.6, 39426.6, 39426.6, 39426.6,
                39426.6, 39426.6, 39426.6, 39426.6, 39426.6, 39426.6, 39426.6, 39426.6, 39426.6, 39426.6, 39426.6,
                39426.6, 39426.6, 39426.6, 39426.6, 39426.6, 39426.6, 39426.6, 39426.6, 39426.6, 39426.6, 39426.6,
                39426.6, 39426.6, 39426.6, 39426.6, 39426.6, 39426.6, 39426.6, 39426.6, 39426.6, 39426.7, 39426.6,
                39426.7, 39426.7, 39426.7, 39426.7, 39426.7, 39426.7, 39426.7, 39426.7, 39426.7, 39426.7, 39426.7,
                39426.7, 39426.7, 39426.7, 39426.7, 39426.7, 39426.7, 39426.7, 39426.7, 39426.7, 39426.7, 39426.7,
                39426.7, 39426.7, 39426.7, 39426.7, 39426.7, 39426.7, 39426.7, 39426.7, 39426.7, 39426.7, 39426.7,
                39426.7, 39426.7, 39426.7, 39426.7, 39426.7, 39426.7, 39426.7, 39426.7, 39426.7, 39426.7, 39426.7,
                39426.7, 39426.7, 39426.7, 39426.7, 39426.7, 39426.7, 39426.7, 39426.7, 39426.7, 39426.7, 39426.7,
                39426.7, 39426.7, 39426.7, 39426.7, 39426.7, 39426.7, 39426.7, 39426.7, 39426.7, 39426.7, 39426.7,
                39426.7, 39426.7, 39426.7, 39426.7, 39426.7, 39426.7, 39426.7, 39426.7, 39426.7, 39426.7, 39426.7]

        plt.plot(kinE, label='Kinetic Energy')
        plt.plot(potE, label='Potential Energy')
        plt.plot(totE, label='Total Energy')
        plt.legend()
        plt.show()


if __name__ == "__main__":
    #MakeLipids

    #bins, counts = read_histogram_data('C:/Users/Daniel/git_repo/LIMA_data/psome2/histogram_data.bin')
    #plot_histogram(bins, counts)

    plt.plot(potE)
    plt.plot(kinE)
    plt.plot(totE)
    plt.show()


    #pot_energy = [38976.5, 33643.6, 26578.7, 18757.6, 11260.7, 5123.36, 1193.4, 13.6142, 1746.96, 6154.03, 12626.1, 20269.2, 28027.6, 34829.7, 39736, 42068.7, 41505.7, 38124.7, 32392.8, 25101.5, 17258.2, 9946.05, 4175.13, 742.553, 122.428, 2400.4, 7261.83, 14035.2, 21785, 29440.7, 35945, 40399.3, 42188.5, 41065.4, 37185, 31083.6, 23603.7, 15778.5, 8688.96, 3314.26, 396.775, 339.488, 3150.32, 8441.02, 15480.8, 23297.3, 30810.9, 36983.7, 40963.1, 42199.5, 40522, 36162.5, 29722.9, 22092.9, 14326.4, 7495.93, 2545.15, 157.815, 663.688, 3992.9, 9685.6, 16955.4, 24798.3, 32131, 37940.5, 41424.5, 42101.7, 39878.5, 35062.2, 28317.9, 20577.1, 12909.2, 6373.15, 1871.82, 26.9351, 1093.33, 4923.71, 10989, 18451.5, 26280.3, 33394.2, 38810.4, 41781, 41895.5, 39138.2, 33890, 26875.7, 19064.1, 11534.3, 5326.36, 1297.69, 4.79429, 1626.23, 5938.05, 12344.7, 19961.2, 27735.5, 34593.9]
    #plt.plot(pot_energy)
    #plt.show()


   # plotStuff()