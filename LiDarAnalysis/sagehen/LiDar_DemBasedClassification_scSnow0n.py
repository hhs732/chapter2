import laspy as ls
import numpy as np
import pandas as pd

def readLasFile (lasFilePath):
    infile = ls.file.File(lasFilePath, mode="r")
    coords = np.vstack((infile.x, infile.y, infile.z)).T
    coords_df0 = pd.DataFrame(coords,columns=['x','y','z'])
    coords_df = pd.concat([coords_df0[['x','y']].astype(int),coords_df0['z']], axis=1)
    coords_df.sort_values(by=['x','y'],inplace=True)
    coords_df.index=np.arange(0,len(coords_df)) 
    
    return coords_df

coordSnow0nScM1 = readLasFile("C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/LiDarAnalysis/sagehen/las2016snow0n/26032016/WGS84_G1762_to_NAD83_NAVD88/USCASH20160326f2a1 - Channel 2 - 160326_223932_2 - originalpoints_dem_filter.las")
coordSnow0nScM2 = readLasFile("C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/LiDarAnalysis/sagehen/las2016snow0n/26032016/WGS84_G1762_to_NAD83_NAVD88/USCASH20160326f2a1 - Channel 1 - 160326_221954_1 - originalpoints_dem_filter.las")
coordSnow0nScM3 = readLasFile("C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/LiDarAnalysis/sagehen/las2016snow0n/26032016/WGS84_G1762_to_NAD83_NAVD88/USCASH20160326f2a1 - Channel 1 - 160326_222224_1 - originalpoints_dem_filter.las")
coordSnow0nScM4 = readLasFile("C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/LiDarAnalysis/sagehen/las2016snow0n/26032016/WGS84_G1762_to_NAD83_NAVD88/USCASH20160326f2a1 - Channel 1 - 160326_222504_1 - originalpoints_dem_filter.las")
coordSnow0nScM5 = readLasFile("C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/LiDarAnalysis/sagehen/las2016snow0n/26032016/WGS84_G1762_to_NAD83_NAVD88/USCASH20160326f2a1 - Channel 1 - 160326_222743_1 - originalpoints_dem_filter.las")
coordSnow0nScM6 = readLasFile("C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/LiDarAnalysis/sagehen/las2016snow0n/26032016/WGS84_G1762_to_NAD83_NAVD88/USCASH20160326f2a1 - Channel 1 - 160326_223059_1 - originalpoints_dem_filter.las")
coordSnow0nScM7 = readLasFile("C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/LiDarAnalysis/sagehen/las2016snow0n/26032016/WGS84_G1762_to_NAD83_NAVD88/USCASH20160326f2a1 - Channel 1 - 160326_223357_1 - originalpoints_dem_filter.las")
coordSnow0nScM8 = readLasFile("C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/LiDarAnalysis/sagehen/las2016snow0n/26032016/WGS84_G1762_to_NAD83_NAVD88/USCASH20160326f2a1 - Channel 1 - 160326_223635_1 - originalpoints_dem_filter.las")
coordSnow0nScM9 = readLasFile("C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/LiDarAnalysis/sagehen/las2016snow0n/26032016/WGS84_G1762_to_NAD83_NAVD88/USCASH20160326f2a1 - Channel 1 - 160326_223932_1 - originalpoints_dem_filter.las")
coordSnow0nScM10 = readLasFile("C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/LiDarAnalysis/sagehen/las2016snow0n/26032016/WGS84_G1762_to_NAD83_NAVD88/USCASH20160326f2a1 - Channel 2 - 160326_221955_2 - originalpoints_dem_filter.las")
coordSnow0nScM11 = readLasFile("C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/LiDarAnalysis/sagehen/las2016snow0n/26032016/WGS84_G1762_to_NAD83_NAVD88/USCASH20160326f2a1 - Channel 2 - 160326_222224_2 - originalpoints_dem_filter.las")
coordSnow0nScM12 = readLasFile("C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/LiDarAnalysis/sagehen/las2016snow0n/26032016/WGS84_G1762_to_NAD83_NAVD88/USCASH20160326f2a1 - Channel 2 - 160326_222504_2 - originalpoints_dem_filter.las")
coordSnow0nScM13 = readLasFile("C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/LiDarAnalysis/sagehen/las2016snow0n/26032016/WGS84_G1762_to_NAD83_NAVD88/USCASH20160326f2a1 - Channel 2 - 160326_222743_2 - originalpoints_dem_filter.las")
coordSnow0nScM14 = readLasFile("C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/LiDarAnalysis/sagehen/las2016snow0n/26032016/WGS84_G1762_to_NAD83_NAVD88/USCASH20160326f2a1 - Channel 2 - 160326_223059_2 - originalpoints_dem_filter.las")
coordSnow0nScM15 = readLasFile("C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/LiDarAnalysis/sagehen/las2016snow0n/26032016/WGS84_G1762_to_NAD83_NAVD88/USCASH20160326f2a1 - Channel 2 - 160326_223358_2 - originalpoints_dem_filter.las")
coordSnow0nScM16 = readLasFile("C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/LiDarAnalysis/sagehen/las2016snow0n/26032016/WGS84_G1762_to_NAD83_NAVD88/USCASH20160326f2a1 - Channel 2 - 160326_223635_2 - originalpoints_dem_filter.las")

coordsnow0nSc26M = pd.concat([coordSnow0nScM1,coordSnow0nScM2,coordSnow0nScM3,coordSnow0nScM4,coordSnow0nScM5,coordSnow0nScM6,
                              coordSnow0nScM7,coordSnow0nScM8,coordSnow0nScM9,coordSnow0nScM10,coordSnow0nScM11,coordSnow0nScM12,
                              coordSnow0nScM13,coordSnow0nScM14,coordSnow0nScM15,coordSnow0nScM16])

coordsnow0nSc26M.sort_values(by=['x','y'],inplace=True)
#%%
coordSnow0nSc1 = readLasFile("C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/LiDarAnalysis/sagehen/las2016snow0n/17042016/WGS84_G1762_to_NAD83_NAVD88/mcc_part_b_tile_000_002.las")
coordSnow0nSc2 = readLasFile("C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/LiDarAnalysis/sagehen/las2016snow0n/17042016/WGS84_G1762_to_NAD83_NAVD88/mcc_part_b_tile_000_003.las")
coordSnow0nSc3 = readLasFile("C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/LiDarAnalysis/sagehen/las2016snow0n/17042016/WGS84_G1762_to_NAD83_NAVD88/mcc_part_b_tile_000_004.las")
coordSnow0nSc4 = readLasFile("C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/LiDarAnalysis/sagehen/las2016snow0n/17042016/WGS84_G1762_to_NAD83_NAVD88/mcc_part_b_tile_000_005.las")
coordSnow0nSc5 = readLasFile("C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/LiDarAnalysis/sagehen/las2016snow0n/17042016/WGS84_G1762_to_NAD83_NAVD88/mcc_part_b_tile_001_000.las")
coordSnow0nSc6 = readLasFile("C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/LiDarAnalysis/sagehen/las2016snow0n/17042016/WGS84_G1762_to_NAD83_NAVD88/mcc_part_b_tile_001_001.las")
coordSnow0nSc7 = readLasFile("C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/LiDarAnalysis/sagehen/las2016snow0n/17042016/WGS84_G1762_to_NAD83_NAVD88/mcc_part_b_tile_001_002.las")
coordSnow0nSc8 = readLasFile("C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/LiDarAnalysis/sagehen/las2016snow0n/17042016/WGS84_G1762_to_NAD83_NAVD88/mcc_part_b_tile_001_003.las")
coordSnow0nSc9 = readLasFile("C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/LiDarAnalysis/sagehen/las2016snow0n/17042016/WGS84_G1762_to_NAD83_NAVD88/mcc_part_b_tile_001_004.las")
coordSnow0nSc10 = readLasFile("C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/LiDarAnalysis/sagehen/las2016snow0n/17042016/WGS84_G1762_to_NAD83_NAVD88/mcc_part_b_tile_001_005.las")
coordSnow0nSc11 = readLasFile("C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/LiDarAnalysis/sagehen/las2016snow0n/17042016/WGS84_G1762_to_NAD83_NAVD88/mcc_part_b_tile_002_000.las")
coordSnow0nSc12 = readLasFile("C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/LiDarAnalysis/sagehen/las2016snow0n/17042016/WGS84_G1762_to_NAD83_NAVD88/mcc_part_b_tile_002_001.las")
coordSnow0nSc13 = readLasFile("C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/LiDarAnalysis/sagehen/las2016snow0n/17042016/WGS84_G1762_to_NAD83_NAVD88/mcc_part_b_tile_002_002.las")
coordSnow0nSc14 = readLasFile("C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/LiDarAnalysis/sagehen/las2016snow0n/17042016/WGS84_G1762_to_NAD83_NAVD88/mcc_part_b_tile_002_003.las")
coordSnow0nSc15 = readLasFile("C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/LiDarAnalysis/sagehen/las2016snow0n/17042016/WGS84_G1762_to_NAD83_NAVD88/mcc_part_b_tile_002_004.las")
coordSnow0nSc16 = readLasFile("C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/LiDarAnalysis/sagehen/las2016snow0n/17042016/WGS84_G1762_to_NAD83_NAVD88/mcc_part_b_tile_002_005.las")
coordSnow0nSc17 = readLasFile("C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/LiDarAnalysis/sagehen/las2016snow0n/17042016/WGS84_G1762_to_NAD83_NAVD88/mcc_part_b_tile_003_000.las")
coordSnow0nSc18 = readLasFile("C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/LiDarAnalysis/sagehen/las2016snow0n/17042016/WGS84_G1762_to_NAD83_NAVD88/mcc_part_b_tile_003_001.las")
coordSnow0nSc19 = readLasFile("C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/LiDarAnalysis/sagehen/las2016snow0n/17042016/WGS84_G1762_to_NAD83_NAVD88/mcc_part_b_tile_003_002.las")
coordSnow0nSc20 = readLasFile("C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/LiDarAnalysis/sagehen/las2016snow0n/17042016/WGS84_G1762_to_NAD83_NAVD88/mcc_part_b_tile_003_003.las")
coordSnow0nSc21 = readLasFile("C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/LiDarAnalysis/sagehen/las2016snow0n/17042016/WGS84_G1762_to_NAD83_NAVD88/mcc_part_b_tile_003_004.las")
coordSnow0nSc22 = readLasFile("C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/LiDarAnalysis/sagehen/las2016snow0n/17042016/WGS84_G1762_to_NAD83_NAVD88/mcc_part_b_tile_003_005.las")
coordSnow0nSc23 = readLasFile("C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/LiDarAnalysis/sagehen/las2016snow0n/17042016/WGS84_G1762_to_NAD83_NAVD88/mcc_part_b_tile_004_000.las")
coordSnow0nSc24 = readLasFile("C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/LiDarAnalysis/sagehen/las2016snow0n/17042016/WGS84_G1762_to_NAD83_NAVD88/mcc_part_b_tile_004_001.las")
coordSnow0nSc25 = readLasFile("C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/LiDarAnalysis/sagehen/las2016snow0n/17042016/WGS84_G1762_to_NAD83_NAVD88/mcc_part_b_tile_004_002.las")
coordSnow0nSc26 = readLasFile("C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/LiDarAnalysis/sagehen/las2016snow0n/17042016/WGS84_G1762_to_NAD83_NAVD88/mcc_part_b_tile_004_003.las")
coordSnow0nSc27 = readLasFile("C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/LiDarAnalysis/sagehen/las2016snow0n/17042016/WGS84_G1762_to_NAD83_NAVD88/mcc_part_b_tile_004_004.las")
coordSnow0nSc28 = readLasFile("C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/LiDarAnalysis/sagehen/las2016snow0n/17042016/WGS84_G1762_to_NAD83_NAVD88/mcc_part_b_tile_004_005.las")
coordSnow0nSc29 = readLasFile("C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/LiDarAnalysis/sagehen/las2016snow0n/17042016/WGS84_G1762_to_NAD83_NAVD88/mcc_part_b_tile_005_000.las")
coordSnow0nSc30 = readLasFile("C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/LiDarAnalysis/sagehen/las2016snow0n/17042016/WGS84_G1762_to_NAD83_NAVD88/mcc_part_b_tile_005_001.las")
coordSnow0nSc31 = readLasFile("C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/LiDarAnalysis/sagehen/las2016snow0n/17042016/WGS84_G1762_to_NAD83_NAVD88/mcc_part_b_tile_005_002.las")
coordSnow0nSc32 = readLasFile("C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/LiDarAnalysis/sagehen/las2016snow0n/17042016/WGS84_G1762_to_NAD83_NAVD88/mcc_part_b_tile_005_003.las")
coordSnow0nSc33 = readLasFile("C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/LiDarAnalysis/sagehen/las2016snow0n/17042016/WGS84_G1762_to_NAD83_NAVD88/mcc_part_b_tile_005_004.las")
coordSnow0nSc34 = readLasFile("C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/LiDarAnalysis/sagehen/las2016snow0n/17042016/WGS84_G1762_to_NAD83_NAVD88/mcc_part_b_tile_005_005.las")
coordSnow0nSc35 = readLasFile("C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/LiDarAnalysis/sagehen/las2016snow0n/17042016/WGS84_G1762_to_NAD83_NAVD88/mcc_part_b_tile_000_000.las")
coordSnow0nSc36 = readLasFile("C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/LiDarAnalysis/sagehen/las2016snow0n/17042016/WGS84_G1762_to_NAD83_NAVD88/mcc_part_b_tile_000_001.las")

coordsnow0nSc17A = pd.concat([coordSnow0nSc1,coordSnow0nSc2,coordSnow0nSc3,coordSnow0nSc4,coordSnow0nSc5,coordSnow0nSc6,
                              coordSnow0nSc7,coordSnow0nSc8,coordSnow0nSc9,coordSnow0nSc10,coordSnow0nSc11,coordSnow0nSc12,
                              coordSnow0nSc13,coordSnow0nSc14,coordSnow0nSc15,coordSnow0nSc16,coordSnow0nSc17,coordSnow0nSc18,
                              coordSnow0nSc19,coordSnow0nSc20,coordSnow0nSc21,coordSnow0nSc22,coordSnow0nSc23,coordSnow0nSc24,
                              coordSnow0nSc25,coordSnow0nSc26,coordSnow0nSc27,coordSnow0nSc28,coordSnow0nSc29,coordSnow0nSc30,
                              coordSnow0nSc31,coordSnow0nSc32,coordSnow0nSc33,coordSnow0nSc34,coordSnow0nSc35,coordSnow0nSc36])

coordsnow0nSc17A.sort_values(by=['x','y'],inplace=True)
#%%
coordSnow0nScY1 = readLasFile("C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/LiDarAnalysis/sagehen/las2016snow0n/18052016/WGS84_G1762_to_NAD83_NAVD88/mcc_part_b_tile_002_003.las")
coordSnow0nScY2 = readLasFile("C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/LiDarAnalysis/sagehen/las2016snow0n/18052016/WGS84_G1762_to_NAD83_NAVD88/mcc_part_b_tile_002_004.las")
coordSnow0nScY3 = readLasFile("C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/LiDarAnalysis/sagehen/las2016snow0n/18052016/WGS84_G1762_to_NAD83_NAVD88/mcc_part_b_tile_002_005.las")
coordSnow0nScY4 = readLasFile("C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/LiDarAnalysis/sagehen/las2016snow0n/18052016/WGS84_G1762_to_NAD83_NAVD88/mcc_part_b_tile_003_000.las")
coordSnow0nScY5 = readLasFile("C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/LiDarAnalysis/sagehen/las2016snow0n/18052016/WGS84_G1762_to_NAD83_NAVD88/mcc_part_b_tile_003_001.las")
coordSnow0nScY6 = readLasFile("C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/LiDarAnalysis/sagehen/las2016snow0n/18052016/WGS84_G1762_to_NAD83_NAVD88/mcc_part_b_tile_003_002.las")
coordSnow0nScY7 = readLasFile("C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/LiDarAnalysis/sagehen/las2016snow0n/18052016/WGS84_G1762_to_NAD83_NAVD88/mcc_part_b_tile_003_003.las")
coordSnow0nScY8 = readLasFile("C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/LiDarAnalysis/sagehen/las2016snow0n/18052016/WGS84_G1762_to_NAD83_NAVD88/mcc_part_b_tile_003_004.las")
coordSnow0nScY9 = readLasFile("C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/LiDarAnalysis/sagehen/las2016snow0n/18052016/WGS84_G1762_to_NAD83_NAVD88/mcc_part_b_tile_003_005.las")
coordSnow0nScY10 = readLasFile("C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/LiDarAnalysis/sagehen/las2016snow0n/18052016/WGS84_G1762_to_NAD83_NAVD88/mcc_part_b_tile_004_000.las")
coordSnow0nScY11 = readLasFile("C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/LiDarAnalysis/sagehen/las2016snow0n/18052016/WGS84_G1762_to_NAD83_NAVD88/mcc_part_b_tile_004_001.las")
coordSnow0nScY12 = readLasFile("C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/LiDarAnalysis/sagehen/las2016snow0n/18052016/WGS84_G1762_to_NAD83_NAVD88/mcc_part_b_tile_004_002.las")
coordSnow0nScY13 = readLasFile("C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/LiDarAnalysis/sagehen/las2016snow0n/18052016/WGS84_G1762_to_NAD83_NAVD88/mcc_part_b_tile_004_003.las")
coordSnow0nScY14 = readLasFile("C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/LiDarAnalysis/sagehen/las2016snow0n/18052016/WGS84_G1762_to_NAD83_NAVD88/mcc_part_b_tile_004_004.las")
coordSnow0nScY15 = readLasFile("C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/LiDarAnalysis/sagehen/las2016snow0n/18052016/WGS84_G1762_to_NAD83_NAVD88/mcc_part_b_tile_004_005.las")
coordSnow0nScY16 = readLasFile("C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/LiDarAnalysis/sagehen/las2016snow0n/18052016/WGS84_G1762_to_NAD83_NAVD88/mcc_part_b_tile_005_000.las")
coordSnow0nScY17 = readLasFile("C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/LiDarAnalysis/sagehen/las2016snow0n/18052016/WGS84_G1762_to_NAD83_NAVD88/mcc_part_b_tile_005_001.las")
coordSnow0nScY18 = readLasFile("C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/LiDarAnalysis/sagehen/las2016snow0n/18052016/WGS84_G1762_to_NAD83_NAVD88/mcc_part_b_tile_005_002.las")
coordSnow0nScY19 = readLasFile("C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/LiDarAnalysis/sagehen/las2016snow0n/18052016/WGS84_G1762_to_NAD83_NAVD88/mcc_part_b_tile_005_003.las")
coordSnow0nScY20 = readLasFile("C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/LiDarAnalysis/sagehen/las2016snow0n/18052016/WGS84_G1762_to_NAD83_NAVD88/mcc_part_b_tile_005_004.las")
coordSnow0nScY21 = readLasFile("C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/LiDarAnalysis/sagehen/las2016snow0n/18052016/WGS84_G1762_to_NAD83_NAVD88/mcc_part_b_tile_005_005.las")
coordSnow0nScY22 = readLasFile("C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/LiDarAnalysis/sagehen/las2016snow0n/18052016/WGS84_G1762_to_NAD83_NAVD88/mcc_part_b_tile_000_000.las")
coordSnow0nScY23 = readLasFile("C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/LiDarAnalysis/sagehen/las2016snow0n/18052016/WGS84_G1762_to_NAD83_NAVD88/mcc_part_b_tile_000_001.las")
coordSnow0nScY24 = readLasFile("C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/LiDarAnalysis/sagehen/las2016snow0n/18052016/WGS84_G1762_to_NAD83_NAVD88/mcc_part_b_tile_000_002.las")
coordSnow0nScY25 = readLasFile("C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/LiDarAnalysis/sagehen/las2016snow0n/18052016/WGS84_G1762_to_NAD83_NAVD88/mcc_part_b_tile_000_003.las")
coordSnow0nScY26 = readLasFile("C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/LiDarAnalysis/sagehen/las2016snow0n/18052016/WGS84_G1762_to_NAD83_NAVD88/mcc_part_b_tile_000_004.las")
coordSnow0nScY27 = readLasFile("C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/LiDarAnalysis/sagehen/las2016snow0n/18052016/WGS84_G1762_to_NAD83_NAVD88/mcc_part_b_tile_000_005.las")
coordSnow0nScY28 = readLasFile("C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/LiDarAnalysis/sagehen/las2016snow0n/18052016/WGS84_G1762_to_NAD83_NAVD88/mcc_part_b_tile_001_000.las")
coordSnow0nScY29 = readLasFile("C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/LiDarAnalysis/sagehen/las2016snow0n/18052016/WGS84_G1762_to_NAD83_NAVD88/mcc_part_b_tile_001_001.las")
coordSnow0nScY30 = readLasFile("C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/LiDarAnalysis/sagehen/las2016snow0n/18052016/WGS84_G1762_to_NAD83_NAVD88/mcc_part_b_tile_001_002.las")
coordSnow0nScY31 = readLasFile("C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/LiDarAnalysis/sagehen/las2016snow0n/18052016/WGS84_G1762_to_NAD83_NAVD88/mcc_part_b_tile_001_003.las")
coordSnow0nScY32 = readLasFile("C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/LiDarAnalysis/sagehen/las2016snow0n/18052016/WGS84_G1762_to_NAD83_NAVD88/mcc_part_b_tile_001_004.las")
coordSnow0nScY33 = readLasFile("C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/LiDarAnalysis/sagehen/las2016snow0n/18052016/WGS84_G1762_to_NAD83_NAVD88/mcc_part_b_tile_001_005.las")
coordSnow0nScY34 = readLasFile("C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/LiDarAnalysis/sagehen/las2016snow0n/18052016/WGS84_G1762_to_NAD83_NAVD88/mcc_part_b_tile_002_000.las")
coordSnow0nScY35 = readLasFile("C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/LiDarAnalysis/sagehen/las2016snow0n/18052016/WGS84_G1762_to_NAD83_NAVD88/mcc_part_b_tile_002_001.las")
coordSnow0nScY36 = readLasFile("C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/LiDarAnalysis/sagehen/las2016snow0n/18052016/WGS84_G1762_to_NAD83_NAVD88/mcc_part_b_tile_002_002.las")

coordsnow0nSc18Y = pd.concat([coordSnow0nScY1,coordSnow0nScY2,coordSnow0nScY3,coordSnow0nScY4,coordSnow0nScY5,coordSnow0nScY6,
                              coordSnow0nScY7,coordSnow0nScY8,coordSnow0nScY9,coordSnow0nScY10,coordSnow0nScY11,coordSnow0nScY12,
                              coordSnow0nScY13,coordSnow0nScY14,coordSnow0nScY15,coordSnow0nScY16,coordSnow0nScY17,coordSnow0nScY18,
                              coordSnow0nScY19,coordSnow0nScY20,coordSnow0nScY21,coordSnow0nScY22,coordSnow0nScY23,coordSnow0nScY24,
                              coordSnow0nScY25,coordSnow0nScY26,coordSnow0nScY27,coordSnow0nScY28,coordSnow0nScY29,coordSnow0nScY30,
                              coordSnow0nScY31,coordSnow0nScY32,coordSnow0nScY33,coordSnow0nScY34,coordSnow0nScY35,coordSnow0nScY36])

coordsnow0nSc18Y.sort_values(by=['x','y'],inplace=True)




















