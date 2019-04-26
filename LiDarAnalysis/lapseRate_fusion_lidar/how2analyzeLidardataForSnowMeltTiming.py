1) obtaining Lidar data snow-0n and snow-0ff and dem raster file 
2) tiling lidar data, if needed
3) generating allCoverCount (2m<), totalCount (allReturns), count (0.15m<) using FUSION
4) running "lidar_fusion_siteName.py" to calculate fsca for snow under tree and in open
5) executing "lapseRate6_elevNorthVegDClassification_allSites.py" to create final figures