#-------------------------------------------------------------------------------
# Name:        donut_maker.py
# Purpose:     Convert footprint shapefiles to raster dataset
#
# Author:      moostang
#
# Created:     17/10/2018
# Copyright:   (c) moostang 2018
# Licence:     Free
# Version:     20181103
#-------------------------------------------------------------------------------
#
#   classDat   List of Raster Dataset in dirUnsClass
#     RasDat   Single instance of Raster Dataset from classDat
#weighted averaging
print "Running donut_make.py program"

# Import modules
print "Importing arcpy module and others"
import arcpy
import os
import datetime as dt

# Check if Spatial Analyst license is available
print "Checking license for Spatial Extension"
chk = arcpy.CheckOutExtension('Spatial')
print "Spatial Analyst is " + chk
del chk

# Do not show output to ArcMap (if code is executed in ArcMap)
# ------------------------------------------------------------
arcpy.env.addOutputsToMap = False

# ---------------------
# Directories and Files
# ---------------------
print "Reading input paramters"
dirBase     = "C:/GURUNG/MGIS-FINAL-PROJECT-WORK"
dirShp      = dirBase + "/DATA/SHAPEFILES"
dirCompRas  =  dirBase + "/DATA/IMAGERY/SENTINEL"
dirUnsClass = dirBase + "/DATA/CLASSIFICATION"
dirReclass  = dirBase + '/DATA/IMAGERY/FOOTPRINT'
#
ScratchGDB   = dirBase + "/temp/Cambay.gdb"
FootprintGDB = dirBase + "/temp/Footprint.gdb"
ReclassGDB   = dirBase + "/temp/Reclass.gdb"
DonutGDB     = dirBase + "/temp/Donut.gdb"
#
#ClipFile = dirShp + "/ClipRegion.shp"
ClipFile = dirBase + "/DATA/IMAGERY/clipregion"
#mxdFile  = dirBase + "/ARCGIS/Classification.mxd"
maskFile =  dirBase + "/DATA/IMAGERY/sea_mask.tif"
#
compFileExt = ".tif"
cellSize = 10
cellSzPr = 1
tStart = 190000 # HHMMSS

# Get Raster Dataset names
print "Getting list of SENTINEL-2 Composite images"
#files = []
#files = os.listdir(dirCompRas)
classDat = ['T13W_170622','T13W_170627','T13W_170629','T13W_170630','T13W_170702','T13W_170704','T13W_170705','T13W_170709']
##for compFile in files:
##    if compFile.endswith(compFileExt):
##        RasDat = compFile[0:4] + "_" + compFile[9:-11]
##        outRasDat = dirUnsClass + "/" + RasDat + ".tif"
##        classDat.append(RasDat)
##        if os.path.isfile(outRasDat) == False: #  or os.path.isfile(outRasFile) == True:
##            print "Unsupervised classification not yet done for " + compFile
##            print "Run Preprocessing program"
##        else:
##            print "Unsupervised classification for " + compFile + " Exists"
##
##del compFile, RasDat, outRasDat
#
# Create directories for final output
# -----------------------------------
for j in range(0,len(classDat)):
    path = dirReclass + '/' + classDat[j]
    if os.path.exists(path) == False:
        os.makedirs(path)
        print "Making directory " + path
    else:
        print "Output directorie are ready"

del path

# --------------------
# Open ArcMap Document (COMMENTED OUT)
# --------------------
# Define MXD Mapping Document and Dataframes
# mxd = arcpy.mapping.MapDocument(mxdFile)
## mxd = arcpy.mapping.MapDocument("CURRENT") # Change mxdFile to "CURRENT" when
                                             # scripting inside ArcMap.
#df = arcpy.mapping.ListDataFrames(mxd)[0]

# ------------------
# Set up Environment
# ------------------
print "Preparing Workspace Environment settings for Raster Analysis"
arcpy.env.workspace = dirUnsClass
# Define Extent
arcpy.env.extent = ClipFile
# Set cell size to .1 meters (to make smooth rasters)
arcpy.env.cellSize = cellSzPr
# Set the pyramid environment to build all pyramids levels with
#   cubic convolution resampling, LZ77 compression.
arcpy.env.pyramid = "PYRAMIDS -1 CUBIC LZ77 NO_SKIP"
# Statistics using a skip factor of 100 for x and y, and
# ignore values of 0 and 255.
arcpy.env.rasterStatistics = 'STATISTICS 100 100 (0 255)'
# Set the resampling method environment to bilinear interpolation.
arcpy.env.resamplingMethod = "BILINEAR"
#Set the tileSize environment to 128 by 128
arcpy.env.tileSize = "128 128"
# Set the nodata mapping method environment to NONE.
arcpy.env.nodata = "NONE"
# Set mask to exclude land areas
arcpy.env.mask = maskFile
# Overwrite TRUE
arcpy.env.overwriteOutput = True
# M Flag
arcpy.env.outputMFlag = "Disabled"
arcpy.env.outputZFlag = "Disabled"

# --------------------------
# Loop over Satellite Images
# --------------------------
for j in range(6,8):

    # ----------------------------------------
    # Read, convert, break footprint shapefile
    # ----------------------------------------
    # The footprint shapefiles resulting from the matlab program are polylines.
    # They have to be converted to polygon prior to raster processing.
    # Usually most SENTINEL images were taken at 1900 hours, i.e. 7PM
    # If the footprint at 1900 hours does not exist, look for another one that exist
    # before midnight, or 2340 hours.
    # If no footprint does not exists between tStart and 2340 hours, then search for
    # footprint before time tStart

    print "Reading RS image " + classDat[j]
    print "Get Footprints corresponding to +/- 12 hrs of of RS image"

    tYY = '20' + classDat[j][5: 7]
    tMM =        classDat[j][7: 9]
    tDD =        classDat[j][9:11]

    tBase = dt.datetime(int(tYY),int(tMM),int(tDD),6,40,0)
    tEnd  = tBase + dt.timedelta(hours= 24)

    del tYY, tMM, tDD

    while tBase < tEnd:
        tBase = tBase + dt.timedelta(minutes = 20)

        [date, time] = str(tBase).split()
        tagA = date.split('-')
        tagB = time.split(':')
        tagC = 'T'     + tagB[0] + tagB[1] + tagB[2]
        tag  = tagA[0] + tagA[1] + tagA[2] + tagC

        fcTarget = dirShp + '/FOOTPRINT/ALL/footprint_' + tag + '.shp'

        while os.path.isfile(fcTarget) == False and tBase < tEnd:
            tBase = tBase + dt.timedelta(minutes = 20)

            [date, time] = str(tBase).split()
            tagA = date.split('-')
            tagB = time.split(':')
            tagC = 'T'     + tagB[0] + tagB[1] + tagB[2]
            tag  = tagA[0] + tagA[1] + tagA[2] + tagC

            fcTarget = dirShp + '/FOOTPRINT/ALL/footprint_' + tag + '.shp'

        fc = fcTarget
        tagDT = tag

        # Sometimes program reads file at time tEnd, which may not exist.
        # This if clause will prevent it from reading it.
        if os.path.isfile(fc) == False:
            break
        else:
            print "Reading footprint footprint_" + tag + ".shp"
            print "--------------------------------------------------"

        del tagA, tagB, tag, fcTarget

        # -----------------------------
        # ArcMap: Show output in ArcMap (COMMENTED OUT)
        # -----------------------------
        #fcLyr = arcpy.mapping.Layer(fc)
        #arcpy.mapping.AddLayer(df, fcLyr, 'TOP')

        # -------------------------------------------
        # Loop over CONTOURS: 90%, 80%, ..., 20%, 10%
        # -------------------------------------------
        # Reset arcpy environment settings
        arcpy.env.workspace = dirUnsClass
        arcpy.env.mask = maskFile
        arcpy.env.extent = ClipFile
        arcpy.env.snapRaster = "NONE"

        tmpRas = []
        for i in range(1,10):
            ContourValue = str(i * cellSize) # ArcGIS only accept numeric as text

            fcOutPline = ScratchGDB + '/fcOutPline/fcOutPline' + ContourValue
            fcOutPoly = ScratchGDB + '/fcOutPoly/fcOutPoly' + ContourValue
            rsOutPoly = FootprintGDB + '/T' + classDat[j][5:] + 'p' + ContourValue
            tmpRas.append(rsOutPoly)

            # Execute Select
            where_clause = '"CONTOUR" = ' + ContourValue
            arcpy.Select_analysis(fc, fcOutPline, where_clause)
            print "Selecting contour " + ContourValue + "for " + fcOutPline

            # Several footprint shapefiles do not have the polylin for
            # 90% contribution. It exist in the table, BUT its geometry does
            # not exist. This causes a M-aware error, i.e. the geometry is
            # not M-aware. Since there is not polyline for 90% contribution,
            # this causes an error when converting an empty polyline to polygon.
            #
            try:
                # Convert polyline to polygon
                arcpy.FeatureToPolygon_management(fcOutPline,fcOutPoly,"#","ATTRIBUTES","#")
                print "Converting Polyline to Polygon for " + ContourValue

                # Change ID field value to match value of CONTOUR
                arcpy.AddField_management(fcOutPoly,'CONTOUR', "FLOAT", '2', "", "", "", "NULLABLE")
                with arcpy.da.UpdateCursor(fcOutPoly, 'CONTOUR') as cursor:
                    for row in cursor:
                        row[0] = i * cellSize
                        cursor.updateRow(row)

                # Convert polygon to raster
                arcpy.PolygonToRaster_conversion(fcOutPoly,'CONTOUR',rsOutPoly,"CELL_CENTER",'CONTOUR',cellSzPr)
                print "Converting Polygon to Raster for " + ContourValue

            except arcpy.ExecuteError:

                # Skip M-aware error and continue with computations
                arcpy.AddError(arcpy.GetMessages(2))
                print 'footprint_' + tagDT + " does not have contour value " + ContourValue
                continue


            del fcOutPline, fcOutPoly, rsOutPoly

        del fc, ContourValue, where_clause

        print "Finished converting polyline to raster"
        print "--------------------------------------"

        # ---------------------------------------
        # Make Doughnuts: 90%, 80%, ..., 20%, 10%
        # ---------------------------------------
        # Erase contribution raster with the next smaller contribution raster
        #. eg. Erase 90% with 80%.

        donut = []
        for i in range(8,0,-1):
            tag = str((i+1)*cellSize)

            # Set mask and extents: Use larger contribution
            arcpy.env.mask = tmpRas[i]
            arcpy.env.extent = tmpRas[i]
            arcpy.env.snapRaster = tmpRas[i]
            donought = DonutGDB + '/donut' + tag
            reclass = dirReclass + '/' + classDat[j] + '/' + classDat[j] +  tagC + '_' + tag + '.tif'
            donut.append(reclass)

            # Set Null Values to Zero
            exp = "Con(IsNull('" + tmpRas[i-1] + "'),0,'" + tmpRas[i-1] + "')"
            arcpy.gp.RasterCalculator_sa(exp,donought)
            print "Making Donut " + donought

            # Reclassify 0 to 1, and contour values to NoData
            arcpy.gp.Reclassify_sa(donought,"VALUE","0 1;"+str((i)*cellSize)+" NODATA",reclass,"NODATA")
            print "Reclassing " + donought + " into Reclass " +  tag

            del tag, donought, reclass

        # Donut for 10% contribution
        # --------------------------
        # Set mask and extents: Use 10%
        arcpy.env.mask = tmpRas[0]
        arcpy.env.extent = tmpRas[0]
        arcpy.env.snapRaster = tmpRas[0]

        donought = DonutGDB + "/donut10"
        reclass = dirReclass + '/' + classDat[j] + '/' + classDat[j] +  tagC + '_10.tif'
        donut.append(reclass)
        # Set Null Values to Zero
        exp = "Con(IsNull('" + tmpRas[0] + "'),0,'" + tmpRas[0] + "')"
        arcpy.gp.RasterCalculator_sa(exp,donought)
        # Reclassify 10% to 1. NoData values does not exist.
        arcpy.gp.Reclassify_sa(donought,"VALUE","10 1",reclass,"DATA")
        print "Reclassing " + donought + " into Reclass10"
        del donought, reclass

        print "----------------------"
        print "Finished making donuts"
        print "----------------------"

        # Create Fractions by multiplying with aggregate raster
        arcpy.env.mask = ''
        arcpy.env.snapRaster = ''

        AggDat = dirUnsClass + '/AGGREGATE/' + classDat[j] + '.tif'

        # Open textfile to write tabular data
        tFile = open(dirBase + '/temp/' + tagDT + '.csv', 'w')
        # https://community.esri.com/thread/84698

        for k in range(9,0,-1):
            ContourValue = str((10-k)*10) # Unfortunately the donut list starts
                                          # from 90, instead of 10
            print "Extract Sea-ice footprint"
##            arcpy.env.mask = tmpRas[k-1]
##            arcpy.env.extent = tmpRas[k-1]
##            arcpy.env.snapRaster = tmpRas[k-1]
            FptDat = (arcpy.sa.Raster(AggDat) * 1 )* ( arcpy.sa.Raster(donut[k-1]) *1)
            FptName = FootprintGDB + "/" + classDat[j] + "fpt" + ContourValue
            FptDat.save(FptName)

            # Create new field called "Contribution"
            arcpy.AddField_management(FptName, "Contribution", "SHORT", "", "", "", "", "NULLABLE")
            arcpy.CalculateField_management(FptName, field="Contribution", expression=ContourValue, expression_type="VB", code_block="")

            # Create new field "Shape_Area" and calculate area of shape
            print "Calculating Area for Contribution " + ContourValue + "%"
            arcpy.AddField_management(FptName, "Shape_Area", "FLOAT", '2', "", "", "", "NULLABLE")
            arcpy.CalculateField_management(FptName, field="Shape_Area", expression="[Count]*" + str(cellSzPr) + "", expression_type="VB", code_block="")

            # Create new field called "Percent"
            arcpy.AddField_management(FptName, "Percent", "FLOAT", "", "", "", "", "NULLABLE")

            # Create empty list to store the values of all aggregated fields
            Shape_List = []
            rows = arcpy.SearchCursor(FptName)
            for row in rows:
                area = row.getValue("Shape_Area")
                Shape_List.append(area)

            # Calculate total area of all aggregated fields
            Total_Area= sum(Shape_List)

            print "Calculating Percentage of Contribution for " + ContourValue + "% Contribution"
            arcpy.CalculateField_management(FptName, field="Percent", expression="([Shape_Area]/" + str(Total_Area) + ")*100", expression_type="VB", code_block="")

            # Write header information
            flds = arcpy.ListFields(FptName)
            header = ','.join(out.name for out in flds)
            if k == 9:
                tFile.write(header + '\n')

            rows = arcpy.SearchCursor(FptName) # Always read rows again to RESET
            for row in rows:
                lst  = [row.getValue(fld.name) for fld in flds]
                line = ','.join([str(a) for a in lst])
                tFile.write(line + '\n')

            del FptDat, Shape_List, rows, Total_Area

        tFile.close()

        del tFile, header, flds, FptName, AggDat

        # Delete files
        arcpy.env.workspace = ScratchGDB
        fc = arcpy.ListFeatureClasses()
        for file in fc:
            if arcpy.Exists(file):
                arcpy.Delete_management(file)






# test writing to gdb and deleting it.
# https://www.programcreek.com/python/example/107192/arcpy.Delete_management









