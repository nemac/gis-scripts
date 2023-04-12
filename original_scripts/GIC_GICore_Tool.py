##=================================================================================
## Ian's Green Infrastructure Model, no part 4, ocean & background
## For ArcGIS 10.3 or greater
## Last Updated: September 2, 2020, 2:11pm
## This script combines various datasets available nationally or locally
## to identify and rank intact habitat cores, including forests, wetlands,
## aquatic habitat, and beach/dunes. Also includes background values to extend into deeper
## areas of the ocean. 
##
## The script accomplishes four major tasks:
##---------------------------------------------------------------------------------
##  Part 1: Setup
##  Scope: Make sure everything is set to run the model
##    1. Creates a scratch workspace for intermediate data if one does not exist.
##    2. Checks that the user has the necessary license level.
##    3. Checks that the user has the Spatial Analyst extension.
##    4. Checks the spatial reference information of the input layers, warning the
##       user if layers do not use the recommened projection system, and stopping
##      the script if any layers are unprojected.
##    5. Checks the shape type (point, line, polygon) of input layers.
##    6. If the user provides custom weights for ranking cores, the script
##       validates them.
##---------------------------------------------------------------------------------
##  Part 2: Identify Core Habitat Area
##  Scope: Identify natural areas based on land cover data, including forests,
##         wetlands, aquatic habitat, and beach/dunes.
##    1. Buffers the Area of Interest to user specified distance (if provided).
##    2. Clips the land cover and wetland data to the AOI.
##    3. Selects the following land cover types from the Cropland Data Layer or NLCD
##       ( NLCD values in brackets [] ):
##          111 [11]: OPEN WATER
##          141 [41]: DECIDUOUS FOREST
##          142 [42]: EVERGREEN FOREST
##          143 [43]: MIXED FOREST
##          190 [90]: WOODY WETLANDS
##          195 [95]: HERBACEOUS WETLANDS
##    4. Converts the wetlands data to raster
##    5. Uses NOAA C-Cap coastal land use data to identify beaches/dunes and
##       unconsolidated shore
##    5. Combines the layers to create a wetland, forest, and aquatic core habitat
##       raster
##---------------------------------------------------------------------------------
##  Part 3: Assess Fragmentation of Core Habitat
##  Scope: Combine layers of fragmenting features, and buffer them to identify edge
##         habitat in the Core Habitat raster from Part 2
##    1. Converts features that fragment the landscape into raster format, including:
##          * Paved Roads
##          * Buildings
##          * 'Developed' areas from Land Cover data
##          * Railroad tracks
##    2. Buffers fragmenting features.
##    3. Removes fragmenting features and their buffers (edge habitat) from core
##       habitat raster.
##    4. Identifies unique cores and removes cores less than the specified minimum
##       size.
##    5. Creates a separate raster for habitat fragments
##---------------------------------------------------------------------------------
##  Part 4: Calculate Core Metrics and Rank Cores
##  Scope: Calculate a variety of metrics for each core and use them to rank each
##         core based on a composite score.
##    1. Calculate metrics and add them as new fields in both a raster and a polygon
##       dataset.
##    2. Based on the metrics, calcualate the Core Quality Index and add to attribute
##       table.
##===================================================================================

import arcpy, os, sys, math, numpy, time
from datetime import datetime
from arcpy.sa import *
import gi_config as cfg

arcpy.env.overwriteOutput  = True

# Set workspace
out_workspace           = arcpy.GetParameterAsText(0)
arcpy.env.workspace     = out_workspace
# Collect parameters
area_of_interest        = arcpy.GetParameterAsText(1)
AOI_BUFFER_DISTANCE     = arcpy.GetParameter(2)      ## In miles. Optional.
LULC_SOURCE             = arcpy.GetParameterAsText(3) ## Cropland Data Layer or NLCD or C-CAP
lulc_layer              = arcpy.GetParameterAsText(4)
aoi_is_coastal          = arcpy.GetParameter(5) ## Boolean
noaa_ccap               = arcpy.GetParameterAsText(6)
ssurgo_soils            = arcpy.GetParameterAsText(7)
nwi_wetlands            = arcpy.GetParameterAsText(8)
nhd_flowlines           = arcpy.GetParameterAsText(9)
nhd_waterbodies         = arcpy.GetParameterAsText(10)
nhd_area                = arcpy.GetParameterAsText(11)
specrich_ras            = arcpy.GetParameterAsText(12)
ecosysred_ras             = arcpy.GetParameterAsText(13)
biodiv_ras              = arcpy.GetParameterAsText(14)
dem                     = arcpy.GetParameterAsText(15)
local_roads             = arcpy.GetParameterAsText(16)
railroads               = arcpy.GetParameterAsText(17)
building_locations      = arcpy.GetParameterAsText(18) ## Can accept point, polygon, or line features
smoothing_passes        = arcpy.GetParameter(19) ## Long
simplify_polygons       = arcpy.GetParameter(20) ## SIMPLIFY or NO_SIMPLIFY

# Globals
BUFFER_DISTANCE   = cfg.BUFFER_DISTANCE
MIN_CORE_SIZE     = cfg.MIN_CORE_SIZE
MIN_FRAGMENT_SIZE = cfg.MIN_FRAGMENT_SIZE
DEFAULT_WEIGHTS   = cfg.DEFAULT_WEIGHTS
OUTPUT_NAME_POLY  = cfg.OUTPUT_NAME_POLY
OUTPUT_NAME       = cfg.OUTPUT_NAME
HAB_FRAG_NAME     = cfg.HAB_FRAG_NAME
CORE_EDGE         = cfg.CORE_EDGE
lulc_aoi          = None ## global to use if a buffer is applied to the AOI
#-------------------------------------------------------------------------------------------------------------------------------------------#
# Helper functions -- Module #1
def checkShapeType(fc, expected_type):
    shapeType = arcpy.Describe(fc).shapeType
    if shapeType.upper() != expected_type.upper():
        errorhelp = (fc, expected_type, shapeType)
        return errorhelp
    else:
        return "PASS"

    raise CustomWeightsError #Stop Script

def EndModule(module_name):
    arcpy.AddMessage(("*"*20) + " Finished " + module_name + " " + ("*"*20) +"\n")

# Helper functions -- Module #2
def ScratchName(name):
    return os.path.join(scratchws, name)

def CheckBufferDistance(buffer):  ## Handle various cases of buffer distances
    global AOI_BUFFER_DISTANCE
    if buffer:
        if buffer == 0:
            arcpy.AddWarning("Hey dingbat, the Area of Interest buffer is set to 0. It is recommended to use a buffer for your study area.")
        if buffer > 50:
            arcpy.AddWarning("Hey dingbat, the Area of Interest buffer distance is greater than 50 miles. This is a very large buffer and will slow processing.")
    else:
        arcpy.AddWarning("Hey dingbat, the Area of Interest buffer distance was left blank. The AOI will not be buffered.")
        AOI_BUFFER_DISTANCE = 0

def WetlandsRasterField(wetlands): ## Find the field to use for feature to raster conversion. Input must be gdb. Prefer to use NWI 'ATTRIBUTE' field, use OID field otherwise.
    fields = arcpy.ListFields(wetlands)
    fieldnames = [field.name for field in fields]
    if ("ATTRIBUTE" in fieldnames):
        return "ATTRIBUTE"
    else:
        desc = arcpy.Describe(wetlands)
        if desc.hasOID:
            return desc.OIDFieldName
        else:
            arcpy.AddWarning("The wetlands layer must have a OID field.")
            raise WetlandsOIDError


# Helper Functions -- Module #3
def BuildSQL(data_path, value_list, where_field, operator):
    sql = arcpy.AddFieldDelimiters(data_path, where_field) + "='" + value_list.pop(0) + "'"
    if len(value_list):
        for value in value_list:
            sql = sql + " " + operator + " " + arcpy.AddFieldDelimiters(data_path, where_field) + "='" + value + "'"
    return sql


def FindMinCoreCells(core_size, data):
    if (arcpy.Describe(data).spatialReference.linearUnitName == "Meter"):
        meters_sq = core_size * 4046.87 ## convert the minimum core size from acres to square meters
        return meters_sq / math.pow(float(arcpy.env.cellSize), 2) ## divide the minimum core size (in sq. meters) by the cell size of the raster to obtain the number of cells that equals the minimum core size
    elif (arcpy.Describe(data).spatialReference.linearUnitName == "Foot"):
        feet_sq = core_size * 43560 ## convert the minimum core size from acres to square feet
        return feet_sq / math.pow(float(arcpy.env.cellSize), 2)
    else:
        arcpy.AddError("While analyzing the size of the cores, the script found that the input data is not using meters or feet as the linear unit of the projection. This will cause problems with area measurement.")
        raise LinearUnitError


def CreateCentroids(fc, out_path):
    multipart_fc = arcpy.MultipartToSinglepart_management(fc, os.path.join("in_memory", "multi"))
    in_count = int(arcpy.GetCount_management(fc).getOutput(0))
    test_count = int(arcpy.GetCount_management(multipart_fc).getOutput(0))
    x = fc if in_count >= test_count else multipart_fc
    sr = arcpy.Describe(fc).spatialReference
    centroid_coords = []
    with arcpy.da.SearchCursor(x, "SHAPE@TRUECENTROID") as cursor:
        for feature in cursor:
            centroid_coords.append(feature[0])
    point = arcpy.Point()
    pointGeometryList = []
    for pt in centroid_coords:
        if isinstance(pt[0], float):
            point.X = pt[0]
            point.Y = pt[1]
            pointGeometry = arcpy.PointGeometry(point, sr)
            pointGeometryList.append(pointGeometry)
    arcpy.Delete_management(os.path.join("in_memory", "multi"))
    return arcpy.CopyFeatures_management(pointGeometryList, out_path)


# Helper Functions -- Module #4
def ScratchName(name):
    return os.path.join(arcpy.env.scratchGDB, name)


def BuildSQL2(data_path, value_list, where_field, operator):
    sql = arcpy.AddFieldDelimiters(data_path, where_field) + "=" + str(value_list.pop(0))
    if len(value_list):
        for value in value_list:
            sql = sql + " " + operator + " " + arcpy.AddFieldDelimiters(data_path, where_field) + "=" + str(value)
    return sql


def GetMetric(field, data_dict, coreID):
    if ( (field =="AREA") or (field =="PERIMETER") or (field =="THICKNESS")):# or (field =="XCENTROID") or (field =="YCENTROID") ):
        value = 0
        for i in data_dict:
            if i["VALUE"] == coreID:
                value = i[field]
        if field == "AREA":
            return value / 4046.86
        elif ( (field == "PERIMETER") or (field == "THICKNESS") ):
            return value * 3.28084
        else:
            return value
    elif (field == "WATER_AREA" or field == "WET_AREA"):
        return (data_dict[coreID] * math.pow(float(arcpy.env.cellSize), 2)) / 4046.86
    else:
        if coreID in data_dict.keys():
            return data_dict[coreID]
        else:
            return 0


# -- Functions for Calculating Metrics
def Calc_Geometric_Stats(cores): ## Creates a list of dictionaries representing each core and its geometric attributes
    arcpy.AddMessage("Calculating some dumb geometric statistics for cores...")
    geostat_list    = []
    arcpy.AddField_management(cores, 'AREA','DOUBLE')
    arcpy.CalculateField_management(cores, 'AREA',"!shape.area!",'PYTHON_9.3')
    arcpy.AddField_management(cores, 'PERIMETER','DOUBLE')
    arcpy.CalculateField_management(cores, 'PERIMETER',"!shape.length!",'PYTHON_9.3')
    corelines = arcpy.FeatureToLine_management(cores, ScratchName('corelines'))
    coreEuc = EucDistance(corelines)
    arcpy.Delete_management(corelines)
    geostat_table = ZonalStatisticsAsTable(cores,"Value", coreEuc,  ScratchName("geostats"),"DATA", "MAXIMUM" )
    arcpy.Delete_management(coreEuc)
    arcpy.AddField_management(geostat_table, 'THICKNESS','DOUBLE')
    arcpy.CalculateField_management(geostat_table, 'THICKNESS',"!MAX!",'PYTHON_9.3')
    arcpy.JoinField_management(geostat_table,"VALUE",cores,'VALUE',"VALUE;AREA;PERIMETER;THICKNESS")
    fields          = ["VALUE",
                       "AREA",
                       "PERIMETER",
                       "THICKNESS"]
    with arcpy.da.SearchCursor(geostat_table, fields) as cursor:
        for row in cursor:
            obj_list = []
            for i in range(0, len(fields)):
                obj_list.append( (fields[i], row[i]) )
            geostat_list.append(dict(obj_list))

    arcpy.Delete_management(geostat_table)
    arcpy.AddMessage("Geometric statistics calculated.")
    return geostat_list  ## In the format: [ {'VALUE': 0, 'AREA': 00.0, ... }, {'VALUE': 1, 'AREA': 00.0, ... }, ... ]


def Calc_Compactness(data_dict, coreID, modifier):
    # Find radius of a circle with a circumference the same as the core's circumference
    for i in data_dict:
        if i["VALUE"] == coreID:
            p = float(i["PERIMETER"]) * 3.28084  # perimeter of the core, in feet
            a = float(i["AREA"]) / 4046.86 # area of the core, in acres
            radius     = float(p/(2*math.pi))
    # Find the areas of a circle with such a radius
    area_circle = math.pi * math.pow(radius, 2)
    compactness = a/area_circle
    return compactness * modifier


def Calc_PA_Ratio(data_dict, coreID):
    for i in data_dict:
        if i["VALUE"] == coreID:
            p = float(i["PERIMETER"]) * 3.28084  # perimeter of the core, in feet
            a = float(i["AREA"]) / 4046.86 # area of the core, in acres
            return p / a


def Calc_Topo_Diversity(cores, dem):
    arcpy.AddMessage("Calculating topographic diversity of cores...")
    topostat_dict = {}
    topostat_table = ZonalStatisticsAsTable(cores, "VALUE", dem, ScratchName("topostats"), "DATA", "STD")
    with arcpy.da.SearchCursor(topostat_table, ["VALUE", "STD"]) as cursor:
        for row in cursor:
            topostat_dict[row[0]] = row[1]

    arcpy.Delete_management(topostat_table)
    arcpy.AddMessage("Topographic diversity calculated.")
    return topostat_dict


def Count_Water_Cells(cores, nhd_area, nhd_waterbodies):
    arcpy.AddMessage("Finding the amount of aquatic habitat in each core...")
    # Convert the NHD Area data to raster
    nhd_area_path = arcpy.Describe(nhd_area).path
    WATER_AREA_FTYPES = [336, ## CanalDitch
                         445, ## SeaOcean
                         460] ## StreamRiver
    area_sql          = BuildSQL2(nhd_area_path, WATER_AREA_FTYPES, "FTYPE", "OR")
    arcpy.MakeFeatureLayer_management(nhd_area, "nhd_area", area_sql)
    nhd_area_r        = arcpy.PolygonToRaster_conversion("nhd_area", "FTYPE", ScratchName("nhd_area_r"))
    nhd_area_remap    = RemapRange([[0, 999999, 1], ["NODATA", 0]]) ## to convert NODATA cells to '0'
    nhd_area_reclass  = Reclassify(nhd_area_r, "VALUE", nhd_area_remap)

    # Convert the NHD Waterbody data to raster
    nhd_waterbodies_path = arcpy.Describe(nhd_waterbodies).path
    WATERBODY_FTYPES = [390, ## LakePond
                        436] ## Reservoir
    waterbodies_sql          = BuildSQL2(nhd_waterbodies_path, WATERBODY_FTYPES, "FTYPE", "OR")
    arcpy.MakeFeatureLayer_management(nhd_waterbodies, "nhd_waterbodies", waterbodies_sql)
    nhd_waterbodies_r        = arcpy.PolygonToRaster_conversion("nhd_waterbodies", "FTYPE", ScratchName("nhd_waterbodies_r"))
    nhd_waterbodies_remap    = RemapRange([[0, 999999, 1], ["NODATA", 0]]) ## to convert NODATA cells to '0'
    nhd_waterbodies_reclass  = Reclassify(nhd_waterbodies_r, "VALUE", nhd_waterbodies_remap)

    # Combine water features and count cells
    water_features = nhd_area_reclass | nhd_waterbodies_reclass
    aquastat_table = ZonalStatisticsAsTable(cores, "VALUE", water_features, ScratchName("aquastats"), "DATA", "SUM")
    aquastat_dict = {}
    with arcpy.da.SearchCursor(aquastat_table, ["VALUE", "SUM"]) as cursor:
        for row in cursor:
            aquastat_dict[row[0]] = row[1]

    arcpy.Delete_management(nhd_area_r)
    arcpy.Delete_management(nhd_waterbodies_r)
    arcpy.Delete_management(aquastat_table)
    arcpy.Delete_management("nhd_area")
    arcpy.Delete_management(nhd_area_reclass)
    arcpy.Delete_management("nhd_waterbodies")
    arcpy.Delete_management(nhd_waterbodies_reclass)
    arcpy.Delete_management(water_features)

    arcpy.AddMessage("Aquatic habitat found.")
    return aquastat_dict


def Count_Wetland_Cells(cores, wetlands):
    arcpy.AddMessage("Finding the amount of wetland habitat in each core...")
    # Convert the NWI wetlands data to raster
    wetlands_path = arcpy.Describe(wetlands).path
    WETLAND_TYPES = ["'Estuarine and Marine Wetland'",
                     "'Freshwater Emergent Wetland'" ,
                     "'Freshwater Forested/Shrub Wetland'"]
    wetland_sql       = BuildSQL2(wetlands_path, WETLAND_TYPES, "WETLAND_TYPE", "OR")
    arcpy.MakeFeatureLayer_management(wetlands, "wetlands", wetland_sql)
    wetlands_r        = arcpy.PolygonToRaster_conversion("wetlands", "ATTRIBUTE", ScratchName("wetlands_r"))
    wetlands_remap    = RemapRange([[0, 9999999, 1], ["NODATA", 0]]) ## to convert NODATA cells to '0'
    wetlands_reclass  = Reclassify(wetlands_r, "VALUE", wetlands_remap)

    wetlandstat_table = ZonalStatisticsAsTable(cores, "VALUE", wetlands_reclass, ScratchName("wetlandstats"), "DATA", "SUM")
    wetlandstat_dict = {}
    with arcpy.da.SearchCursor(wetlandstat_table, ["VALUE", "SUM"]) as cursor:
        for row in cursor:
            wetlandstat_dict[row[0]] = row[1]

    arcpy.Delete_management(wetlands_r)
    arcpy.Delete_management("wetlands")
    arcpy.Delete_management(wetlands_reclass)
    arcpy.Delete_management(wetlandstat_table)

    arcpy.AddMessage("Wetland habitat found.")
    return wetlandstat_dict

def Calc_Soil_Diversity(cores, soils_r):
    arcpy.AddMessage("Calculating soil diversity of cores...")
    soilstat_dict   = {}
    soilstat_table  = ZonalStatisticsAsTable(cores, "VALUE", soils_r, ScratchName("soilstats"), "DATA", "ALL")
    with arcpy.da.SearchCursor(soilstat_table, ["VALUE", "VARIETY"]) as cursor:
        for row in cursor:
            soilstat_dict[row[0]] = row[1]
    arcpy.AddMessage("Soil diversity calculated.")
    arcpy.Delete_management(soilstat_table)
    return soilstat_dict


def Calc_Biodiveristy_Stats(cores, bdras):
    arcpy.AddMessage("Calculating Biodiversity Priority of cores...")
    biodiv_dict ={}
    biodiv_table  = ZonalStatisticsAsTable(cores, "VALUE", bdras, ScratchName("biodiv"), "DATA", "MAXIMUM")
    arcpy.AlterField_management(ScratchName("biodiv"),"MAX", "BiodiversityPriorityIndex", "Biodiversity Priority Index" )
    with arcpy.da.SearchCursor(ScratchName("biodiv"), ["VALUE", "BiodiversityPriorityIndex"]) as cursor:
        for row in cursor:
            biodiv_dict[row[0]] = row[1]
    arcpy.AddMessage("Biodiversity Priority Max calculated.")

    arcpy.Delete_management(biodiv_table)

    return biodiv_dict


def Calc_EcoSys_Redundancy(cores, inras):
    arcpy.AddMessage("Calculating ecoregion redundancy of cores...")
    ecosyst_table = ZonalStatisticsAsTable(cores, "VALUE", inras, ScratchName('erstats'), "DATA", "ALL")
    arcpy.AlterField_management(ScratchName('erstats'), 'MIN', 'EcolSystem_Redundancy', 'Ecological System Redundancy')
    ecolredund_dict = {}
    with arcpy.da.SearchCursor(ScratchName('erstats'), ["VALUE", 'EcolSystem_Redundancy']) as cursor:
        for row in cursor:
            ecolredund_dict[row[0]] = row[1]

    arcpy.Delete_management(ecosyst_table)
    return ecolredund_dict


def Calc_Speciesrichness_Stats(cores, srras):
    arcpy.AddMessage("Calculating Endemic Species Max of cores...")
    specrich_table = ZonalStatisticsAsTable(cores, "VALUE", srras, ScratchName('srstats'), "DATA", "MAXIMUM")
    arcpy.AlterField_management(ScratchName('srstats'), "MAX", "EndemicSpeciesMax", "Endemic Species Max")
    arcpy.AddMessage("Endemic Species Max calculated.")
    specrich_dict = {}
    with arcpy.da.SearchCursor(ScratchName('srstats'), ["VALUE", "EndemicSpeciesMax"]) as cursor:
        for row in cursor:
            specrich_dict[row[0]] = row[1]

    arcpy.Delete_management(specrich_table)
    return specrich_dict


def SWDivIdx(incoreras, corepoly, divras, divfld):

    skip = 200
    min = int(arcpy.GetRasterProperties_management(incoreras,"MINIMUM").getOutput(0))
    max = int(arcpy.GetRasterProperties_management(incoreras,"MAXIMUM").getOutput(0))
    combtypes = numpy.dtype([('Value',numpy.int32),
                             ('{}_SWI'.format(divfld), numpy.float)])
    arraylist = []
    useThis = r'in_memory\raz'
    tempfc = r'in_memory\tempFC'
    loop = 0
    for i in range(min,max,skip):
        arcpy.env.extent = incoreras
        bot = i
        top = i + skip
        if top >= max :
            top = max + 1

        if bot == 0:
            valstr = "VALUE >= {} AND VALUE < {}".format(1, top)
        else:
            valstr = "VALUE >= {} AND VALUE < {}".format(bot, top)

        valstr = "VALUE >= {} AND VALUE < {}".format(bot, top)
        arcpy.MakeFeatureLayer_management(corepoly, 'subpoly', valstr )
        result = arcpy.GetCount_management('subpoly')
        count = int(result.getOutput(0))
        if count == 0:
            continue
        tempCoresPoly = arcpy.CopyFeatures_management('subpoly', tempfc)
        arcpy.env.extent = tempfc
        arcpy.CopyRaster_management(incoreras, useThis)
        arcpy.Delete_management(tempfc)
        coreras = arcpy.sa.ExtractByAttributes(useThis,valstr)
        arcpy.Delete_management(useThis)
        tempsoils = ScratchName('soilsrassub')
        arcpy.CopyRaster_management(divras, tempsoils)
        corename = arcpy.Describe(coreras).name #os.path.basename(coreras)
        divname = os.path.basename(tempsoils)
        combras = arcpy.sa.Combine([coreras, tempsoils])
        if arcpy.Exists(coreras):
            arcpy.Delete_management(coreras)
        if arcpy.Exists(tempsoils):
            arcpy.Delete_management(tempsoils)

        fldinfo =arcpy.FieldInfo()
        combflds =arcpy.ListFields(combras)
        for field in combflds:
            fldinfo.addField(field.name, field.name, "VISIBLE", "")
        combrastblview = arcpy.MakeTableView_management(combras,"temprastab", "", "", fldinfo)
        flds = [corename, divname, 'COUNT']
        combarray = arcpy.da.TableToNumPyArray(combrastblview, flds, skip_nulls=True)

        if combarray.shape[0]== 0:
            continue
        numcores = numpy.unique(combarray[corename]).shape[0]
        shi_arr = numpy.zeros(numcores,combtypes)
        coreindex = 0
        for coreval in numpy.unique(combarray[corename]):
            corearray = combarray[numpy.where(combarray[corename]==coreval)]
            totarea = sum(corearray['COUNT'])
            percArea = corearray['COUNT']/totarea
            shi_arr['{}_SWI'.format(divfld)][coreindex] = -1 * sum(percArea * numpy.log(percArea))
            shi_arr['Value'][coreindex] = coreval
            coreindex += 1
        arraylist.append(shi_arr)
        #arcpy.Delete_management(combras)

    partarray = numpy.concatenate(tuple(arraylist),1)
    scrtbl = ScratchName('swi_tbl')
    if arcpy.Exists( scrtbl):
        arcpy.Delete_management(scrtbl)
    arcpy.da.NumPyArrayToTable(partarray, scrtbl)
    arcpy.env.extent = incoreras
    swi_dict = {}
    swifld = "{}_SWI".format(divfld)
    with arcpy.da.SearchCursor(scrtbl, ["VALUE", swifld]) as cursor:
        for row in cursor:
            swi_dict[row[0]] = row[1]

    arcpy.Delete_management(scrtbl)
    return swi_dict


def Calc_Sum_Field(data_dict, coreID):
    if coreID in data_dict.keys():
        return data_dict[coreID]
    else:
        return 0


def Measure_Streams_Lengths(cores_poly, nhd_flowlines):
    arcpy.AddMessage("Adding up stream lengths for each core...")
    slength_dict = {}
    nhd_flow_path = arcpy.Describe(nhd_flowlines).path
    arcpy.MakeFeatureLayer_management(nhd_flowlines, "stream_lengths", """{0} = 460 OR {0} = 558""".format(arcpy.AddFieldDelimiters(nhd_flow_path, "FType")))
    core_stream_intersect = arcpy.Intersect_analysis([cores_poly, "stream_lengths"], ScratchName("Cstream_int"), "ALL")
    length_table = arcpy.Statistics_analysis(core_stream_intersect, ScratchName("Slengthtable"), [["Shape_Length", "SUM"]], "Value")
    with arcpy.da.SearchCursor(length_table, ["Value", "Sum_Shape_Length"]) as cursor:
        for row in cursor:
            slength_dict[row[0]] = row[1]

    arcpy.AddMessage("Stream lengths calculated.")

    arcpy.Delete_management(core_stream_intersect)
    arcpy.Delete_management(length_table)
    return slength_dict


def Populate_Stats(table, fields):
    with arcpy.da.UpdateCursor(table, fields) as cursor:
        for row in cursor:
            coreID  = row[0]
            row[1]  = GetMetric("AREA", geo_stats, coreID)
            row[2]  = GetMetric("PERIMETER", geo_stats, coreID)
            row[3]  = GetMetric("THICKNESS", geo_stats, coreID)
            row[4]  = GetMetric("TOPO_STD", topo_stats, coreID)
            row[17]  = GetMetric('BiodiversityPriorityIndex', biodiv_stats, coreID)
            row[6]  = GetMetric("WATER_AREA", water_stats, coreID)
            row[7]  = GetMetric("WET_AREA", wetland_stats, coreID)
            row[8]  = GetMetric("SOIL_DIV", soil_stats, coreID)
            row[14] = GetMetric("SOIL_SWI", soilSWI_stats, coreID)
            row[9] = Calc_Compactness(geo_stats, coreID, 100000) # Use a modifier of 100,000 to make the values larger
            row[10] = Calc_PA_Ratio(geo_stats, coreID)
            row[11] = Calc_Sum_Field(stream_stats, coreID) * 3.28084 # convert to feet
            row[16] = GetMetric("EcolSystem_Redundancy", ecocsysred_stats, coreID)
            row[15] = GetMetric("EndemicSpeciesMax", specrich_stats, coreID)
            cursor.updateRow(row)


# Helper Functions -- Module #5
def Score(number, break_points, invert = False): # Assign a score 1-5 based in which quintile the metric falls
    if number <= break_points[0]:
        return 1 if not invert else 5
    elif number <= break_points[1]:
        return 2 if not invert else 4
    elif number <= break_points[2]:
        return 3 if not invert else 3
    elif number <= break_points[3]:
        return 4 if not invert else 2
    else:
        return 5 if not invert else 1


def Calc_Score_Field(dataset, fields, weights):
    arcpy.AddField_management(dataset, "SCORE", "DOUBLE")
    with arcpy.da.UpdateCursor(dataset, ["SCORE"] + fields) as cursor:
        for row in cursor:
            # Look up the metric and assign a score (1-5) based on quintiles
            area_score          = Score(row[1], quintile_breaks["AREA"])
            thick_score         = Score(row[2], quintile_breaks["THICKNESS"])
            topo_score          = Score(row[3], quintile_breaks["TOPO_STD"])
            species_score       = Score(row[4], quintile_breaks["BiodiversityPriorityIndex"])
            wetland_score       = Score(row[5]/row[1], quintile_breaks["WET_AREA"])
            soil_score          = Score(row[6], quintile_breaks["SOIL_SWI"])
            compact_score       = Score(row[7], quintile_breaks["COMPACT"])
            stream_score        = Score(row[8]/row[1], quintile_breaks["STRM_LEN"]) # Normalized by area (i.e. feet of stream per acre)
            rte_abund_score     = Score(row[9], quintile_breaks["EcolSystem_Redundancy"],True) # Normalize by area??
            rte_div_score       = Score(row[10], quintile_breaks["EndemicSpeciesMax"]) # Normalize by area??

            # Apply weights to scores and sum them
            core_scores         = []
            core_scores.append(area_score      * weights[0])
            core_scores.append(thick_score     * weights[1])
            core_scores.append(topo_score      * weights[2])
            core_scores.append(species_score   * weights[3])
            core_scores.append(wetland_score   * weights[4])
            core_scores.append(soil_score      * weights[5])
            core_scores.append(compact_score   * weights[6])
            core_scores.append(stream_score    * weights[7])
            core_scores.append(rte_abund_score * weights[8])
            core_scores.append(rte_div_score   * weights[9])
            row[0] = sum(core_scores)
            cursor.updateRow(row)


def Core_Class(cores):
    arcpy.AddField_management(cores, "CLASS", "TEXT", "", "", 20)
    with arcpy.da.UpdateCursor(cores, ["AREA", "WATER_AREA", "CLASS"]) as cursor:
        for row in cursor:
            if row[0]-row[1] < 100:
                row[2] = "Habitat Fragment"
            elif row[0]-row[1] < 1000:
                row[2] = "Small"
            elif row[0]-row[1] < 10000:
                row[2] = "Medium"
            else:
                row[2] = "Large"
            cursor.updateRow(row)
#-------------------------------------------------------------------------------------------------------------------------------------------#
# Custom exceptions
class LicenseError(Exception):
    pass
class ProjectionError(Exception):
    pass
class ShapeTypeError(Exception):
    pass
class WetlandsOIDError(Exception):
    pass
class LinearUnitError(Exception):
    pass
class CoreOIDError(Exception):
    pass
class CustomSQLError(Exception):
    pass
class CustomWeightsError(Exception):
    pass
#-------------------------------------------------------------------------------------------------------------------------------------------#

'''==================================='''
'''               MAIN                '''
'''==================================='''

## PART 1: SETUP ##

# Check for scratch workspace and create one if the default does not exist
scriptpath  = sys.path[0]
toolpath    = os.path.dirname(scriptpath)

scratchws = arcpy.env.scratchGDB
arcpy.AddMessage("Scratch Workspace:\n" + scratchws)

arcpy.env.workspace     = out_workspace

# Check product info
product_info = arcpy.ProductInfo()
arcpy.AddMessage("Product Info: " + product_info)

# Check Spatial Analyst availability
sa_license_availability = arcpy.CheckExtension("Spatial")
if sa_license_availability == "Available":
    arcpy.CheckOutExtension("Spatial")
    arcpy.AddMessage("The Spatial Analyst extension is available and has been successfully checked out.")
else:
    arcpy.AddError("The script detected the availability of the Spatial Analyst extension as: " + sa_license_availability)
    raise LicenseError # Stops script

# Check if any input layers are not in the NAD 1983 UTM Zone 17N projection (warning), or are not projected at all(error)
input_data_list = [area_of_interest,
                   lulc_layer,
                   noaa_ccap,
                   nwi_wetlands,
                   ssurgo_soils,
                   nhd_flowlines,
                   nhd_waterbodies,
                   nhd_area,
                   dem,
                   local_roads,
                   railroads,
                   building_locations,
                   specrich_ras,
                   ecosysred_ras,
                   biodiv_ras ]

red_flags = {"nonmeter":[]}
for layer in input_data_list:
    if layer:
        desc = arcpy.Describe(layer)
        if desc.spatialReference.linearUnitCode != 0: # If linear units are not undefined
            if desc.spatialReference.linearUnitName != "Meter":
                red_flags["nonmeter"].append(desc.name)


if len(red_flags["nonmeter"]) == 0:
    arcpy.AddMessage("All input layers have linear units of Meters.")
else:
    if len(red_flags["nonmeter"]) > 0:
        nonmeter_layers = ""
        for layer in red_flags["nonmeter"]:
            nonmeter_layers = nonmeter_layers + layer + "\n"
        arcpy.AddError("All input layers should use a map projection that uses meters as its linear unit (to prevent errors in area calculations).\nPlease project the following layers using meters as the linear units.:\n" + nonmeter_layers)
        raise LinearUnitError # Stop script

EndModule("Part 1: Setup")

## PART 2: IDENTIFY CORE HABITAT ##

# Buffer AOI, if required
CheckBufferDistance(AOI_BUFFER_DISTANCE)
if AOI_BUFFER_DISTANCE > 0:  ## If a buffer is specified...
    aoi_buffer = arcpy.Buffer_analysis(area_of_interest, "aoi_buffer", "{0} Miles".format(AOI_BUFFER_DISTANCE), "", "", "ALL")
    arcpy.AddMessage("Area of Interest buffered")
    arcpy.env.extent = arcpy.Describe('aoi_buffer').extent
    lulc_aoi   = ExtractByMask(lulc_layer, aoi_buffer)
    arcpy.AddMessage("Land cover data for Area of Interest extracted.")
else: ## Use the AOI boundary to clip
    arcpy.env.extent = arcpy.Describe(area_of_interest).extent
    lulc_aoi   = ExtractByMask(lulc_layer, area_of_interest)
    arcpy.AddMessage("Land cover data for Area of Interest extracted.")

arcpy.env.cellSize   = lulc_aoi
arcpy.env.extent     = lulc_aoi
arcpy.env.snapRaster = lulc_aoi

# Convert wetlands data to raster
arcpy.AddMessage("Finding wetlands in area of interest...")
wetlands_nat_sql = """{0} NOT LIKE '%d' AND {0} NOT LIKE '%f' AND {0} NOT LIKE '%r' AND {0} NOT LIKE '%x' AND {1} <> 'Estuarine and Marine Deepwater'""".format(arcpy.AddFieldDelimiters(nwi_wetlands, "ATTRIBUTE"), arcpy.AddFieldDelimiters(nwi_wetlands, "WETLAND_TYPE")) 
## Use only natural wetlands, based on NWI modifiers
arcpy.MakeFeatureLayer_management(nwi_wetlands, "wetlands_nat", wetlands_nat_sql)
wetlands_r = arcpy.FeatureToRaster_conversion("wetlands_nat", WetlandsRasterField(nwi_wetlands), ScratchName("wetlands_r"))
arcpy.AddMessage("Wetlands in area of interest found.")
wetlands_remap  = RemapRange([[0, 9999999, 1], ["NODATA", 0]])

# If the county is coastal, use C-CAP data to find natural landscapes along the coast (but fragmenting features from the LULC layer can override these)
if aoi_is_coastal:
    arcpy.AddMessage("Extracting coastal land uses from NLCD data...")
    seaocean_sql     = """{0} = 445""".format(arcpy.AddFieldDelimiters(nhd_area, "FType")) ## Find 'SeaOcean' polygons
    arcpy.MakeFeatureLayer_management(nhd_area, "SeaOcean", seaocean_sql)
    ocean_buf        = EucDistance("SeaOcean", 5001) ## Only include C-CAP data within 5 kilometers of the shoreline
    ccap_shore_zone  = ExtractByMask(noaa_ccap, ocean_buf)
    arcpy.Delete_management(ocean_buf)
    # Disregard classes: Background, Unclassified, Developed, Pasture/Hay, and Cultivated 
    coastal_lu  = Reclassify(ccap_shore_zone, "VALUE", RemapRange([
        [0, 12, 1],
        [13, 29, 0],
        [31, 95, 1],
        ["NODATA", 0] ]))

    arcpy.AddMessage("Coastal data has been extracted, you chump!")
# Remaps for reclassifying both CDL and NCLD data
if (LULC_SOURCE == "Cropland Data Layer"):
    lulc_remap       = RemapRange(
        [[0, 110, 0],
        [110, 112, 1],
        [112, 130, 0],
        [130, 131, 1],
        [131, 140, 0],
        [140, 143, 1],
        [143, 151, 0],
        [151, 152, 1],
        [152, 189, 0],
        [189, 190, 1],
        [190, 194, 0],
        [194, 195, 1],
        [195, 256, 0],
        ["NODATA", 0]])
elif (LULC_SOURCE == "NLCD"):
    lulc_remap = RemapRange(
        [[0, 10, 0], # Disregard unclassified values
         [11, 19, 1], # Includes Open water, ice/snow 
         [19, 29, 0], # Disregards Developed
         [29, 95, 1], # Includes Barren, Forested, Scrub, Herbaceous, Cultivated, Wetlands
       # [79, 89, 0],
       # [89, 99, 1],
         [96, 256, 0],
         ["NODATA", 0]])
elif (LULC_SOURCE == "C-CAP"):
    lulc_remap = RemapRange(
       [[0, 0, 1], # Include classes: Background
        [1, 5, 0], # Disregard classes: Unclassified, Developed
        [6, 23, 1], # Include Pasture/Hay, Cultivated, Grassland, Forested, Scrub
        #[20, 20, 1], # Disregard Barren
        #[21, 21, 1], # Include open water
        # [22, 23, 1], # Includes Palustrine and Estuarine Beds
        [24, 256, 0], # All other unclassified values disregarded
        ["NODATA", 0] ])
else: ## User provided LULC
    lulc_remap       = RemapValue(
        [[1, 1],
        [2, 0],
        [3, 0],
        ["NODATA", 0]])

# Combine all habitats that can contribute to a core in a single raster
arcpy.AddMessage("Combining the core habitats, you fool of a Took!...")
forest_and_water          = Reclassify(lulc_aoi, "VALUE", lulc_remap)
wetlands                  = Reclassify(wetlands_r, "VALUE", wetlands_remap)
if aoi_is_coastal:
    core_habitat          = forest_and_water | wetlands | coastal_lu
else:
    core_habitat          = forest_and_water | wetlands


arcpy.AddMessage("Dumb core habitats combined")

# Save the output core habitat areas
core_habitat.save("core_habitat")
if aoi_is_coastal:
    arcpy.Delete_management(coastal_lu)


EndModule("Part 2: Identifying Core Habitat, you stuck up, half-witted, scruffy-looking nerf herder!")

## PART 3: FRAGMENTATION ##
if AOI_BUFFER_DISTANCE > 0:
    aoipoly = aoi_buffer
else:
    aoipoly = area_of_interest
FRAGMENT_RASTERS = []  ## To hold the various rasters that represent fragmenting features (roads, urban areas, etc.). Used for merging them together

# Clip fragmenting features (roads, buildings, rail lines)
arcpy.AddMessage("Clipping local roads data to Area of Interest...")
lroads_aoi     = arcpy.Clip_analysis(local_roads, aoipoly, ScratchName("lroads_clip"))
arcpy.AddMessage("Local roads data extracted for Area of Interest")
if local_roads:
    rdcount = int(arcpy.GetCount_management(lroads_aoi)[0])
if rdcount > 0:
    # Find paved local roads and convert to raster
    arcpy.AddMessage("Finding paved local roads...")
    arcpy.MakeFeatureLayer_management(lroads_aoi,"paved_lroads")#, lpaved_sql)
    arcpy.AddField_management ("paved_lroads", "RASNUM", "LONG")
    arcpy.CalculateField_management ("paved_lroads", "RASNUM", "1", 'PYTHON_9.3')
    paved_lroads_r           = arcpy.PolylineToRaster_conversion("paved_lroads", "RASNUM", ScratchName("lroadspvd"))
    lroads_remap             = RemapValue([[0, 0], [1, 1], ["NODATA", 0]]) ## to convert NODATA cells to '0'
    paved_lroads_reclass     = Reclassify(paved_lroads_r, "VALUE", lroads_remap)
    arcpy.Delete_management(paved_lroads_r)
    arcpy.Delete_management(lroads_aoi)
    FRAGMENT_RASTERS.append(paved_lroads_reclass)
    arcpy.AddMessage("Fragmenting local roads found.")



# Covert railroads to raster
if railroads:
    rrtestlyr = arcpy.MakeFeatureLayer_management(railroads, 'rrlyr')
    arcpy.SelectLayerByLocation_management('rrlyr','intersect', aoipoly)
    rcount = int(arcpy.GetCount_management('rrlyr')[0])
if rcount > 0:
    arcpy.AddMessage("Finding rail lines...")
    desc_rr     = arcpy.Describe(railroads)
    railroads_r = arcpy.PolylineToRaster_conversion(railroads, desc_rr.OIDFieldName, ScratchName("RR_r"))
    rr_remap    = RemapRange([[0, 999999, 1], ["NODATA", 0]])
    rr_reclass  = Reclassify(railroads_r, "VALUE", rr_remap)
    arcpy.Delete_management(railroads_r)
    FRAGMENT_RASTERS.append(rr_reclass)
    arcpy.AddMessage("Rail lines found.")

# Convert building locations to raster
building_count = 0
if building_locations:
    bldtestlyr = arcpy.MakeFeatureLayer_management(building_locations, 'blyr')
    arcpy.SelectLayerByLocation_management('blyr','intersect', aoipoly)
    building_count = int(arcpy.GetCount_management('blyr')[0])
if building_locations and building_count > 0:
    arcpy.AddMessage("Finding building locations")
    bldg_layer = arcpy.Copy_management(building_locations, ScratchName("bldg_layer"))
    if arcpy.Describe(building_locations).shapeType == "Polygon":
        if product_info == "ArcInfo":
            bldg_points = arcpy.FeatureToPoint_management(bldg_layer, ScratchName("bldg_points"), "CENTROID")
            arcpy.AddField_management(bldg_points, "RASNUM", "SHORT")
            arcpy.CalculateField_management(bldg_points, "RASNUM", "1", 'PYTHON_9.3')
            bldg_r      = arcpy.FeatureToRaster_conversion(bldg_points, "RASNUM", ScratchName("bldg_r"))
        else:
            bldg_points = CreateCentroids(bldg_layer, ScratchName("bldg_points"))
            arcpy.AddField_management(bldg_points, "RASNUM", "SHORT")
            arcpy.CalculateField_management(bldg_points, "RASNUM", "1", 'PYTHON_9.3')
            bldg_r      = arcpy.FeatureToRaster_conversion(bldg_points, "RASNUM", ScratchName("bldg_r"))
    else:
        arcpy.AddField_management(bldg_layer, "RASNUM", "SHORT")
        arcpy.CalculateField_management(bldg_layer, "RASNUM", "1", 'PYTHON_9.3')
        bldg_r      = arcpy.FeatureToRaster_conversion(bldg_layer, "RASNUM", ScratchName("bldg_r"))
    bldg_remap      = RemapValue([[0, 0], [1, 1], ["NODATA", 0]])
    bldg_reclass    = Reclassify(bldg_r, "VALUE", bldg_remap)
    arcpy.Delete_management(bldg_r)
    FRAGMENT_RASTERS.append(bldg_reclass)
    arcpy.AddMessage("Building locations found.")

# Extract Fragmenting Features from LULC
arcpy.AddMessage("Finding fragmenting features from land cover data...")
if (LULC_SOURCE == "Cropland Data Layer"):
    lulc_frag_remap     = RemapRange(
    [[0, 110, 1],
    [110, 112, 0],
    [112, 130, 1],
    [130, 152, 0],
    [152, 176, 1],
    [190, 195, 0],
    [195, 256, 1],
    ["NODATA", 0]])
elif (LULC_SOURCE == "NLCD"):
    lulc_frag_remap = RemapRange(
    [[0, 19, 0], # Open Water & Ice/Snow not fragmenting
    [19, 29, 1], # Development is fragmenting
    [29, 95, 0], # Barren, Forested, Scrub, Grassland, Cultivated, Wetlands not fragmenting
  # [79, 89, 1],
    [96, 255, 0], # All other unclassified values not fragmenting
    ["NODATA", 0]])
elif (LULC_SOURCE == "C-CAP"):
    lulc_frag_remap = RemapRange(
    [[0, 0, 0], # Disregards Background
    [1, 5, 1], # Includes Unclassified & Developed as fragmenting features 
    [6, 19, 0],
    [20, 20, 0],
    [21, 255, 0],
    #[6, 255, 0], # Cultivated, Grassland, Forested, Scrub, Barren, Wetlands, Barren, Water not fragmenting  
    ["NODATA", 0]])
else: ## User provided LULC
    lulc_frag_remap     = RemapValue(
        [[1, 0],
        [2, 1],
        [3, 1],
        ["NODATA", 0]])


frag_lu         = Reclassify(lulc_aoi, "VALUE", lulc_frag_remap)
frag_lu.save('frag{}'.format(OUTPUT_NAME))

# Extract 'Developed' areas on their own so they can override coastal land use
if (LULC_SOURCE == "Cropland Data Layer"):
    lulc_developed_remap    = RemapRange(
    [[0, 120, 0],
    [120, 124, 1],
    [124, 256, 0],
    ["NODATA", 0]])
elif (LULC_SOURCE == "NLCD"):
    lulc_developed_remap    = RemapRange(
    [[0, 20, 0],
    [20, 25, 1],
    [26, 256, 0],
    ["NODATA", 0]])
#below added for ccap#
elif (LULC_SOURCE == "C-CAP"):
    lulc_developed_remap    = RemapRange(
    [[0, 0, 0],
    [1, 5, 1],
    [6, 256, 0],
    ["NODATA", 0]])
#above added for ccap#        
else: ## User provided LULC
    lulc_developed_remap    = RemapValue(
        [[1, 0],
        [2, 1],
        [3, 0],
        ["NODATA", 0]])
lulc_developed          = Reclassify(lulc_aoi, "VALUE", lulc_developed_remap)
#arcpy.Delete_management(lulc_aoi)

if aoi_is_coastal:
    ## If an area is natural in the coastal LU layer AND not 'developed' in the main LU layer, then consider it non-fragmenting
    frag_lu_mod = Con (coastal_lu,Con(lulc_developed,0,frag_lu,'Value <> 1'),frag_lu,'Value=1')
    FRAGMENT_RASTERS.append(frag_lu_mod)
else:
    FRAGMENT_RASTERS.append(frag_lu)
arcpy.AddMessage("Fragmenting features found.")

# Combine fragmenting features rasters
arcpy.AddMessage("Combining all fragmenting features...")
frag_raster = CellStatistics(FRAGMENT_RASTERS, "MAXIMUM", "NODATA")
arcpy.AddMessage("Fragmenting features combined.")


# Buffer (euclidean distance) the fragmentation layer
arcpy.AddMessage("Buffering fragmenting features to identify edge habitat...")
frag_raster_extract    = ExtractByAttributes(frag_raster, "VALUE = 1")
frag_dist              = EucDistance(frag_raster_extract)
edge_habitat           = frag_dist <= BUFFER_DISTANCE
arcpy.AddMessage("Edge habitat identified.")



# subtract out edge habitat from core habitat to find interior habitat
arcpy.AddMessage("Removing edge habitat from core habitat...")
cores_raw              = Con((core_habitat - edge_habitat) < 0, 0, core_habitat - edge_habitat)
arcpy.AddMessage("Edge habitat removed.")

# remove any cells outside the area of interest (created from raster conversions)
if AOI_BUFFER_DISTANCE > 0:  ## If a buffer is specified...
    cores_aoi   = ExtractByMask(cores_raw, aoi_buffer)
else: ## Use the AOI boundary to clip
    cores_aoi   = ExtractByMask(cores_raw, area_of_interest)

# region group and set up cores
arcpy.AddMessage("Grouping cores...")
cores                  = RegionGroup(cores_aoi, "FOUR", "WITHIN", "NO_LINK", 0)
cores_culled           = ExtractByAttributes(cores, "COUNT >= {0}".format(FindMinCoreCells(MIN_CORE_SIZE, cores))) ## Remove regions of less than the minimum core size
core_groups_reclass    = SetNull(cores_culled, cores_culled, "VALUE = 0")

# Add back in edge habitat
cores_ca = EucAllocation(core_groups_reclass, CORE_EDGE)
cores_final = Con((cores_ca & core_habitat),1) # cores_ca)
arcpy.AddMessage('Merging touching and overlapping cores')
regcores = arcpy.sa.RegionGroup(cores_final, 'FOUR', "WITHIN", "NO_LINK", 0)
arcpy.BuildRasterAttributeTable_management(regcores)
mincorecells = math.ceil(FindMinCoreCells(MIN_CORE_SIZE, regcores))
mincores = arcpy.sa.ExtractByAttributes(regcores, "COUNT >= {}".format(mincorecells))
arcpy.Delete_management(regcores)

arcpy.AddMessage('Saving HabitatFragments raster. Minimum fragment size: {} acres'.format(MIN_FRAGMENT_SIZE))
habminfrag = Con(IsNull(frag_raster_extract), core_habitat)
habfragmask = ExtractByMask(habminfrag,area_of_interest)
fragras = arcpy.sa.Con(IsNull(mincores),habfragmask)
fragrasrg = arcpy.sa.RegionGroup(fragras, "FOUR", "WITHIN", "NO_LINK", 0)
fragsculled = ExtractByAttributes(fragrasrg, "COUNT >= {0}".format(FindMinCoreCells(MIN_FRAGMENT_SIZE, fragrasrg)))
arcpy.Delete_management(fragras)
arcpy.Delete_management(fragrasrg)
arcpy.Delete_management(habfragmask)
arcpy.Delete_management(habminfrag)
fragsculled.save(HAB_FRAG_NAME)


# Assign a new (sequential) unique value to each core
arcpy.AddMessage("Creating a unique ID for each core...")
arcpy.BuildRasterAttributeTable_management(mincores, "Overwrite")
desc = arcpy.Describe(mincores)
if desc.hasOID:
    oid_field = desc.OIDFieldName
else:
    arcpy.AddError("The cores raster must have a OID field to calculate unique statistics.")
    raise CoreOIDError
arcpy.AddField_management(mincores, "CoreID", "LONG")
arcpy.CalculateField_management(mincores, "CoreID", '!{0}!'.format(oid_field), "PYTHON_9.3")
cores_unique = Lookup(mincores, "CoreID")

arcpy.AddMessage("Unique IDs created for cores.")

# Save cores
full_core_name = os.path.join(arcpy.env.workspace, OUTPUT_NAME)
full_core_poly_name = os.path.join(arcpy.env.workspace, OUTPUT_NAME_POLY)
cores_unique.save(full_core_name)


# Clean up intermediate outputs before Module 4
arcpy.AddMessage("Deleting intermediate data...")
arcpy.Delete_management(frag_raster_extract)
arcpy.Delete_management(frag_dist)
arcpy.Delete_management("wetlands_nat")
if aoi_is_coastal:
    arcpy.Delete_management("SeaOcean")
    arcpy.Delete_management(ccap_shore_zone)
    arcpy.Delete_management(coastal_lu)
    arcpy.Delete_management(frag_lu_mod)
arcpy.Delete_management(forest_and_water)
arcpy.Delete_management(wetlands)
arcpy.Delete_management(core_habitat)
if rdcount > 0:
    arcpy.Delete_management(paved_lroads_reclass)
if rcount > 0:
    arcpy.Delete_management(rr_reclass)
if building_count > 0:
    arcpy.Delete_management(bldg_reclass)
arcpy.Delete_management(frag_lu)
arcpy.Delete_management(frag_raster)
arcpy.Delete_management(frag_raster_extract)
arcpy.Delete_management(frag_dist)
arcpy.Delete_management(edge_habitat)
arcpy.Delete_management(cores_raw)
arcpy.Delete_management(cores_culled)
arcpy.Delete_management(cores_aoi)
arcpy.Delete_management(core_groups_reclass)
arcpy.Delete_management(cores_ca)
arcpy.Delete_management(cores_final)
arcpy.Delete_management(cores)
arcpy.Delete_management(mincores)
arcpy.Delete_management(lulc_developed)
if building_locations:
    arcpy.Delete_management(bldg_reclass)
if AOI_BUFFER_DISTANCE > 0 and arcpy.Exists(aoi_buffer):
    arcpy.Delete_management(aoi_buffer)

arcpy.AddMessage("Intermediate data deleted...")

EndModule("Part 3: Assessing Fragmentation")

## PART 4: CALCULATE CORE METRICS AND RANKINGS ##

# Create a polygon version of cores (with optional simplification i.e. removing holes)
#TEMP_SMOOTH_NAME = "cores_smooth"
#if smoothing_passes > 0:
#    smoothing_counter = 0
#    while smoothing_counter < smoothing_passes:
#        if smoothing_counter == 0:
#            smoothed_cores = MajorityFilter(cores_unique, "FOUR", "HALF")
#            smoothed_cores.save(ScratchName(TEMP_SMOOTH_NAME + str(smoothing_counter)))
#        else:
#            smoothed_cores = MajorityFilter(ScratchName(TEMP_SMOOTH_NAME + str(smoothing_counter-1)), "FOUR", "HALF")
#            smoothed_cores.save(ScratchName(TEMP_SMOOTH_NAME + str(smoothing_counter)))
#        smoothing_counter += 1
#    cores_poly_tmp = arcpy.RasterToPolygon_conversion(ScratchName(TEMP_SMOOTH_NAME + str(smoothing_counter-1)), ScratchName("cores_poly_tmp"), simplify_polygons)
#    if arcpy.Exists(full_core_poly_name):
#        arcpy.Delete_management(full_core_poly_name)
#    arcpy.Dissolve_management(cores_poly_tmp, full_core_poly_name, "gridcode")
#else:
#    cores_poly_tmp = arcpy.RasterToPolygon_conversion(cores_unique, ScratchName("cores_poly_tmp"), simplify_polygons)
#    if arcpy.Exists(full_core_poly_name):
#        arcpy.Delete_management(full_core_poly_name)
#    results = arcpy.Dissolve_management(cores_poly_tmp, full_core_poly_name, "gridcode")
#    flds = arcpy.ListFields(full_core_poly_name)
#    for fld in flds:
#        print(fld.name)
#
#arcpy.Delete_management(cores_poly_tmp)
#
#result = arcpy.AlterField_management(full_core_poly_name,'gridcode', 'Value', 'Value')

# Calculate statistics
#geo_stats       = Calc_Geometric_Stats(full_core_poly_name)
#topo_stats      = Calc_Topo_Diversity(cores_unique, dem)
#biodiv_stats  = Calc_Biodiveristy_Stats(full_core_poly_name,biodiv_ras)
#water_stats     = Count_Water_Cells(cores_unique, nhd_area, nhd_waterbodies)
#wetland_stats   = Count_Wetland_Cells(cores_unique, nwi_wetlands)
#soil_stats      = Calc_Soil_Diversity(cores_unique, ssurgo_soils)
#stream_stats    = Measure_Streams_Lengths(full_core_poly_name, nhd_flowlines)
#ecocsysred_stats = Calc_EcoSys_Redundancy(full_core_poly_name, ecosysred_ras)
#specrich_stats = Calc_Speciesrichness_Stats(full_core_poly_name,specrich_ras)
#soilSWI_stats = SWDivIdx(full_core_name,full_core_poly_name,ssurgo_soils,'SOIL')

#new_fields = ["AREA",       #-----1
#              "PERIMETER",  #-----2
#              "THICKNESS",  #-----3
#              "TOPO_STD",   #-----4
#              "S_RICHNESS", #-----5
#              "WATER_AREA", #-----6
#              "WET_AREA",   #-----7
#              "SOIL_DIV",   #-----8
#              "COMPACT",    #-----9
#              "PA_RATIO",   #-----10
#              "STRM_LEN",   #-----11
#              "RTE_ABUND",  #-----12
#              "RTE_DIV",    #-----13
#              "SOIL_SWI",   #-----14
#              "EndemicSpeciesMax",    #-----15
#              'EcolSystem_Redundancy',  #-----16
#              'BiodiversityPriorityIndex']    #-----17
#
#
#for field in new_fields:
#    arcpy.AddMessage("Adding new field: {0}".format(field))
#    arcpy.AddField_management (full_core_name, field, "DOUBLE")
#    arcpy.AddField_management (full_core_poly_name, field, "DOUBLE")
#
#arcpy.AddMessage("Populating fields...")
#Populate_Stats(full_core_name, ["VALUE"] + new_fields)
#Populate_Stats(full_core_poly_name, ["VALUE"] + new_fields)
#arcpy.AddMessage("Fields populated with core statistics.")
#
# Rank Cores
#arcpy.AddMessage("Calculating Core Quality Index...")
#
#ranking_fields  = ["AREA",
#                   "THICKNESS",
#                   "TOPO_STD",
#                   "BiodiversityPriorityIndex",
#                   "WET_AREA",
#                   "SOIL_SWI",
#                   "COMPACT",
#                   "STRM_LEN",
#                   "EcolSystem_Redundancy",
#                   "EndemicSpeciesMax"]
# Create and populate a dictionary with all values for each attribute field
#stats_master    = {}
#quintile_breaks = {}
#for i in ranking_fields:
#    stats_master[i]     = []
#    quintile_breaks[i]  = []
#with arcpy.da.SearchCursor(full_core_poly_name, ranking_fields) as cursor:
#    for row in cursor:
#        stats_master["AREA"].append(row[0])
#        stats_master["THICKNESS"].append(row[1])
#        stats_master["TOPO_STD"].append(row[2])
#        stats_master["BiodiversityPriorityIndex"].append(row[3])
#        stats_master["WET_AREA"].append(row[4]/row[0]) # Normalized by area (% Wetlands)
#        stats_master["SOIL_SWI"].append(row[5])
#        stats_master["COMPACT"].append(row[6])
#        stats_master["STRM_LEN"].append(row[7]/row[0]) # Normalized by area (i.e. feet of stream per acre)
#        stats_master["EcolSystem_Redundancy"].append(row[8])  # Normalize by area??
#        stats_master["EndemicSpeciesMax"].append(row[9])  # Normalize by area??
## Find the break points for each attribute field to create quintiles
#for i in ranking_fields:
#    arcpy.AddMessage(i)
#    values = numpy.array(stats_master[i])
#    for break_point in range(20, 101, 20):
#        quintile_breaks[i].append(numpy.percentile(values, break_point))
#    quintile_breaks[i].sort()
#
# If user has provided custom weights, use those
#
#weights = DEFAULT_WEIGHTS
#
## Create a new 'Score' field and calculate the score using weights
#Calc_Score_Field(full_core_poly_name, ranking_fields, weights)
#Calc_Score_Field(full_core_name, ranking_fields, weights)
#arcpy.AddMessage("Core Quality Index calculated.")
#
## Create a core class field and populate with either: Habitat Fragment, Small, Medium, or Large
#Core_Class(full_core_poly_name)
#Core_Class(full_core_name)
#arcpy.AddMessage("Cores classified into Small, Medium, Large, and Habitat Fragments.")
#
#arcpy.AddMessage ('Output saved in output workspace: {}'.format(out_workspace))
#arcpy.AddMessage ('Intact habitat core polygons saved to: {}'.format(OUTPUT_NAME_POLY))
#arcpy.AddMessage ('Intact habitat core raster saved to: {}'.format(OUTPUT_NAME))
#arcpy.AddMessage ('Habitat fragments not within intact habitat cores saved to: {}'.format(HAB_FRAG_NAME))
#
#
#EndModule("Part 4: Calculating Core Statistics")
#
#arcpy.AddMessage("Attempting to delete intermediate data...")
#if not arcpy.env.scratchWorkspace:
#    arcpy.Delete_management(scratchws)
#
## Create a log file for the model run
#logtime = str(int(time.time()))
#arcpy.AddMessage("Logging details of model run here: "+os.path.join(toolpath, 'Log_'+logtime+'.txt'))
#f = open(os.path.join(toolpath, 'Log_'+logtime+'.txt'), 'w')
#f.write("Log for Model Run at "+str(datetime.now())+'\n'+'='*30+'\n')
#f.write("\nMODEL PARAMETERS\n")
#f.write("   Workspace: " + out_workspace +'\n')
#f.write('   Output saved in output workspace: {}'.format(out_workspace))
#f.write('   Intact habitat core polygons saved to: {}'.format(OUTPUT_NAME_POLY))
#f.write('   Intact habitat core raster saved to: {}'.format(OUTPUT_NAME))
#f.write('   Habitat fragments not with intact habitat cores saved to: {}'.format(HAB_FRAG_NAME))
#f.write("   Area of Interest: " + area_of_interest +'\n')
#f.write("   AOI Buffer Distance: " +str(AOI_BUFFER_DISTANCE) +'\n')
#f.write("   Land Use / Land Cover Data: " +lulc_layer +'\n')
#f.write("   Land Use / Land Cover Data Source: " +LULC_SOURCE +'\n')
#f.write("   AOI is Coastal?: " +str(aoi_is_coastal) +'\n')
#f.write("   NWI Wetlands: " +nwi_wetlands +'\n')
#f.write("   SSURGO Soils: " +ssurgo_soils +'\n')
#f.write("   NHD Flowlines: " +nhd_flowlines +'\n')
#f.write("   NHD Waterbodies: " +nhd_waterbodies +'\n')
#f.write("   NHD Areas: " +nhd_area +'\n')
#f.write("   Species Richness: " + specrich_ras +'\n')
#f.write("   Ecological System Redundancy: " + ecosysred_ras +'\n')
#f.write("   Digital Elevation Model: " +dem +'\n')
#f.write("   Local Roads File: " +local_roads +'\n')
#f.write("   Railroads: " +railroads +'\n')
#f.write("   Building Locations: " +building_locations +'\n')
#f.write("\nRANKING WEIGHTS\n")
#for i in ranking_fields:
#    f.write("   "+i+': '+str(weights[ranking_fields.index(i)])+'\n')
#f.write("\nCORE ATTRIBUTE PERCENTILES\n")
#for i in ranking_fields:
#    f.write("   "+i+': \n')
#    x = 20
#    for value in quintile_breaks[i]:
#        f.write('       {0}th percentile: '.format(x)+str(value)+'\n')
#        x+=20
#f.close()
