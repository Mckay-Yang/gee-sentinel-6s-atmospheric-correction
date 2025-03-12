import ee
from Py6S import *
import datetime
import math
import os
import sys
sys.path.append(os.path.join(os.path.dirname(os.getcwd()),"bin"))
from atmospheric import Atmospheric

ee.Authenticate()
ee.Initialize(project="ee-yangluhao990714")

def py6s_main(S2: ee.Image, index: int):
    # S2 = ee.Image(
    # ee.ImageCollection("COPERNICUS/S2")
    #     .filterBounds(geom)
    #     .filterDate(date, date.advance(3, "month"))
    #     .sort("system:time_start")
    #     .first()
    # )

    geom = ee.Geometry(S2.geometry().centroid())
    # top of atmosphere reflectance
    toa = S2.divide(10000)

    info = S2.getInfo()["properties"]
    scene_date = datetime.datetime.fromtimestamp(info["system:time_start"] / 1000, tz = datetime.timezone.utc)# i.e. Python uses seconds, EE uses milliseconds
    date = ee.Date(info["system:time_start"])
    solar_z = info["MEAN_SOLAR_ZENITH_ANGLE"]

    h2o = Atmospheric.water(geom,date).getInfo()
    o3 = Atmospheric.ozone(geom,date).getInfo()
    aot = Atmospheric.aerosol(geom,date).getInfo()

    SRTM = ee.Image("CGIAR/SRTM90_V4") # Shuttle Radar Topography mission covers *most* of the Earth
    alt = SRTM.reduceRegion(reducer = ee.Reducer.mean(), geometry = geom.centroid()).get("elevation").getInfo()
    km = alt / 1000 # i.e. Py6S uses units of kilometers

    # Instantiate
    s = SixS()

    # Atmospheric constituents
    s.atmos_profile = AtmosProfile.UserWaterAndOzone(h2o, o3)
    s.aero_profile = AeroProfile.Continental
    s.aot550 = aot

    # Earth-Sun-satellite geometry
    s.geometry = Geometry.User()
    s.geometry.view_z = 0               # always NADIR (I think..)
    s.geometry.solar_z = solar_z        # solar zenith angle
    s.geometry.month = scene_date.month # month and day used for Earth-Sun distance
    s.geometry.day = scene_date.day     # month and day used for Earth-Sun distance
    s.altitudes.set_sensor_satellite_level()
    s.altitudes.set_target_custom_altitude(km)


    def spectralResponseFunction(bandname):
        """
        Extract spectral response function for given band name
        """
        bandSelect = {
            "B1":PredefinedWavelengths.S2A_MSI_01,
            "B2":PredefinedWavelengths.S2A_MSI_02,
            "B3":PredefinedWavelengths.S2A_MSI_03,
            "B4":PredefinedWavelengths.S2A_MSI_04,
            "B5":PredefinedWavelengths.S2A_MSI_05,
            "B6":PredefinedWavelengths.S2A_MSI_06,
            "B7":PredefinedWavelengths.S2A_MSI_07,
            "B8":PredefinedWavelengths.S2A_MSI_08,
            "B8A":PredefinedWavelengths.S2A_MSI_8A,
            "B9":PredefinedWavelengths.S2A_MSI_09,
            "B10":PredefinedWavelengths.S2A_MSI_10,
            "B11":PredefinedWavelengths.S2A_MSI_11,
            "B12":PredefinedWavelengths.S2A_MSI_12,
            }
        return Wavelength(bandSelect[bandname])


    def toa_to_rad(bandname):
        """
        Converts top of atmosphere reflectance to at-sensor radiance
        """

        # solar exoatmospheric spectral irradiance
        ESUN = info["SOLAR_IRRADIANCE_" + bandname]
        solar_angle_correction = math.cos(math.radians(solar_z))

        # Earth-Sun distance (from day of year)
        doy = scene_date.timetuple().tm_yday
        d = 1 - 0.01672 * math.cos(0.9856 * (doy - 4))# http://physics.stackexchange.com/questions/177949/earth-sun-distance-on-a-given-day-of-the-year

        # conversion factor
        multiplier = ESUN * solar_angle_correction/(math.pi * d ** 2)

        # at-sensor radiance
        rad = toa.select(bandname).multiply(multiplier)

        return rad


    def surface_reflectance(bandname):
        """
        Calculate surface reflectance from at-sensor radiance given waveband name
        """

        # run 6S for this waveband
        s.wavelength = spectralResponseFunction(bandname)
        s.run()

        # extract 6S outputs
        Edir = s.outputs.direct_solar_irradiance             #direct solar irradiance
        Edif = s.outputs.diffuse_solar_irradiance            #diffuse solar irradiance
        Lp   = s.outputs.atmospheric_intrinsic_radiance      #path radiance
        absorb  = s.outputs.trans["global_gas"].upward       #absorption transmissivity
        scatter = s.outputs.trans["total_scattering"].upward #scattering transmissivity
        tau2 = absorb * scatter                                #total transmissivity

        # radiance to surface reflectance
        rad = toa_to_rad(bandname)
        ref = rad.subtract(Lp).multiply(math.pi).divide(tau2 * (Edir + Edif))

        return ref


    # all wavebands
    output = S2.select("QA60").double()
    for band in ["B1", "B2", "B3", "B4", "B5", "B6", "B7", "B8", "B8A", "B9", "B10", "B11", "B12"]:
        output = output.addBands(surface_reflectance(band).double())

    region = geom.buffer(5000).bounds().getInfo()["coordinates"]
    channels = ["B4","B3","B2"]

    # set some properties for export
    dateString = scene_date.strftime("%Y-%m-%d-%H-%M-%S-%f")
    ref = output.set({"satellite":"Sentinel 2",
                "fileID":info["system:index"],
                "date":dateString,
                "aerosol_optical_thickness":aot,
                "water_vapour":h2o,
                "ozone":o3,
                "system:time_start":info["system:time_start"],
                "system:time_end":info["system:time_end"]}
                )

    assetID = "projects/ee-yangluhao990714/assets/sentinel2_l1c_atmcorr/image" + str(info["system:time_start"]) + str(index)

    # export
    export = ee.batch.Export.image.toAsset(
        image = ref,
        description = "sentinel2_atmcorr_export" + str(index),
        # folder = "projects/ee-yangluhao990714/assets/sentinel2_l1c_atmcorr",
        # fileNamePrefix = "image"+str(info["system:time_start"]) + str(index),
        crs = "EPSG:4326",
        assetId = assetID,
        region = region,
        scale = 10
    )

    # uncomment to run the export
    export.start()

start_date = ee.Date("2015-01-01T00:00:00")
end_date = ee.Date("2016-01-01T00:00:00")

aoi = ee.FeatureCollection("projects/ee-yangluhao990714/assets/downstream_aoi")
sentinel2_l1c_col = ee.ImageCollection("COPERNICUS/S2_HARMONIZED")\
    .filterBounds(aoi).filterDate(start_date, end_date)

sentinel2_l1c_col_list = sentinel2_l1c_col.toList(sentinel2_l1c_col.size())
length = sentinel2_l1c_col.size().getInfo()
print("Total number of images:", length)
for i in range(length):
    if i == 0:
        S2 = ee.Image(sentinel2_l1c_col_list.get(i))
        py6s_main(S2, i)
        print("Exporting image", i, "done.")
