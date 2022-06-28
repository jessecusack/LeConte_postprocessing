library(ncdf4)

# Netcdf parameters
FillValue <- NaN
velocity_units <- "m s-1"
a_units <- "dB"  # Echo intensity units
g_units <- ""  # Percent good units
q_units <- ""  # Correlation units
time_units <- "s"
temperature_units <- "degree_C"
salinity_units <- ""  # Practical salinity is technically unitless (PSU)
turbidity_units <- ""  # Turbidity is technically unitless (NTU)
pressure_units <- "dbar"  # Pressure units
conductivity_units <- "S m-1"  # Conductivity units
allowed_coordinates <- c("enu", "xyz", "beam")  # Defined by OCE

# Metadata parameters
disallowed_metadata_names_adp <- c("units", "flags", "ensembleNumber", "ensembleInFile", "oceCoordinate", "oceBeamUnspreaded", "depthMean", "filename", "longitude", "latitude")
disallowed_metadata_names_ctd <- c("units", "flags", "filename", "longitude", "latitude", "ship", "scientist", 
                                   "institute", "address", "cruise", "date", "waterDepth", "header", "deploymentType", "dataNamesOriginal", "hexfilename", 
                                   "systemUploadTime", "startTime", "recoveryTime")
disallowed_metadata_names_rsk <- c("units", "flags", "filename", "longitude", "latitude", "dataNamesOriginal")


adp_write <- function(adp, file_path, overwrite=FALSE){
  # Writes data from an oce ADP object to a netcdf file.
  #
  # Parameters
  # ----------
  #    adp : An adp object.
  #    file_path : The file path, which should include the .nc suffix.
  #
  if (!inherits(adp, "adp")){
    stop("Must be an object of class adp")
  }

  if (file.exists(file_path) & !overwrite) {
    stop(paste("Write aborted. File", file_path, "exists and overwrite = False.", sep=" "))
  }
    
  oceCoord <- adp[['oceCoordinate']]
    
  if (!is.element(oceCoord, allowed_coordinates)) {
    stop(paste("Unknown coordinate system ", oceCoord, ", only '", paste(allowed_coordinates, collapse = ", "), "' are allowed.", sep=""))
  }

  # Logical checks
  g_exists <- !is.null(adp[["g"]])  # Percent good variable exists
  q_exists <- !is.null(adp[["q"]])  # Correlation variable exists
  is_sentinelV <- adp[['instrumentSubtype']] == "sentinelV"  # ADCP has a fifth vertical beam
  if (length(unique(adp[['salinity']])) == 1) {
    is_salinity <- FALSE  # Only one unique salinity value means it was not measured
  } else {
    is_salinity <- TRUE
  }
  if (length(unique(adp[['pressure']])) == 1) {
    is_pressure <- FALSE  # Only one unique pressure value means it was not measured
  } else {
    is_pressure <- TRUE
  }
  if (length(unique(adp[['temperature']])) == 1) {
    is_temperature <- FALSE  # Only one unique temperature value means it was not measured
  } else {
    is_temperature <- TRUE
  }

  # Create dimensions
  time <- adp[['time']]
  dist <- adp[['distance', 'numeric']]
  timedim <- ncdim_def("time", "seconds since 1970-01-01 00:00:00", as.double(time), longname="Time")
  depthdim <- ncdim_def("distance", "m", as.double(dist), longname="Instrument Z distance")

  # Define latitude and longitude
  lon_def <- ncvar_def("lon", "degree_east", list(), FillValue, "Longitude", prec="float")
  lat_def <- ncvar_def("lat", "degree_north", list(), FillValue, "Latitude", prec="float")
    
  # Define ensemble number
  ens_def <- ncvar_def("ensemble_number", "", list(timedim), FillValue, "Ensemble number", prec="integer")

  vars <- list(lon_def, lat_def, ens_def)

  # Define velocity variables
  if (oceCoord == 'beam') {
    v1_def <- ncvar_def("v1", velocity_units, list(timedim, depthdim), FillValue, "Beam 1 velocity", prec="float")
    v2_def <- ncvar_def("v2", velocity_units, list(timedim, depthdim), FillValue, "Beam 2 velocity", prec="float")
    v3_def <- ncvar_def("v3", velocity_units, list(timedim, depthdim), FillValue, "Beam 3 velocity", prec="float")
    v4_def <- ncvar_def("v4", velocity_units, list(timedim, depthdim), FillValue, "Beam 4 velocity", prec="float")
  } else if (oceCoord == 'enu') {
    v1_def <- ncvar_def("u", velocity_units, list(timedim, depthdim), FillValue, "Eastward velocity", prec="float")
    v2_def <- ncvar_def("v", velocity_units, list(timedim, depthdim), FillValue, "Northward velocity", prec="float")
    v3_def <- ncvar_def("w", velocity_units, list(timedim, depthdim), FillValue, "Vertical velocity", prec="float")
    v4_def <- ncvar_def("err", velocity_units, list(timedim, depthdim), FillValue, "Error velocity", prec="float")
  } else if (oceCoord == 'xyz') {
    v1_def <- ncvar_def("u", velocity_units, list(timedim, depthdim), FillValue, "Instrument X velocity", prec="float")
    v2_def <- ncvar_def("v", velocity_units, list(timedim, depthdim), FillValue, "Instrument Y velocity", prec="float")
    v3_def <- ncvar_def("w", velocity_units, list(timedim, depthdim), FillValue, "Instrument Z velocity", prec="float")
    v4_def <- ncvar_def("err", velocity_units, list(timedim, depthdim), FillValue, "Instrument error velocity", prec="float")
  }

  vars <- append(vars, list(v1_def, v2_def, v3_def, v4_def))

  # Define echo intensity variables
  a1_def <- ncvar_def("a1", a_units, list(timedim, depthdim), FillValue, "Beam 1 echo intensity", prec="float")
  a2_def <- ncvar_def("a2", a_units, list(timedim, depthdim), FillValue, "Beam 2 echo intensity", prec="float")
  a3_def <- ncvar_def("a3", a_units, list(timedim, depthdim), FillValue, "Beam 3 echo intensity", prec="float")
  a4_def <- ncvar_def("a4", a_units, list(timedim, depthdim), FillValue, "Beam 4 echo intensity", prec="float")

  vars <- append(vars, list(a1_def, a2_def, a3_def, a4_def))

  # Define percent good, if it exists
  if (g_exists) {
    g1_def <- ncvar_def("g1", g_units, list(timedim, depthdim), FillValue, "Beam 1 percent good", prec="float")
    g2_def <- ncvar_def("g2", g_units, list(timedim, depthdim), FillValue, "Beam 2 percent good", prec="float")
    g3_def <- ncvar_def("g3", g_units, list(timedim, depthdim), FillValue, "Beam 3 percent good", prec="float")
    g4_def <- ncvar_def("g4", g_units, list(timedim, depthdim), FillValue, "Beam 4 percent good", prec="float")

    vars <- append(vars, list(g1_def, g2_def, g3_def, g4_def))
  }

  # Define correlation, if it exists
  if (q_exists) {
    q1_def <- ncvar_def("q1", q_units, list(timedim, depthdim), FillValue, "Beam 1 correlation", prec="float")
    q2_def <- ncvar_def("q2", q_units, list(timedim, depthdim), FillValue, "Beam 2 correlation", prec="float")
    q3_def <- ncvar_def("q3", q_units, list(timedim, depthdim), FillValue, "Beam 3 correlation", prec="float")
    q4_def <- ncvar_def("q4", q_units, list(timedim, depthdim), FillValue, "Beam 4 correlation", prec="float")

    vars <- append(vars, list(q1_def, q2_def, q3_def, q4_def))
  }

  # Define 5th beam variables, if they exist
  if (is_sentinelV) {   
    if (oceCoord == 'beam') {
      vv_def <- ncvar_def("vv", velocity_units, list(timedim, depthdim), FillValue, "Beam 5 velocity", prec="float")
    } else if (oceCoord == 'enu') {
      vv_def <- ncvar_def("vv", velocity_units, list(timedim, depthdim), FillValue, "Vertical velocity", prec="float")
    } else if (oceCoord == 'xyz') {
      vv_def <- ncvar_def("vv", velocity_units, list(timedim, depthdim), FillValue, "Instrument Z velocity", prec="float")
    }
      
    va_def <- ncvar_def("va", a_units, list(timedim, depthdim), FillValue, "Beam 5 echo intensity", prec="float")

    vars <- append(vars, list(vv_def, va_def))

    if (g_exists) {
      vg_def <- ncvar_def("vg", g_units, list(timedim, depthdim), FillValue, "Beam 5 percent good", prec="float")
      vars <- append(vars, list(vg_def))
    }

    if (q_exists) {
      vq_def <- ncvar_def("vq", q_units, list(timedim, depthdim), FillValue, "Beam 5 correlation", prec="float")
      vars <- append(vars, list(vq_def))
    }

  }

  # Orientation variables
  pitch_def <- ncvar_def("pitch", "degree", list(timedim), FillValue, "Pitch", prec="float")
    # Confusingly, but deliberately, we name the variable 'rol' below because otherwise
    # we get conflicts with the xarray roll operation.
  roll_def <- ncvar_def("rol", "degree", list(timedim), FillValue, "Roll", prec="float")
  heading_def <- ncvar_def("heading", "degree", list(timedim), FillValue, "Heading", prec="float")

  vars <- append(vars, list(pitch_def, roll_def, heading_def))

  # Physical variables, if they exist
  if (is_temperature) {
    t_def <- ncvar_def("t", temperature_units, list(timedim), FillValue, "Temperature", prec="float")
    vars <- append(vars, list(t_def))
  }
  if (is_salinity) {
    SP_def <- ncvar_def("SP", salinity_units, list(timedim), FillValue, "Practical salinity", prec="float")
    vars <- append(vars, list(SP_def))
  }
  if (is_pressure) {
    p_def <- ncvar_def("p", pressure_units, list(timedim), FillValue, "Pressure", prec="float")
    vars <- append(vars, list(p_def))
  }

  # Logical check at the top looks for the case where file exists and overwrite is FALSE
  if (file.exists(file_path) & overwrite) {
    write(paste("Overwriting (deleting):", file_path, sep = " "), stdout())
    file.remove(file_path)
  }

  # Create file and directories if needed
  dir.create(file.path(dirname(file_path)))
  ncout <- nc_create(file_path, vars, force_v4 = TRUE)
    
  ncatt_put(ncout, "time", "standard_name", "time")
  # This no standard name for this variable
  # ncatt_put(ncout, depthdim, "standard_name", "distance_away_from_instrument")

  # Input latitude and longitude
  ncvar_put(ncout, lon_def, adp[['longitude']])
  ncvar_put(ncout, lat_def, adp[['latitude']])
  ncatt_put(ncout, lon_def, "standard_name", "longitude")
  ncatt_put(ncout, lat_def, "standard_name", "latitude")
    
  # Input ensemble number
  ncvar_put(ncout, ens_def, adp[['ensembleNumber']])

  # Input velocity
  ncvar_put(ncout, v1_def, adp[['v']][,,1])
  ncvar_put(ncout, v2_def, adp[['v']][,,2])
  ncvar_put(ncout, v3_def, adp[['v']][,,3])
  ncvar_put(ncout, v4_def, adp[['v']][,,4])
  if (oceCoord == 'beam') {
    ncatt_put(ncout, v1_def, "standard_name", "radial_sea_water_velocity_away_from_instrument")
    ncatt_put(ncout, v2_def, "standard_name", "radial_sea_water_velocity_away_from_instrument")
    ncatt_put(ncout, v3_def, "standard_name", "radial_sea_water_velocity_away_from_instrument")
    ncatt_put(ncout, v4_def, "standard_name", "radial_sea_water_velocity_away_from_instrument")
  } else if (oceCoord == 'enu') {
    ncatt_put(ncout, v1_def, "standard_name", "eastward_sea_water_velocity")
    ncatt_put(ncout, v2_def, "standard_name", "northward_sea_water_velocity")
    ncatt_put(ncout, v3_def, "standard_name", "upward_sea_water_velocity")
    ncatt_put(ncout, v4_def, "standard_name", "indicative_error_from_multibeam_acoustic_doppler_velocity_profiler_in_sea_water")
  }

  # Input echo intensity
  ncvar_put(ncout, a1_def, adp[['a', 'numeric']][,,1])
  ncvar_put(ncout, a2_def, adp[['a', 'numeric']][,,2])
  ncvar_put(ncout, a3_def, adp[['a', 'numeric']][,,3])
  ncvar_put(ncout, a4_def, adp[['a', 'numeric']][,,4])
  ncatt_put(ncout, a1_def, "standard_name", "signal_intensity_from_multibeam_acoustic_doppler_velocity_sensor_in_sea_water")
  ncatt_put(ncout, a2_def, "standard_name", "signal_intensity_from_multibeam_acoustic_doppler_velocity_sensor_in_sea_water")
  ncatt_put(ncout, a3_def, "standard_name", "signal_intensity_from_multibeam_acoustic_doppler_velocity_sensor_in_sea_water")
  ncatt_put(ncout, a4_def, "standard_name", "signal_intensity_from_multibeam_acoustic_doppler_velocity_sensor_in_sea_water")

  # Input percent good, if it exists
  if (g_exists) {
    ncvar_put(ncout, g1_def, adp[['g', 'numeric']][,,1])
    ncvar_put(ncout, g2_def, adp[['g', 'numeric']][,,2])
    ncvar_put(ncout, g3_def, adp[['g', 'numeric']][,,3])
    ncvar_put(ncout, g4_def, adp[['g', 'numeric']][,,4])
    ncatt_put(ncout, g1_def, "standard_name", "proportion_of_acceptable_signal_returns_from_acoustic_instrument_in_sea_water")
    ncatt_put(ncout, g2_def, "standard_name", "proportion_of_acceptable_signal_returns_from_acoustic_instrument_in_sea_water")
    ncatt_put(ncout, g3_def, "standard_name", "proportion_of_acceptable_signal_returns_from_acoustic_instrument_in_sea_water")
    ncatt_put(ncout, g4_def, "standard_name", "proportion_of_acceptable_signal_returns_from_acoustic_instrument_in_sea_water")
  }

  # Input correlation, if it exists
  if (q_exists) {
    ncvar_put(ncout, q1_def, adp[['q', 'numeric']][,,1])
    ncvar_put(ncout, q2_def, adp[['q', 'numeric']][,,2])
    ncvar_put(ncout, q3_def, adp[['q', 'numeric']][,,3])
    ncvar_put(ncout, q4_def, adp[['q', 'numeric']][,,4])
    ncatt_put(ncout, q1_def, "standard_name", "beam_consistency_indicator_from_multibeam_acoustic_doppler_velocity_profiler_in_sea_water")
    ncatt_put(ncout, q2_def, "standard_name", "beam_consistency_indicator_from_multibeam_acoustic_doppler_velocity_profiler_in_sea_water")
    ncatt_put(ncout, q3_def, "standard_name", "beam_consistency_indicator_from_multibeam_acoustic_doppler_velocity_profiler_in_sea_water")
    ncatt_put(ncout, q4_def, "standard_name", "beam_consistency_indicator_from_multibeam_acoustic_doppler_velocity_profiler_in_sea_water")
  }

  # Input 5th beam variables, if they exist
  if (is_sentinelV) {
    ncvar_put(ncout, vv_def, adp[['vv']])
    if (oceCoord == 'beam') {
      ncatt_put(ncout, vv_def, "standard_name", "radial_sea_water_velocity_away_from_instrument")
    } else if (oceCoord == 'enu') {
      ncatt_put(ncout, vv_def, "standard_name", "upward_sea_water_velocity")
    }
    ncvar_put(ncout, va_def, adp[['va', 'numeric']])
    ncatt_put(ncout, va_def, "standard_name", "signal_intensity_from_multibeam_acoustic_doppler_velocity_sensor_in_sea_water")
    if (g_exists) {
      ncvar_put(ncout, vg_def, adp[['vg', 'numeric']])
      ncatt_put(ncout, vg_def, "standard_name", "proportion_of_acceptable_signal_returns_from_acoustic_instrument_in_sea_water")
    }

    if (q_exists) {
      ncvar_put(ncout, vq_def, adp[['vq', 'numeric']])
      ncatt_put(ncout, vq_def, "standard_name", "beam_consistency_indicator_from_multibeam_acoustic_doppler_velocity_profiler_in_sea_water")
    }
  }

  # Input orientation variables
  ncvar_put(ncout, pitch_def, adp[['pitch']])
  ncatt_put(ncout, pitch_def, "standard_name", "platform_pitch")
  ncvar_put(ncout, roll_def, adp[['roll']])
  ncatt_put(ncout, roll_def, "standard_name", "platform_roll")
  ncvar_put(ncout, heading_def, adp[['heading']])
  ncatt_put(ncout, heading_def, "standard_name", "platform_orientation")

  # Input physical variables, if they exist
  if (is_temperature) {
    ncvar_put(ncout, t_def, adp[['temperature']])
    ncatt_put(ncout, t_def, "standard_name", "sea_water_temperature")
  }
  if (is_salinity) {
    ncvar_put(ncout, SP_def, adp[['salinity']])
    ncatt_put(ncout, SP_def, "standard_name", "sea_water_practical_salinity")
  }
  if (is_pressure) {
    ncvar_put(ncout, p_def, adp[['pressure']])
    ncatt_put(ncout, p_def, "standard_name", "sea_water_pressure")
  }

  # Copy metadata
  copy_metadata(adp@metadata, ncout, disallowed_metadata_names=disallowed_metadata_names_adp)

  write(paste("Saving:", file_path, "\n"), stdout())
  nc_close(ncout)

}


write_moored_ctd <- function(ctd, file_path, overwrite=FALSE) {
  # Writes data from an oce SBE CTD object obtained from a moored sensor to a
  # netcdf file.
  #
  # Parameters
  # ----------
  #    ctd : A SBE ctd object.
  #    file_path : The file path, which should include the .nc suffix.
  #
  if (!inherits(ctd, "ctd")){
    stop("Must be an object of class ctd")
  }

  if (file.exists(file_path) & !overwrite) {
    stop(paste("Write aborted. File", file_path, "exists and overwrite = False.", sep=" "))
  }
    
  data_names = names(ctd[["data"]])

  # Logical checks
  is_salinity <- "salinity" %in% data_names
  is_pressure <- "pressure" %in% data_names
  is_temperature <- "temperature" %in% data_names
  is_conductivity <- "conductivity" %in% data_names
  is_turbidity <- "turbidity" %in% data_names

  # Sometimes the time is a timestep rather than a POSIX time
  is_timestep <- is.null(ctd[['time']])
  if (is_timestep) {
    ctd[['time']] <- ctd[['timeS']] + ctd[['startTime']]
  }

  # Create dimensions
  time <- ctd[['time']]
  timedim <- ncdim_def("time", "seconds since 1970-01-01 00:00:00", as.double(time), longname="Time")

  # Define latitude and longitude
  lon_def <- ncvar_def("lon", "degree_east", list(), FillValue, "Longitude", prec="float")
  lat_def <- ncvar_def("lat", "degree_north", list(), FillValue, "Latitude", prec="float")

  vars <- list(lon_def, lat_def)

  # Physical variables, if they exist
  if (is_temperature) {
    t_def <- ncvar_def("t", temperature_units, list(timedim), FillValue, "Temperature", prec="float")
    vars <- append(vars, list(t_def))
  }
  if (is_salinity) {
    SP_def <- ncvar_def("SP", salinity_units, list(timedim), FillValue, "Practical salinity", prec="float")
    vars <- append(vars, list(SP_def))
  }
  if (is_pressure) {
    p_def <- ncvar_def("p", pressure_units, list(timedim), FillValue, "Pressure", prec="float")
    vars <- append(vars, list(p_def))
  }
  if (is_conductivity) {
    C_def <- ncvar_def("C", conductivity_units, list(timedim), FillValue, "Conductivity", prec="float")
    vars <- append(vars, list(C_def))
  }
  if (is_turbidity) {
    turb_def <- ncvar_def("turb", turbidity_units, list(timedim), FillValue, "Turbidity", prec="float")
    vars <- append(vars, list(turb_def))
  }

  # Logical check at the top looks for the case where file exists and overwrite is FALSE
  if (file.exists(file_path) & overwrite) {
    write(paste("Overwriting (deleting):", file_path, sep = " "), stdout())
    file.remove(file_path)
  }

  # Create file and directories if needed
  dir.create(file.path(dirname(file_path)))
  ncout <- nc_create(file_path, vars, force_v4 = TRUE)

  ncatt_put(ncout, "time", "standard_name", "time")

  # Input latitude and longitude
  ncvar_put(ncout, lon_def, ctd[['longitude']])
  ncvar_put(ncout, lat_def, ctd[['latitude']])
  ncatt_put(ncout, lon_def, "standard_name", "longitude")
  ncatt_put(ncout, lat_def, "standard_name", "latitude")

  # Input physical variables, if they exist
  if (is_temperature) {
    ncvar_put(ncout, t_def, ctd[['temperature']])
    ncatt_put(ncout, t_def, "standard_name", "sea_water_temperature")
  }
  if (is_salinity) {
    ncvar_put(ncout, SP_def, ctd[['salinity']])
    ncatt_put(ncout, SP_def, "standard_name", "sea_water_practical_salinity")
  }
  if (is_pressure) {
    ncvar_put(ncout, p_def, ctd[['pressure']])
    ncatt_put(ncout, p_def, "standard_name", "sea_water_pressure")
  }
  if (is_conductivity) {
    ncvar_put(ncout, C_def, ctd[['conductivity']])
    ncatt_put(ncout, C_def, "standard_name", "sea_water_electrical_conductivity")
  }
  if (is_turbidity) {
    ncvar_put(ncout, turb_def, ctd[['turbidity']])
    ncatt_put(ncout, turb_def, "standard_name", "sea_water_turbidity")
  }

  # Copy metadata
  copy_metadata(ctd@metadata, ncout, disallowed_metadata_names=disallowed_metadata_names_ctd)
    
  # Add header metadata, if it exists
  if ("header" %in% names(ctd@metadata)) {
    write("Adding header metadata", stdout())
    ncatt_put(ncout, 0, "fileHeader", paste(ctd[['header']], collapse = "\n"))
  }
              
  write(paste("Saving:", file_path, "\n"), stdout())
  nc_close(ncout)

}


write_moored_rsk <- function(rsk, file_path, overwrite=FALSE) {
  # Writes data from an oce RSK object obtained from a moored sensor to a
  # netcdf file.
  #
  # Parameters
  # ----------
  #    rsk : An rsk object.
  #    file_path : The file path, which should include the .nc suffix.
  #
  if (!inherits(rsk, "rsk")){
    stop("Must be an object of class rsk")
  }

  if (file.exists(file_path) & !overwrite) {
    stop(paste("Write aborted. File", file_path, "exists and overwrite = False.", sep=" "))
  }
    
  data_names = names(rsk[["data"]])

  # Logical checks
  is_salinity <- "salinity" %in% data_names
  is_pressure <- "pressure" %in% data_names
  is_temperature <- "temperature" %in% data_names
  is_conductivity <- "conductivity" %in% data_names
  is_turbidity <- "turbidity" %in% data_names

  # Create dimensions
  time <- rsk[['time']]
  timedim <- ncdim_def("time", "seconds since 1970-01-01 00:00:00", as.double(time), longname="Time")

  # Define latitude and longitude
  lon_def <- ncvar_def("lon", "degree_east", list(), FillValue, "Longitude", prec="float")
  lat_def <- ncvar_def("lat", "degree_north", list(), FillValue, "Latitude", prec="float")

  vars <- list(lon_def, lat_def)

  # Physical variables, if they exist
  if (is_temperature) {
    t_def <- ncvar_def("t", temperature_units, list(timedim), FillValue, "Temperature", prec="float")
    vars <- append(vars, list(t_def))
  }
  if (is_salinity) {
    SP_def <- ncvar_def("SP", salinity_units, list(timedim), FillValue, "Practical salinity", prec="float")
    vars <- append(vars, list(SP_def))
  }
  if (is_pressure) {
    p_def <- ncvar_def("p", pressure_units, list(timedim), FillValue, "Pressure", prec="float")
    vars <- append(vars, list(p_def))
  }
  if (is_conductivity) {
    C_def <- ncvar_def("C", conductivity_units, list(timedim), FillValue, "Conductivity", prec="float")
    vars <- append(vars, list(C_def))
  }
  if (is_turbidity) {
    turb_def <- ncvar_def("turb", turbidity_units, list(timedim), FillValue, "Turbidity", prec="float")
    vars <- append(vars, list(turb_def))
  }

  # Logical check at the top looks for the case where file exists and overwrite is FALSE
  if (file.exists(file_path) & overwrite) {
    write(paste("Overwriting (deleting):", file_path, sep = " "), stdout())
    file.remove(file_path)
  }

  # Create file and directories if needed
  dir.create(file.path(dirname(file_path)))
  ncout <- nc_create(file_path, vars, force_v4 = TRUE)

  ncatt_put(ncout, "time", "standard_name", "time")

  # Input latitude and longitude
  ncvar_put(ncout, lon_def, rsk[['longitude']])
  ncvar_put(ncout, lat_def, rsk[['latitude']])
  ncatt_put(ncout, lon_def, "standard_name", "longitude")
  ncatt_put(ncout, lat_def, "standard_name", "latitude")

  # Input physical variables, if they exist
  if (is_temperature) {
    ncvar_put(ncout, t_def, rsk[['temperature']])
    ncatt_put(ncout, t_def, "standard_name", "sea_water_temperature")
  }
  if (is_salinity) {
    ncvar_put(ncout, SP_def, rsk[['salinity']])
    ncatt_put(ncout, SP_def, "standard_name", "sea_water_practical_salinity")
  }
  if (is_pressure) {
    ncvar_put(ncout, p_def, rsk[['pressure']])
    ncatt_put(ncout, p_def, "standard_name", "sea_water_pressure")
  }
  if (is_conductivity) {
    ncvar_put(ncout, C_def, rsk[['conductivity']])
    ncatt_put(ncout, C_def, "standard_name", "sea_water_electrical_conductivity")
  }
  if (is_turbidity) {
    ncvar_put(ncout, turb_def, rsk[['turbidity']])
    ncatt_put(ncout, turb_def, "standard_name", "sea_water_turbidity")
  }

  # Copy metadata
  copy_metadata(rsk@metadata, ncout, disallowed_metadata_names=disallowed_metadata_names_rsk)
              
  write(paste("Saving:", file_path, "\n"), stdout())
  nc_close(ncout)

}


copy_metadata <- function(meta, ncout, namePrefix="", disallowed_metadata_names=c()){
  # Recursively copies oce object metadata to a netcdf file.
  #
  # Parameters
  # ----------
  #    meta : A named list of metadata.
  #    ncout : The netcdf file handle.
  #    namePrefix : When metadata object is a list of lists, a prefix is
  #                 appended the attribute name.
  #    disallowed_metadata_names : List of disallowed names.
  #
  if (namePrefix == "") {write("Adding metadata attributes", stdout())}

  for (name in names(meta)) {

    if (is.element(name, disallowed_metadata_names)) {
      write(paste("Skipping:", name), stdout())
      next
    }
      
    dat <- meta[[name]]

    if (is.list(dat)) {  # Recursive copy of other metadata
      write(paste("Stepping into:", name), stdout())
      copy_metadata(dat, ncout, paste(namePrefix, name, "_", sep=""), disallowed_metadata_names)
      next
    }

    if (class(dat)[1] == "matrix" | class(dat)[1] == "array") {
      write(paste("Skipping:", name), stdout())
      next
    }

    if (class(dat) == "logical") {
      dat <- as.integer(dat)
      write("Converting logical", stdout())
    }
      
    finalName <- paste(namePrefix, name, sep="")
    write(paste("Inputting:", name, "as", finalName), stdout())
    ncatt_put(ncout, 0, finalName, dat)
  }
}
