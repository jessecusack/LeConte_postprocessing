library(ncdf4)

v_units <- "m s-1"  # Velocity units
a_units <- "dB"  # Echo intensity units
g_units <- "percent"  # Percent good units
q_units <- ""  # Correlation units

adp_write <- function(adp, filePath, overwrite=FALSE){
  # Writes data from an oce ADP object to a netcdf file.
  #
  # Parameters
  # ----------
  #    adp : An adp object.
  #    filePath : The file path, which should include the .nc suffix.
  #
  if (!inherits(adp, "adp")){
    stop("Must be an object of class adp")
  }

  if (file.exists(filePath) & !overwrite) {
    stop(paste("Write aborted. File", filePath, "exists and overwrite = False.", sep=" "))
  }

  # Parameters
  FillValue <- NaN

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
  timedim <- ncdim_def("time", "POSIXct", as.double(time))    #time formatting FIX
  depthdim <- ncdim_def("distance", "m", as.double(dist))

  # Define latitude and longitude
  lon_def <- ncvar_def("lon", "degrees", list(), FillValue, "longitude_degrees_east", prec="float")
  lat_def <- ncvar_def("lat", "degrees", list(), FillValue, "latitude_degrees_north", prec="float")

  vars <- list(lon_def, lat_def)

  # Define velocity variables
  dlname <- "velocity_beam_1"
  v1_def <- ncvar_def("v1", v_units, list(timedim, depthdim), FillValue, dlname, prec="float")
  dlname <- "velocity_beam_2"
  v2_def <- ncvar_def("v2", v_units, list(timedim, depthdim), FillValue, dlname, prec="float")
  dlname <- "velocity_beam_3"
  v3_def <- ncvar_def("v3", v_units, list(timedim, depthdim), FillValue, dlname, prec="float")
  dlname <- "velocity_beam_4"
  v4_def <- ncvar_def("v4", v_units, list(timedim, depthdim), FillValue, dlname, prec="float")

  vars <- append(vars, list(v1_def, v2_def, v3_def, v4_def))

  # Define echo intensity variables
  dlname <- "echo_intensity_beam_1"
  a1_def <- ncvar_def("a1", a_units, list(timedim, depthdim), FillValue, dlname, prec="float")
  dlname <- "echo_intensity_beam_2"
  a2_def <- ncvar_def("a2", a_units, list(timedim, depthdim), FillValue, dlname, prec="float")
  dlname <- "echo_intensity_beam_3"
  a3_def <- ncvar_def("a3", a_units, list(timedim, depthdim), FillValue, dlname, prec="float")
  dlname <- "echo_intensity_beam_4"
  a4_def <- ncvar_def("a4", a_units, list(timedim, depthdim), FillValue, dlname, prec="float")

  vars <- append(vars, list(a1_def, a2_def, a3_def, a4_def))

  # Define percent good, if it exists
  if (g_exists) {
    dlname <- "percent_good_beam_1"
    g1_def <- ncvar_def("g1", g_units, list(timedim, depthdim), FillValue, dlname, prec="float")
    dlname <- "percent_good_beam_2"
    g2_def <- ncvar_def("g2", g_units, list(timedim, depthdim), FillValue, dlname, prec="float")
    dlname <- "percent_good_beam_3"
    g3_def <- ncvar_def("g3", g_units, list(timedim, depthdim), FillValue, dlname, prec="float")
    dlname <- "percent_good_beam_4"
    g4_def <- ncvar_def("g4", g_units, list(timedim, depthdim), FillValue, dlname, prec="float")

    vars <- append(vars, list(g1_def, g2_def, g3_def, g4_def))
  }

  # Define correlation, if it exists
  if (q_exists) {
    dlname <- "correlation_beam_1"
    q1_def <- ncvar_def("q1", q_units, list(timedim, depthdim), FillValue, dlname, prec="float")
    dlname <- "correlation_beam_2"
    q2_def <- ncvar_def("q2", q_units, list(timedim, depthdim), FillValue, dlname, prec="float")
    dlname <- "correlation_beam_3"
    q3_def <- ncvar_def("q3", q_units, list(timedim, depthdim), FillValue, dlname, prec="float")
    dlname <- "correlation_beam_4"
    q4_def <- ncvar_def("q4", q_units, list(timedim, depthdim), FillValue, dlname, prec="float")

    vars <- append(vars, list(q1_def, q2_def, q3_def, q4_def))
  }

  # Define 5th beam variables, if they exist
  if (is_sentinelV) {
    dlname <- "velocity_beam_5"
    vv_def <- ncvar_def("vv", v_units, list(timedim, depthdim), FillValue, dlname, prec="float")
    dlname <- "echo_intensity_beam_5"
    va_def <- ncvar_def("va", a_units, list(timedim, depthdim), FillValue, dlname, prec="float")

    vars <- append(vars, list(vv_def, va_def))

    if (g_exists) {
      dlname <- "percent_good_beam_5"
      vg_def <- ncvar_def("vg", g_units, list(timedim, depthdim), FillValue, dlname, prec="float")
      vars <- append(vars, list(vg_def))
    }

    if (q_exists) {
      dlname <- "correlation_beam_5"
      vq_def <- ncvar_def("vq", q_units, list(timedim, depthdim), FillValue, dlname, prec="float")
      vars <- append(vars, list(vq_def))
    }

  }

  # Orientation variables
  dlname <- "pitch"
  pitch_def <- ncvar_def("pitch", "degrees", list(timedim), FillValue, dlname, prec="float")
  dlname <- "roll"
    # Confusingly, but deliberately, we name the variable 'rol' below because otherwise
    # we get conflicts with the xarray roll operation.
  roll_def <- ncvar_def("rol", "degrees", list(timedim), FillValue, dlname, prec="float")
  dlname <- "heading"
  heading_def <- ncvar_def("heading", "degrees", list(timedim), FillValue, dlname, prec="float")

  vars <- append(vars, list(pitch_def, roll_def, heading_def))

  # Physical variables, if they exist
  if (is_temperature) {
    dlname <- "temperature"
    t_def <- ncvar_def("t", "deg C", list(timedim), FillValue, dlname, prec="float")
    vars <- append(vars, list(t_def))
  }
  if (is_salinity) {
    dlname <- "salinity"
    SP_def <- ncvar_def("SP", "PSU", list(timedim), FillValue, dlname, prec="float")
    vars <- append(vars, list(SP_def))
  }
  if (is_pressure) {
    dlname <- "pressure"
    p_def <- ncvar_def("p", "dbar", list(timedim), FillValue, dlname, prec="float")
    vars <- append(vars, list(p_def))
  }

  # Logical check at the top looks for the case where file exists and overwrite is FALSE
  if (file.exists(filePath) & overwrite) {
    write(paste("Overwriting (deleting):", filePath, sep = " "), stdout())
    file.remove(filePath)
  }

  # Create file
  ncout <- nc_create(filePath, vars, force_v4 = TRUE)

  # Input latitude and longitude
  ncvar_put(ncout, lon_def, adp[['longitude']])
  ncvar_put(ncout, lat_def, adp[['latitude']])

  # Input velocity
  ncvar_put(ncout, v1_def, adp[['v']][,,1])
  ncvar_put(ncout, v2_def, adp[['v']][,,2])
  ncvar_put(ncout, v3_def, adp[['v']][,,3])
  ncvar_put(ncout, v4_def, adp[['v']][,,4])

  # Input echo intensity
  ncvar_put(ncout, a1_def, adp[['a', 'numeric']][,,1])
  ncvar_put(ncout, a2_def, adp[['a', 'numeric']][,,2])
  ncvar_put(ncout, a3_def, adp[['a', 'numeric']][,,3])
  ncvar_put(ncout, a4_def, adp[['a', 'numeric']][,,4])

  # Input percent good, if it exists
  if (g_exists) {
    ncvar_put(ncout, g1_def, adp[['g', 'numeric']][,,1])
    ncvar_put(ncout, g2_def, adp[['g', 'numeric']][,,2])
    ncvar_put(ncout, g3_def, adp[['g', 'numeric']][,,3])
    ncvar_put(ncout, g4_def, adp[['g', 'numeric']][,,4])
  }

  # Input correlation, if it exists
  if (q_exists) {
    ncvar_put(ncout, q1_def, adp[['q', 'numeric']][,,1])
    ncvar_put(ncout, q2_def, adp[['q', 'numeric']][,,2])
    ncvar_put(ncout, q3_def, adp[['q', 'numeric']][,,3])
    ncvar_put(ncout, q4_def, adp[['q', 'numeric']][,,4])
  }

  # Input 5th beam variables, if they exist
  if (is_sentinelV) {
    ncvar_put(ncout, vv_def, adp[['vv']])
    ncvar_put(ncout, va_def, adp[['va', 'numeric']])
    if (g_exists) {
      ncvar_put(ncout, vg_def, adp[['vg', 'numeric']])
    }

    if (q_exists) {
      ncvar_put(ncout, vq_def, adp[['vq', 'numeric']])
    }
  }

  # Input orientation variables
  ncvar_put(ncout, pitch_def, adp[['pitch']])
  ncvar_put(ncout, roll_def, adp[['roll']])
  ncvar_put(ncout, heading_def, adp[['heading']])

  # Input physical variables, if they exist
  if (is_temperature) {
    ncvar_put(ncout, t_def, adp[['temperature']])
  }
  if (is_salinity) {
    ncvar_put(ncout, SP_def, adp[['salinity']])
  }
  if (is_pressure) {
    ncvar_put(ncout, p_def, adp[['pressure']])
  }

  # Copy metadata
  copy_metadata(adp@metadata, ncout)

  write(paste("Saving:", filePath, "\n"), stdout())
  nc_close(ncout)

}

enu_write <- function(adp, filePath, overwrite=FALSE){
  # Writes data from an oce ADP object to a netcdf file.
  #
  # Parameters
  # ----------
  #    adp : An adp object.
  #    filePath : The file path, which should include the .nc suffix.
  #
  if (!inherits(adp, "adp")){
    stop("Must be an object of class adp")
  }

  if (file.exists(filePath) & !overwrite) {
    stop(paste("Write aborted. File", filePath, "exists and overwrite = False.", sep=" "))
  }

  # Parameters
  FillValue <- NaN

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
  timedim <- ncdim_def("time", "POSIXct", as.double(time))    #time formatting FIX
  depthdim <- ncdim_def("distance", "m", as.double(dist))

  # Define latitude and longitude
  lon_def <- ncvar_def("lon", "degrees", list(), FillValue, "longitude_degrees_east", prec="float")
  lat_def <- ncvar_def("lat", "degrees", list(), FillValue, "latitude_degrees_north", prec="float")

  vars <- list(lon_def, lat_def)

  # Define velocity variables
  dlname <- "eastward_sea_water_velocity"
  u_def <- ncvar_def("u", v_units, list(timedim, depthdim), FillValue, dlname, prec="float")
  dlname <- "northward_sea_water_velocity"
  v_def <- ncvar_def("v", v_units, list(timedim, depthdim), FillValue, dlname, prec="float")
  dlname <- "upward_sea_water_velocity"
  w_def <- ncvar_def("w", v_units, list(timedim, depthdim), FillValue, dlname, prec="float")
  dlname <- "error_velocity_in_sea_water"
  e_def <- ncvar_def("err", v_units, list(timedim, depthdim), FillValue, dlname, prec="float")

  vars <- append(vars, list(u_def, v_def, w_def, e_def))

  # Define echo intensity variables
  dlname <- "echo_intensity_beam_1"
  a1_def <- ncvar_def("a1", a_units, list(timedim, depthdim), FillValue, dlname, prec="float")
  dlname <- "echo_intensity_beam_2"
  a2_def <- ncvar_def("a2", a_units, list(timedim, depthdim), FillValue, dlname, prec="float")
  dlname <- "echo_intensity_beam_3"
  a3_def <- ncvar_def("a3", a_units, list(timedim, depthdim), FillValue, dlname, prec="float")
  dlname <- "echo_intensity_beam_4"
  a4_def <- ncvar_def("a4", a_units, list(timedim, depthdim), FillValue, dlname, prec="float")

  vars <- append(vars, list(a1_def, a2_def, a3_def, a4_def))

  # Define percent good, if it exists
  if (g_exists) {
    dlname <- "percent_good_beam_1"
    g1_def <- ncvar_def("g1", g_units, list(timedim, depthdim), FillValue, dlname, prec="float")
    dlname <- "percent_good_beam_2"
    g2_def <- ncvar_def("g2", g_units, list(timedim, depthdim), FillValue, dlname, prec="float")
    dlname <- "percent_good_beam_3"
    g3_def <- ncvar_def("g3", g_units, list(timedim, depthdim), FillValue, dlname, prec="float")
    dlname <- "percent_good_beam_4"
    g4_def <- ncvar_def("g4", g_units, list(timedim, depthdim), FillValue, dlname, prec="float")

    vars <- append(vars, list(g1_def, g2_def, g3_def, g4_def))
  }

  # Define correlation, if it exists
  if (q_exists) {
    dlname <- "correlation_beam_1"
    q1_def <- ncvar_def("q1", q_units, list(timedim, depthdim), FillValue, dlname, prec="float")
    dlname <- "correlation_beam_2"
    q2_def <- ncvar_def("q2", q_units, list(timedim, depthdim), FillValue, dlname, prec="float")
    dlname <- "correlation_beam_3"
    q3_def <- ncvar_def("q3", q_units, list(timedim, depthdim), FillValue, dlname, prec="float")
    dlname <- "correlation_beam_4"
    q4_def <- ncvar_def("q4", q_units, list(timedim, depthdim), FillValue, dlname, prec="float")

    vars <- append(vars, list(q1_def, q2_def, q3_def, q4_def))
  }

  # Define 5th beam variables, if they exist
  if (is_sentinelV) {
    dlname <- "upward_sea_water_velocity_beam_5"
    vv_def <- ncvar_def("vv", v_units, list(timedim, depthdim), FillValue, dlname, prec="float")
    dlname <- "echo_intensity_beam_5"
    va_def <- ncvar_def("va", a_units, list(timedim, depthdim), FillValue, dlname, prec="float")

    vars <- append(vars, list(vv_def, va_def))

    if (g_exists) {
      dlname <- "percent_good_beam_5"
      vg_def <- ncvar_def("vg", g_units, list(timedim, depthdim), FillValue, dlname, prec="float")
      vars <- append(vars, list(vg_def))
    }

    if (q_exists) {
      dlname <- "correlation_beam_5"
      vq_def <- ncvar_def("vq", q_units, list(timedim, depthdim), FillValue, dlname, prec="float")
      vars <- append(vars, list(vq_def))
    }

  }

  # Orientation variables
  dlname <- "pitch"
  pitch_def <- ncvar_def("pitch", "degrees", list(timedim), FillValue, dlname, prec="float")
  dlname <- "roll"
    # Confusingly, but deliberately, we name the variable 'rol' below because otherwise
    # we get conflicts with the xarray roll operation.
  roll_def <- ncvar_def("rol", "degrees", list(timedim), FillValue, dlname, prec="float")
  dlname <- "heading"
  heading_def <- ncvar_def("heading", "degrees", list(timedim), FillValue, dlname, prec="float")

  vars <- append(vars, list(pitch_def, roll_def, heading_def))

  # Physical variables, if they exist
  if (is_temperature) {
    dlname <- "temperature"
    t_def <- ncvar_def("t", "deg C", list(timedim), FillValue, dlname, prec="float")
    vars <- append(vars, list(t_def))
  }
  if (is_salinity) {
    dlname <- "salinity"
    SP_def <- ncvar_def("SP", "PSU", list(timedim), FillValue, dlname, prec="float")
    vars <- append(vars, list(SP_def))
  }
  if (is_pressure) {
    dlname <- "pressure"
    p_def <- ncvar_def("p", "dbar", list(timedim), FillValue, dlname, prec="float")
    vars <- append(vars, list(p_def))
  }

  # Logical check at the top looks for the case where file exists and overwrite is FALSE
  if (file.exists(filePath) & overwrite) {
    write(paste("Overwriting (deleting):", filePath, sep = " "), stdout())
    file.remove(filePath)
  }

  # Create file
  ncout <- nc_create(filePath, vars, force_v4 = TRUE)

  # Input latitude and longitude
  ncvar_put(ncout, lon_def, adp[['longitude']])
  ncvar_put(ncout, lat_def, adp[['latitude']])

  # Input velocity
  ncvar_put(ncout, u_def, adp[['v']][,,1])
  ncvar_put(ncout, v_def, adp[['v']][,,2])
  ncvar_put(ncout, w_def, adp[['v']][,,3])
  ncvar_put(ncout, e_def, adp[['v']][,,4])

  # Input echo intensity
  ncvar_put(ncout, a1_def, adp[['a', 'numeric']][,,1])
  ncvar_put(ncout, a2_def, adp[['a', 'numeric']][,,2])
  ncvar_put(ncout, a3_def, adp[['a', 'numeric']][,,3])
  ncvar_put(ncout, a4_def, adp[['a', 'numeric']][,,4])

  # Input percent good, if it exists
  if (g_exists) {
    ncvar_put(ncout, g1_def, adp[['g', 'numeric']][,,1])
    ncvar_put(ncout, g2_def, adp[['g', 'numeric']][,,2])
    ncvar_put(ncout, g3_def, adp[['g', 'numeric']][,,3])
    ncvar_put(ncout, g4_def, adp[['g', 'numeric']][,,4])
  }

  # Input correlation, if it exists
  if (q_exists) {
    ncvar_put(ncout, q1_def, adp[['q', 'numeric']][,,1])
    ncvar_put(ncout, q2_def, adp[['q', 'numeric']][,,2])
    ncvar_put(ncout, q3_def, adp[['q', 'numeric']][,,3])
    ncvar_put(ncout, q4_def, adp[['q', 'numeric']][,,4])
  }

  # Input 5th beam variables, if they exist
  if (is_sentinelV) {
    ncvar_put(ncout, vv_def, adp[['vv']])
    ncvar_put(ncout, va_def, adp[['va', 'numeric']])
    if (g_exists) {
      ncvar_put(ncout, vg_def, adp[['vg', 'numeric']])
    }

    if (q_exists) {
      ncvar_put(ncout, vq_def, adp[['vq', 'numeric']])
    }
  }

  # Input orientation variables
  ncvar_put(ncout, pitch_def, adp[['pitch']])
  ncvar_put(ncout, roll_def, adp[['roll']])
  ncvar_put(ncout, heading_def, adp[['heading']])

  # Input physical variables, if they exist
  if (is_temperature) {
    ncvar_put(ncout, t_def, adp[['temperature']])
  }
  if (is_salinity) {
    ncvar_put(ncout, SP_def, adp[['salinity']])
  }
  if (is_pressure) {
    ncvar_put(ncout, p_def, adp[['pressure']])
  }

  # Copy metadata
  copy_metadata(adp@metadata, ncout)

  write(paste("Saving:", filePath, "\n"), stdout())
  nc_close(ncout)

}


write_moored_ctd <- function(ctd, filePath, overwrite=FALSE) {
  # Writes data from an oce CTD object obtained from a moored sensor to a
  # netcdf file.
  #
  # Parameters
  # ----------
  #    ctd : An ctd object.
  #    filePath : The file path, which should include the .nc suffix.
  #
  if (!inherits(ctd, "ctd") & !inherits(ctd, "rsk")){
    stop("Must be an object of class ctd")
  }

  if (file.exists(filePath) & !overwrite) {
    stop(paste("Write aborted. File", filePath, "exists and overwrite = False.", sep=" "))
  }

  # Parameters
  FillValue <- NaN
  t_units <- "deg C"  # Temperature units
  SP_units <- "PSU"  # Salinity units
  turb_units <- "NTU"  # Turbidity units
  p_units <- "dbar"  # Pressure units
  C_units <- "mS/cm"  # Conductivity units

  # Logical checks
  is_salinity <- !is.null(ctd[['salinity']])
  is_pressure <- !is.null(ctd[['pressure']])
  is_temperature <- !is.null(ctd[['temperature']])
  is_conductivity <- !is.null(ctd[['conductivity']])
  is_turbidity <- !is.null(ctd[['turbidity']])

  # Sometimes the time is a timestep rather than a POSIX time
  is_timestep <- is.null(ctd[['time']])
  if (is_timestep) {
    ctd[['time']] <- ctd[['timeS']] + ctd[['startTime']]
  }

  # Create dimensions
  time <- ctd[['time']]
  timedim <- ncdim_def("time", "POSIXct", as.double(time))    #time formatting FIX

  # Define latitude and longitude
  lon_def <- ncvar_def("lon", "degrees", list(), FillValue, "longitude_degrees_east", prec="float")
  lat_def <- ncvar_def("lat", "degrees", list(), FillValue, "latitude_degrees_north", prec="float")

  vars <- list(lon_def, lat_def)

  # Physical variables, if they exist
  if (is_temperature) {
    dlname <- "temperature"
    t_def <- ncvar_def("t", t_units, list(timedim), FillValue, dlname, prec="float")
    vars <- append(vars, list(t_def))
  }
  if (is_salinity) {
    dlname <- "salinity"
    SP_def <- ncvar_def("SP", SP_units, list(timedim), FillValue, dlname, prec="float")
    vars <- append(vars, list(SP_def))
  }
  if (is_pressure) {
    dlname <- "pressure"
    p_def <- ncvar_def("p", p_units, list(timedim), FillValue, dlname, prec="float")
    vars <- append(vars, list(p_def))
  }
  if (is_conductivity) {
    dlname <- "conductivity"
    C_def <- ncvar_def("C", C_units, list(timedim), FillValue, dlname, prec="float")
    vars <- append(vars, list(C_def))
  }
  if (is_turbidity) {
    dlname <- "turbidity"
    turb_def <- ncvar_def("turb", turb_units, list(timedim), FillValue, dlname, prec="float")
    vars <- append(vars, list(turb_def))
  }

  # Logical check at the top looks for the case where file exists and overwrite is FALSE
  if (file.exists(filePath) & overwrite) {
    write(paste("Overwriting (deleting):", filePath, sep = " "), stdout())
    file.remove(filePath)
  }

  # Create file
  ncout <- nc_create(filePath, vars, force_v4 = TRUE)

  # Input latitude and longitude
  ncvar_put(ncout, lon_def, ctd[['longitude']])
  ncvar_put(ncout, lat_def, ctd[['latitude']])

  # Input physical variables, if they exist
  if (is_temperature) {
    ncvar_put(ncout, t_def, ctd[['temperature']])
  }
  if (is_salinity) {
    ncvar_put(ncout, SP_def, ctd[['salinity']])
  }
  if (is_pressure) {
    ncvar_put(ncout, p_def, ctd[['pressure']])
  }
  if (is_conductivity) {
    ncvar_put(ncout, C_def, ctd[['conductivity']])
  }
  if (is_turbidity) {
    ncvar_put(ncout, turb_def, ctd[['turbidity']])
  }

  # Copy metadata
  copy_metadata(ctd@metadata, ncout)

  write(paste("Saving:", filePath, "\n"), stdout())
  nc_close(ncout)

}


copy_metadata <- function(meta, ncout, namePrefix=""){
  # Recursively copies oce object metadata to a netcdf file.
  #
  # Parameters
  # ----------
  #    meta : A named list of metadata.
  #    ncout : The netcdf file handle.
  #    namePrefix : When metadata object is a list of lists, a prefix is
  #                 appended the attribute name.
  #
  if (namePrefix == "") {write("Adding metadata attributes", stdout())}

  for (name in names(meta)) {

    dat <- meta[[name]]

    is_units_or_flags <- (name == "units") | (name == "flags")

    if (is_units_or_flags) {
      write(paste("Skipping:", name), stdout())
      next
    }

    if (is.list(dat)) {  # Recursive copy of other metadata
      write(paste("Stepping into:", name))
      copy_metadata(dat, ncout, paste(namePrefix, name, sep=""), stdout())
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
