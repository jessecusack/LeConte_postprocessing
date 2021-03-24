#' @title Net CDF Export,
#'
#' @description Exports an adp object to a netcdf using variables and metadata within adp.
#'
#'
#' @param obj an adp object from the oce class
#' @param name name of the NetCDF file to be produced


##nc metadata
create_nc_from_adp <- function(adp, name){
  if (!inherits(adp, "adp")){
    stop("method is only for objects of class '", "adp", "'")
  }

  #file name and path
  ncpath <- "./"
  ncname <- name
  ncfname <- paste(ncpath, ncname, ".nc", sep = "")
  
  ####setting dimensions and definitions####
  #dimension variables from adp object
  time <- adp[['time']]
  dist <- adp[['distance', 'numeric']]
  lon <- adp[['longitude']]
  lat <- adp[['latitude']]
  
  #create dimensions
  timedim <- ncdim_def("time", "POSIXct", as.double(time))    #time formatting FIX
  depthdim <- ncdim_def("distance", "m", as.double(dist))
  londim <- ncdim_def("lon", "degrees_east" , as.double(lon))
  latdim <- ncdim_def("lat", "degrees_north", as.double(lat))
  
  #set fill value
  FillValue <- 1e35
  
  #define variables
  dlname <- "eastward_sea_water_velocity"
  u_def <- ncvar_def("u", "m/s", list(timedim, depthdim, londim, latdim), FillValue, dlname, prec = "float")
  
  dlname <- "northward_sea_water_velocity"
  v_def <- ncvar_def("v", "m/s", list(timedim, depthdim, londim, latdim), FillValue, dlname, prec = "float")
  
  dlname <- "upward_sea_water_velocity"
  w_def <- ncvar_def("w", "m/s", list(timedim, depthdim, londim, latdim), FillValue, dlname, prec = "float")
  
  dlname <- "error_velocity_in_sea_water"
  e_def <- ncvar_def("err", "m/s", list(timedim, depthdim, londim, latdim), FillValue, dlname, prec = "float")
  
  dlname <- "ADCP_echo_intensity_beam_1"
  b1_def <- ncvar_def("ei1", "a", list(timedim, depthdim, londim, latdim), FillValue, dlname, prec = "float")
  
  dlname <- "ADCP_echo_intensity_beam_2"
  b2_def <- ncvar_def("ei2", "a", list(timedim, depthdim, londim, latdim), FillValue, dlname, prec = "float")
  
  dlname <- "ADCP_echo_intensity_beam_3"
  b3_def <- ncvar_def("ei3", "a", list(timedim, depthdim, londim, latdim), FillValue, dlname, prec = "float")
  
  dlname <- "ADCP_echo_intensity_beam_4"
  b4_def <- ncvar_def("ei4", "a", list(timedim, depthdim, londim, latdim), FillValue, dlname, prec = "float")
  
  dlname <- "percent_good_beam_1"
  pg1_def <- ncvar_def("pg1", "counts", list(timedim, depthdim, londim, latdim), FillValue, dlname, prec = "float")
  
  dlname <- "percent_good_beam_2"
  pg2_def <- ncvar_def("pg2", "counts", list(timedim, depthdim, londim, latdim), FillValue, dlname, prec = "float")
  
  dlname <- "percent_good_beam_3"
  pg3_def <- ncvar_def("pg3", "counts", list(timedim, depthdim, londim, latdim), FillValue, dlname, prec = "float")
  
  dlname <- "percent_good_beam_4"
  pg4_def <- ncvar_def("pg4", "counts", list(timedim, depthdim, londim, latdim), FillValue, dlname, prec = "float")
  
  dlname <- "pitch"
  p_def <- ncvar_def("pitch", "degrees", list(timedim, londim, latdim), FillValue, dlname, prec = "float")
  
  dlname <- "roll"
  r_def <- ncvar_def("roll", "degrees", list(timedim, londim, latdim), FillValue, dlname, prec = "float")
  
  dlname <- "heading"
  head_def <- ncvar_def("heading", "degrees", list(timedim, londim, latdim), FillValue, dlname, prec = "float")
  
  dlname <- "temperature"
  t_def <- ncvar_def("t", "deg C", list(timedim, londim, latdim), FillValue, dlname, prec = "float")
  
  dlname <- "pressure"
  pr_def <- ncvar_def("p", "dbar", list(timedim, londim, latdim), FillValue, dlname, prec = "float")
  
  dlname <- "instrument depth"
  D_def <- ncvar_def("D", "m", list(timedim, londim, latdim), FillValue, dlname, prec = "float")
  
  #write out definitions to new nc file
  ncout <- nc_create(ncfname, list(u_def, v_def, w_def, e_def, b1_def, b2_def, 
                                   b3_def, b4_def, pg1_def, pg2_def, pg3_def, 
                                   pg4_def, p_def, r_def, head_def, t_def, pr_def, 
                                   D_def), force_v4 = TRUE)
  
  # Insert variables into nc file
  ncvar_put(ncout, u_def, adp[['v']][,,1])
  ncvar_put(ncout, v_def, adp[['v']][,,2])
  ncvar_put(ncout, w_def, adp[['v']][,,3])
  ncvar_put(ncout, e_def, adp[['v']][,,4])
  ncvar_put(ncout, b1_def, adp[['a', 'numeric']][,,1])
  ncvar_put(ncout, b2_def, adp[['a', 'numeric']][,,2])
  ncvar_put(ncout, b3_def, adp[['a', 'numeric']][,,3])
  ncvar_put(ncout, b4_def, adp[['a', 'numeric']][,,4])
  ncvar_put(ncout, pg1_def, adp[['g', 'numeric']][,,1])
  ncvar_put(ncout, pg2_def, adp[['g', 'numeric']][,,2])
  ncvar_put(ncout, pg3_def, adp[['g', 'numeric']][,,3])
  ncvar_put(ncout, pg4_def, adp[['g', 'numeric']][,,4])
  ncvar_put(ncout, p_def, adp[['pitch']])
  ncvar_put(ncout, r_def, adp[['roll']])
  ncvar_put(ncout, head_def, adp[['heading']])
  ncvar_put(ncout, t_def, adp[['temperature']])
  ncvar_put(ncout, pr_def, adp[['pressure']])
  ncvar_put(ncout, D_def, adp[['depth']])
  
  # Metadata
  ncatt_put(ncout, 0, "Deployment_date", adp[['deploymentTime']])
  ncatt_put(ncout, 0, "Recovery_date", adp[['recoveryTime']])
  ncatt_put(ncout, 0, "firmware_version", adp[['firmwareVersion']])
  ncatt_put(ncout, 0, "frequency", adp[['frequency']])
  ncatt_put(ncout, 0, "beam_pattern", adp[['beamPattern']])
  ncatt_put(ncout, 0, "janus", adp[['numberOfBeams']])
  ncatt_put(ncout, 0, "instrument_type", adp[['instrumentSubtype']])
  ncatt_put(ncout, 0, "pings_per_ensemble", adp[['pingsPerEnsemble']])
  ncatt_put(ncout, 0, "valid_correlation_range", adp[['lowCorrThresh']])
  ncatt_put(ncout, 0, "minmax_percent_good", adp[['percentGdMinimum']])
  ncatt_put(ncout, 0,"minmax_percent_good", "100")
  ncatt_put(ncout, 0, "error_velocity_threshold", adp[['errorVelocityMaximum']])
  ncatt_put(ncout, 0, "transmit_pulse_length_cm", adp[['xmitPulseLength']])
  ncatt_put(ncout, 0, "false_target_reject_values", adp[['falseTargetThresh']])
  ncatt_put(ncout, 0, "ADCP_serial_number", adp[['serialNumber']])
  ncatt_put(ncout, 0, "transform", adp[['oceCoordinate']])
  ncatt_put(ncout, 0, "DATA_TYPE", adp[['instrumentType']])
  ncatt_put(ncout, 0, "COORD_SYSTEM", adp[['oceCoordinate']])
  ncatt_put(ncout, 0, "longitude", adp[['longitude']])
  ncatt_put(ncout, 0, "latitude", adp[['latitude']])
  ncatt_put(ncout, 0, "magnetic_variation", adp[['magneticVariation']])
  
  ncatt_put(ncout, 0, "pred_accuracy", adp[['velocityResolution']])
  
  ncatt_put(ncout, "lat", "standard_name", "latitude")
  ncatt_put(ncout, "lon", "standard_name", "longitude")
  ncatt_put(ncout, "D", "standard_name", "depth")
  ncatt_put(ncout, "t", "standard_name", "")
  
  ncatt_put(ncout, "distance", "positive", adp[['orientation']])     #direction of depth axis
  ncatt_put(ncout, "distance", "axis", "y")
  ncatt_put(ncout, "time", "axis", "x")
  
  ncatt_put(ncout, "t", "valid_min", min(adp[['temperature']]))
  ncatt_put(ncout, "t", "valid_range", "")
  ncatt_put(ncout, "t", "valid_max", max(adp[['temperature']]))
  
  nc_close(ncout)
  
}