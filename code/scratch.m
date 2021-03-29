clear all, close all

adcp_beam = parse_adcp('~/Downloads/SN_54000.000');
declination = 19.3;
angles = [-adcp_beam.roll(:), -adcp_beam.pitch(:), adcp_beam.heading(:) - declination];
beam2earth = quaternion(angles,'eulerd','YXZ','frame');
adcp_earth = adcp_beam2earth(adcp_beam,beam2earth,'up');
