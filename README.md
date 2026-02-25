# SLAP_Rocket_Example
Rocket finite element model case study for dynamic environment testing. This was first used in "Damage Metric-Based Vibration Testing of a Rocket Component" by Behling, Allen, Bahr, Mayes, & Delima as a case study on which to apply the scaled lab PSD (SLAP) method. This repository contains the files used for the case studies in the paper, as well as code used to simulate the case studies. The component of interest is the removable component from the BARC system, a commonly used case study in the dynamic environments case study. This was chosen because of its common usage and because it has multiple elastic modes below 2000 Hz, making it an interesting case study. 

There is a lot of content in this repository, but the most important is "SLAPScript.m" which was used to simulate the tests whose results are shown in the paper. The folders can be broken down as follows:

Environment: contains flight environment files.

  FlightForces.mat - forcing functions used to create environment
  
  FlightStressPSD.mat - maximum RMS Von Mises stress PSD in flight environment
  
  RocketEnvironment.mat - available at link below because too big for github. Flight environment PSD matrix.

FRFs: contains frequency response function matrices used to simulate tests as well as some files needed to create them.
  FlightAccelNodes.mat - nodes of rocket FEM corresponding to accelerometer locations on baseplate and DUT
  FlightFRF.mat - FRF used to create flight environment (given corresponding forcing functions in FlightForces.mat)
  FlightForceNodes.mat - ditto FlightAccelNodes.mat but with flight force locations
  LabAccelNodes.mat - nodes of barc + baseplate (lab config) FEM where accelerometers are
  LabFRF.mat - FRF used to control environments in simulated tests
  LabShakerNodes.mat - ditto LabAccelNodes with shaker locs
  LabStressFRFs.mat and FlightStressFRFs.mat - relate lab / flight forces to stress at stress element locs. (too big for github - see link below)

Functions: contains functions used in simulations
  GetFatigueRatio.m - given two Von Mises stress PSDs (flight and lab), calculates ratio of fatigue damage in them
  GetModalRMSRatio.m - given RMS modal accelerations, calculates smallest ratio (lab / flight) which controls acc. scaling
  GetPeakStressRatio.m - calculate expected peak stress ratio (lab / flight) given Von Mises stress PSDs
  GetRMSStressRatio.m - same but with RMS stress
  GetStressFunc.m - uses stress FRFs to calculate stress PSDs in flight env / lab test
  GetStressPSDs.m - calculates approximate stress PSDs (eqs. 2 and 3 of paper) given fixed-base modal displacement flight and lab PSDs
  ModalFilterFunc.m - modal filters physical acceleration PSDs
  SLAPfunc.m - implements SLAP method
  cl_model.m - calculates FRF given mode shapes
  get_psd.m - retains diagonal portions of PSD matrix

ModeShapes: contains FEM and reduced sets of modes
  BARCAccelModes.mat - contains mode shapes at BARC accel locations
  BARC_BaseplateModes.mat - FEM modes of lab config (BARC + Baseplate)
  BARC_FixedBaseModes.mat - FEM modes of BARC with fixed boundary condition
  BARC_FixedBaseStressModes.csv - ditto but with stress instead of displacement
  FixedBaseStressModes.mat - above file but sorted into .mat in nicer format
  FlightStressModes.csv - stress modes from abaqus of flight config
  FlightStressModes.csv - sorted into .mat file
  LabStressModes.mat - stress modes of lab config
  FullRocketModes.mat - FEM modes of rocket (too big for github - find at link below)

Results: contains results of simulations from paper
  Results_SLAP_Buzz.mat - results from applying SLAP-Buzz
  Results_SLAP_Buzz_Elastic.mat - results from applying SLAP Buzz excluding modal acceleration metric
  Results_SLAP_Control.mat - results from applying SLAP-Control
  Results_SLAP_Control_4Shakers.mat - results from applying SLAP-Control to 4 shaker case
  Results_SLAP_Control_Elastic.mat - results from applying SLAP-Control to 6 shaker case excluding modal acceleration metric

Scripts: code used to make FRFs, simulate tests, etc.
  FlightEnvironmentScript.m - used to make flight environment
  FlightFRFScript.m - makes FRF relating flight forces to acceleration
  GetAbqModesScript.m - loads modes from Abaqus .inp files into .mat
  GetLabLocationsScript.m - script used to select accel and shaker locations of interest for lab test setup
  GetModesAtAccels.m - script used to retain modes at BARC accels only (for modal filtering)
  GetRocketLocationsScript.m - ditto LabLocationsScript but with flight config
  LabFRFScript.m - make lab control FRF
  ModalFilterScript.m - check how well modal filtering works
  SimulateTestScript.m - simulate MIMO test for initial environment
  StressFRFScript.m - makes FRF relating shaker / lab forces to stress at locations of interest
  StressModesScript.m - converts .csv stress modes output from Abaqus to .mat files for use elsewhere


A few files were too large to include on github, and these are available at: https://byu.box.com/s/xanu0aqfyxedwmop60q8fpvicxrclrxr.


