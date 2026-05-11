# NOTE: We are in currently updating the FEM of the rocket to fix some issues, and also greatly expanding the documentation. Watch for an update by the end of May 2026
# SLAP Rocket Example

This repository contains files and code for a rocket finite element model (FEM) case study on dynamic environment testing, first presented in *"Damage Metric-Based Vibration Testing of a Rocket Component"* by Behling, Allen, Bahr, Mayes, & DeLima. The case study applies the **Scaled Lab PSD (SLAP)** method to the removable component from the BARC system. This component was selected for its common usage and because it has multiple elastic modes below 2000 Hz, making it a rich test case.

> **Note:** A few files were too large for GitHub and are hosted externally. Download them at: https://byu.box.com/s/xanu0aqfyxedwmop60q8fpvicxrclrxr. After the data has been downloaded, unzip the file in the "LargeFiles" folder. There is a text document inside the folder that explains how the file path should be organized.

---

## Getting Started

The most important file in this repository is **`Scripts/SLAPScript.m`**, which was used to simulate the tests whose results appear in the paper. Start there to reproduce the paper's findings.

---

## Repository Structure

### AbaqusFiles
Contains `.cae` and `.inp` files used to extract stress and displacement modes for flight, lab, and fixed-base configurations. Also contains a python script that was used to translate infromation from the Abaqus *.odb files into *.mat files.

---

### Environment
Flight environment files used to define loading conditions.

| File | Description |
|------|-------------|
| `FlightForces.mat` | Forcing functions used to create the flight environment |
| `FlightStressPSD.mat` | Maximum RMS Von Mises stress PSD in the flight environment |

---

### FRFs
Frequency response function (FRF) matrices used to simulate tests, along with supporting files needed to generate them.

| File | Description |
|------|-------------|
| `FlightAccelNodes.mat` | Rocket FEM nodes corresponding to accelerometer locations on the baseplate and DUT |
| `FlightForceNodes.mat` | Same as above, but for flight force locations |
| `FlightFRF.mat` | FRF used to create the flight environment (pairs with `FlightForces.mat`) |
| `LabAccelNodes.mat` | BARC + baseplate (lab config) FEM nodes at accelerometer locations |
| `LabShakerNodes.mat` | Same as above, but for shaker locations |
| `LabFRF.mat` | FRF used to control environments in simulated lab tests |

---

### Functions
MATLAB functions called during simulations.

| File | Description |
|------|-------------|
| `SLAPfunc.m` | Implements the SLAP method |
| `GetFatigueRatio.m` | Given two Von Mises stress PSDs (flight and lab), calculates the ratio of fatigue damage |
| `GetModalRMSRatio.m` | Given RMS modal accelerations, calculates the smallest lab/flight ratio controlling acceleration scaling |
| `GetPeakStressRatio.m` | Calculates expected peak stress ratio (lab/flight) from Von Mises stress PSDs |
| `GetRMSStressRatio.m` | Same as above, but using RMS stress |
| `GetStressFunc.m` | Uses stress FRFs to calculate stress PSDs in the flight environment and lab test |
| `GetStressPSDs.m` | Calculates approximate stress PSDs (Eqs. 2–3 in paper) from fixed-base modal displacement PSDs |
| `ModalFilterFunc.m` | Applies modal filters to physical acceleration PSDs |
| `cl_model.m` | Calculates an FRF from mode shapes |
| `get_psd.m` | Retains the diagonal portions of a PSD matrix |

---

### LargeFiles
Files that are too large to be stored on GitHub *(see link above)*

| File | Description |
|------|-------------|
| `FlightStressFRFs.mat` | Relates flight forces to stress at stress element locations |
| `Full_Rocket_Modes.mat` | FEM modes of the full rocket |
| `FullRocketUpdated.cae` | FEM
| `LabStressFRFs.mat` | Relates lab forces to stress at stress element locations |
| `RocketEnvironment.mat` | Flight environment PSD matrix |

---

### ModeShapes
FEM mode shapes and reduced mode sets for various configurations.

| File | Description |
|------|-------------|
| `BARC_BaseplateModes.mat` | FEM modes of the lab configuration (BARC + baseplate) |
| `BARC_FixedBaseModes.mat` | FEM modes of BARC with a fixed boundary condition |
| `BARC_FixedBaseStressModes.csv` | Same as above, but with stress instead of displacement |
| `FixedBaseStressModes.mat` | `BARC_FixedBaseStressModes.csv` reorganized into `.mat` format |
| `BARCAccelModes.mat` | Mode shapes at BARC accelerometer locations |
| `FlightStressModes.csv` | Stress modes from Abaqus for the flight configuration |
| `FlightStressModes.mat` | `FlightStressModes.csv` sorted into `.mat` format |
| `LabStressModes.mat` | Stress modes for the lab configuration |

---

### Results
Simulation results corresponding to the analyses presented in the paper.

| File | Description |
|------|-------------|
| `Results_SLAP_Buzz.mat` | Results from applying SLAP-Buzz |
| `Results_SLAP_Buzz_Elastic.mat` | Results from SLAP-Buzz excluding the modal acceleration metric |
| `Results_SLAP_Control.mat` | Results from applying SLAP-Control |
| `Results_SLAP_Control_4Shakers.mat` | Results from applying SLAP-Control with 4 shakers |
| `Results_SLAP_Control_Elastic.mat` | Results from SLAP-Control (6 shakers) excluding the modal acceleration metric |

---

### Scripts
MATLAB scripts for generating FRFs, simulating tests, and processing mode shapes.

| File | Description |
|------|-------------|
| `SLAPScript.m` | **Main script** — simulates the tests whose results appear in the paper |
| `SimulateTestScript.m` | Simulates a MIMO test for the initial environment |
| `FlightEnvironmentScript.m` | Generates the flight environment |
| `FlightFRFScript.m` | Builds the FRF relating flight forces to acceleration |
| `LabFRFScript.m` | Builds the lab control FRF |
| `StressFRFScript.m` | Builds FRFs relating shaker/lab forces to stress at locations of interest |
| `GetAbqModesScript.m` | Loads modes from Abaqus `.inp` files into `.mat` format |
| `StressModesScript.m` | Converts `.csv` stress mode outputs from Abaqus to `.mat` files |
| `GetLabLocationsScript.m` | Selects accelerometer and shaker locations for the lab test setup |
| `GetRocketLocationsScript.m` | Same as above, but for the flight configuration |
| `GetModesAtAccels.m` | Retains only the modes at BARC accelerometer locations (for modal filtering) |
| `ModalFilterScript.m` | Checks how well modal filtering performs |