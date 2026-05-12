# We are in currently updating the FEM of the rocket to fix some issues, and also greatly expanding the documentation. Watch for an update by the end of May 2026
# SLAP Rocket Example

This repository contains files and code for a rocket finite element model (FEM) case study on dynamic environment testing, first presented in *"Damage Metric-Based Vibration Testing of a Rocket Component"* by Behling, Allen, Bahr, Mayes, & DeLima. The case study applies the **Scaled Lab PSD (SLAP)** method to the removable component of the BARC system. This component was selected for its common usage and because it has multiple elastic modes below 2000 Hz, making it a rich test case.

> **Note:** A few files were too large for GitHub and are hosted externally. Download them at: https://byu.box.com/s/vbepspsz936w6v2sfjhw7x4hxon92tow. After the data has been downloaded, unzip the file in the "LargeFiles" folder. There is a text document inside the folder that explains how the file path should be organized.

---

## Getting Started

The most important file in this repository is **`Scripts/SLAPScript.m`**, which was used to simulate the tests whose results appear in the paper. Start there to reproduce the paper's findings.

If you would like to try SLAP on a different system, the workflow detailing how to generate the *.mat files used in the paper is given below.

---

## Workflow

The general workflow is given in the image below. The user should refer back to this often to know what step in the process they are on. The documentation following the first flow chart breaks down specific parts of this flow chart where needed.

![Workflow Diagram](Assets/Revised%20SLAP%20Workflow.svg)

### Running Abaqus
![Abaqus Explanation](Assets/Run%20Abaqus%20Explanation.svg)

To start, you will want to run the input decks that are located in "AbaqusFiles". There are three input decks that are needed for this case study:

- Full_Rocket.inp
    - Contains the geometry for the rocket, baseplate, and DUT.
    - Calculates all modes below 3 kHz
    - Records displacement within the *.odb file, as well as stress at specified elements on the DUT
- BARC_Baseplate.inp
    - Contains the geometry for the baseplate and DUT
    - Calculates all modes below 3 kHz
    - Records displacement within the *.odb file, as well as stress at specified elements on the DUT
- BARC_FixedBase.inp
    - Contains the geometry for the DUT
    - Calculates all Fixed Base Modes below 3 kHz
    - Records displacement within the *.odb file, as well as stress at specified elements on the DUT

### Extracting Mode Shapes
![Mode Shape Explanation](Assets/Mode%20Shapes%20Explanation.svg)

After the jobs have finished running, open up a command prompt and navigate to "AbaqusFiles" folder. This is where the .odb files should be saved when the analyses are run. If they are saved elsewhere, move them to the "AbaqusFiles" folder before continuing. In the command prompt, type "abaqus python odb_to_matlab.py --odb <Name_Of_Your_Model>.odb". This will run the python script that tranlates the displacement data in the odb file for each mode into a format that MatLab can access and read. You will need to do this for all of the models that you run. (**Disclaimer: If you change the names of the input files above you will need to change all of the file names in the load commands for the MatLab scripts that are used.**)

### Stress
![Stress Explanation](Assets/Stress%20Explanation.svg)

After the odb file is generated, open up the odb in Abaqus CAE and go to the Visualization Module. Then on the toolbar go to Report -> Field Output... This will bring up a dialog box that will ask how and what you would like exported. To perform SLAP, we need the stress modes. For this, we need to export $\sigma_{xx}$, $\sigma_{yy}$, $\sigma_{zz}$, $\tau_{xy}$, $\tau_{xz}$, and $\tau_{yz}$. However, we don't want to calculate these values for every element in the model. What we did for this is create an element set that we thought was representative enough of the full stress distribution (see Section 3.4 of the paper). In the Abaqus dialogue box that we opened, we will click "All active steps/frames" at the top, and select S11, S22, S33, S12, S13, S23 in the variable tab. Then go to the setup tab, select the csv format option, name the file, and hit apply to export the values to a csv. The next step in the workflow (running "StressModesScript.m"), converts the *.csv files to *.mat files.

### Location Scripts
![Location Explanation](Assets/Location%20Script%20Explanation.svg)

The Locations scripts are used to determine where you want your force inputs and accelerometers on both the rocket and the lab setup. For these scripts, two interactive figures are generated with a grid of points that can be selected. The script does not automatically save the points that you select. We recommend selecting all of the points that you would like for your test on both figures, and then running the two save commands at the bottom of the script before the function definitions. It is also uncertain whether order matters when picking points. For example, if you chose 10 accelerometer locations for flight, you should choose the same 10 locations, in the same order, on the test stand.

### Setup Suite
![Setup Explanation](Assets/Setup%20Suite%20Explanation.svg)

Before you can run SLAP, a variety of FRFs need to be generated. The scripts needed to generate the FRFs and flight environment are what we call the setup suite. There are five scripts that need to be run in order to have all the data needed to perform SLAP. The FRF scripts can be run in any order, but they must all be ran before the FlightEnvironmentScript.m and GetModesAtAccels.m scripts.

### Modal Filter Check
After the Setup Suite has been run, it is recommended that the user runs ModalFilterScript.m, to make sure that the chosen accelerometer locations are enough to capture all of the modes of interest. SLAPScript.m can be run without this step, but it is a good sanity check.

### Simulate Test Script
SLAPScript.m runs SLAP on 100 randomized environments to know if it would work reliably. Before running all 100 cases, it is recommended that the user run SimulateTestScript.m. This simulates one environment and displays the PSDs so that the user can get a feel for what each iteration of the SLAPScript.m looks like. Again, running this script before SLAPScript.m is not necessary, but a good check.

### SLAP 'em Silly
You're ready run SLAPScript.m

---

## Repository Structure

### AbaqusFiles
Contains `.inp` files used to extract stress and displacement modes for flight, lab, and fixed-base configurations. Also contains a python script that was used to translate infromation from the Abaqus *.odb files into *.mat files.

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