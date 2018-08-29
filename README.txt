Author 	Luck Peerlings
Date 	2018-08-22

This repository is a direct copy of commit commit f1dfe83e96f84150a2839e87de2d631d723a0fe4 of the ImpedanceEduction branch. Code can also be found on https://github.com/LuckPeerlings/nPortAnalysis

This repository consists of the following classes

Main Classes
------------
@AreaExpansionModel
	Class to calculate model results from the various simplified models
@MonteCarlo
	Class to calculate the uncertainty using the MonteCarlo technique. It has to be used in conjuction with the UncertainVariable class.
@MultiVariateAnalysis
	Class to calculate the uncertainty on class properties using a linear multivariate analysis. Has to be used in conjuction with the UncertainVariable Class
@NPortAnalysis
	Class to evaluate scattering matrix methods. Various wavenumber models, decomposition techniques and optimization routines are implemented. 
	The possibility to analyze higher order mode measurements is also present. Only the code for the rectangular case has been tested.
@PronyMethod
	Class to determine the wavenumbers from non-equispaced microphone distance using different prony methods
@UncertainVariable
	Class to store information such as the covariance matrix and correlations with other variables. Many functions have been overloaded such that operations based on normal (float, integer) variables also function with variables of the class Uncertain Variable. This class has to be used in conjuction with the MultiVariateAnalysis and the MonteCarlo method to be able to perform the uncertainty analysis
@ZeroNumerical
	Class to calculate the zeros of an analytical function by performing contour integrals in the complex domain and using the argument principle.
	
Main Classes: Not finished
------------
@ImpedanceEduction
	Class to educe the impedance using the mode matching technique. This class is NOT finished	
	
Sub Classes
----------	
@UncertaintyAnalysisNPort
	Subclass of the MultiVariateAnalysis, which allows for the calculation and plotting of uncertainty intervals of information based in an object of NPortAnalysis. At least one input property in the NPortAnalysis object has to be of the class UncertainVariable.

@UncertaintyAnalysisProny
	Subclass of the MultiVariateAnalysis, which allows for the calculation and plotting of uncertainty intervals of information based in an object of PronyImpedance. At least one input property in the PronyImpedance onbject has to be of the class UncertainVariable.
@PronyImpedance
	Subclass of the PronyMethod where the determination of the impedance of liners using the prony method is implemented.	
@MonteCarloProny
	Subclass of the @PronyImpedance to incorporate the @MonteCarlo class	
	
Functions
----------
\Methods
	Functions that are used by the classes
	
	- GasProperties.m, AirProperties.m, GasProperties_chk.m, AirProperties_chk.m
		Functions to calculate the properties of air and other types of gasses, with their respective functions to check if the inputs are correct. These checks have been implemented seperately to avoid doing the checks continiously when for example performing a MonteCarlo simulation
	- MergeStructs.m
		Function to merge different structs, extensively used in the MultiVariate Analysis
	- MicPosCalibration.m, MicPosCalibration_Residual.m 
		Functions to determine the microphone positions when measuring the rigid wall. The input data is based on the VXI-measurement script. For MicPosCalibration_Residual.m the residual of the system of equations is minimized by optimizing the temperature, which may increases the accuracy of the results.

\PhasePlot
	Functions to plot complex variables in the complex domain using phase plots. Code is copied from the Mathworks repository. It is part of the the book "Visual Complex Functions" by Elias Wegert

\TestScripts
	Testcripts have been written for the above classes to show how the classes should be used. The necessary data is present in the folders.
	
\Tools
	- DetermineComplexPressures.m 
		Function to determine complex amplitudes using synchronous demodulation.
	- ComparisonCalibration.m 
		Function to compare different relative calibration measurements for the microphones. The data should be measured with the VXI_Calibration.m program in the VXI_MeasurementFunctions folder
	- fnc_plot_gaussian_ellipsoid.m 
		Function from mathworks to be able to plot uncertainty ellipses in the complex domain.
	- uipickfiles.m
		Function from Mathworks, it is an improved uigetdir/uigetfile
	- WriteToTextFile_PGF.m
		Function to easily write matlab data to a tab delimited file with headers. To use with the LaTeX package PGF plots.

\VXI_MeasurementFunctions
	
	Main Functions
	--------------
	- VXI_Calibration.m
		Function to perform the microphone calibration. Run the function with a computer attached to a VXI function. Follow the instructions and questions in the command window (case sensitive) to start the calibration.
	
	- VXI_nPort_ConfigurationFile.m
		Function to create a configurationfile that can be used in conjuction with the VXI_TwoPort_Measurement.m 
		It is important to double check that the address is correct for the VXI system such that all the modules are initialized.
		Information for each settable parameter is in the comments of the configurationfile.
	
	- VXI_nPort_Measurement.m
		Function to perform the nPortMeasurement. Follow the instructions to start the measurement. If no configuration file is used, the parameters in the function itself are used to perform the measurment. This is helpful to for example speed up the error searching process in the measurement or when implementing new functions.
		
	Sub Functions
	-------------
	- mailgun.m, curl/
		Files necessary for sending an email when the measurement is finished.
	- vt1432.dll, vt1432_check_status.m 
		Necessary files to communicate with the VXI system
	- ScrambleFrequency.m
		Short function to scramble the frequencies and make certain that the N sources do not excite tones that are multiple of each other.
		
	