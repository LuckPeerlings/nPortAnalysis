Author 	Luck Peerlings
Date 	2018-08-22

To perform a (nPort) measurement, the microphones have to be calibrated first using the function situated at "VXI_MeasurementFunctions\VXI_Calibration.m"
Thereafter a configuration file has to be created using the "VXI_MeasurementFunctions\VXI_ConfigurationFile.m"
The measurement is started with the function "VXI_MeasurementFunctions\VXI_Measurement.m" and the raw time data andn auxilliary information is stored in the specified data path.
With help of the function "Tools\DetermineComplexPressures.m" the time data is postprocessed to obtain the complex pressures (transfer function) as function of frequency and source. The post processing time is significantly increased when determining the covariance matrix of the measurement results.
