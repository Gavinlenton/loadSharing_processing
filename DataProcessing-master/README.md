# loadSharing_processing
The main processing script and functions for processing the load sharing data

First run LS_Pipeline to extract the outputs for use in OpenSIM (e.g., IK/ID).
	This calls the acquisitionInterface and c3d2mat from MOtoNMS and 

Run LinScale.m to scale the generic model bodies 
	Includes staticElaboration function from MOtoNMS to process the static trial, a requirement for LinScale
	Fixes the talus position
	
Run Full_CAST.m 
	Places markers on the model after solving one frame of IK
	Optimises muscle parameters
	
Run BOPS_LS.m to process IK and ID for dynamic trials.
	Runs IK first
	Takes IK results and the model spline values to generate abd/add and int/ext 		rotation values and puts them into the .mot file from IK
	Open the knee DOFS of the model - convert splines to linear functions
	Run ID with the modified model and prescribed motion
	Run Point kinematics to get the COM data from he pelvis


