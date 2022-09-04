# D.K. Gaskins' Portfolio

## Codes
Note: My coding eniviroment is (MATLAB R2022a with the [curve_fitting_toolbox](https://www.mathworks.com/products/curvefitting.html) and the statistics_toolbox).
To run my code, you can download it from GitHub by forking this repository ([How to Fork a Repository](https://docs.github.com/en/get-started/quickstart/fork-a-repo))
and ensuring you have the necessary toolboxes from MATLAB. I am always happy to receive any Bug Findings or Feature Requests, though I'm not in any way expecting them. 

### Main Scripts
Located in deloragaskins/pharm/
1. noncompartmental_analysis
2. two_comp
3. two_dosing_regimes
4. pk_plus_hill_pd

### Supporting Functions and Templates
Located in deloragaskins/pharm/Toolbox

* A frequent objective in the main scripts is to obtain a concentraction vs. time profile for a given model and set of parameters. 
  In each case, a function is needed to represent the differential equations that will be integrated in the model. These types of functions
  have the string "DERIVS" as their prefix which is then further specified by model type and followed by a suffix which intdicates their parameterization (currently
  only "rate_constants" although this could be easlity expanded to "clearances"). In some cases, the derivatives themselves call a helper function
  because this expression is used in multiple models (so far there is only one : Hill_Effect.m). Additionally, when multiple does of drug are adminstered during the
  time course of interest there is a need to handle the integration in a piece-wise fashion. The files that do this have the string "run_dosing_course" as their prefix   which is then further specified by model type. 
  
* The main scripts often utilize plots to be used at runtime by the user (e.g. to deterimine constants) or for data display purposes. The code   
  "PLOT_sample_figure_codeblock." can be used as a reference when building customized figures. 

* The following list indicates the main files and the files from the toolbox which support them. 
1. noncompartmental_analysis.m 
   1.DERIVS_1comp_LinearAbs_Linear_Elim__rate_constants
2.  two_dosing_regimes.m
   1. run_dosing_course__two_compartment.m
   2. DERIVS_2comp_LinearAbs_Linear_Elim__rate_constants.m 
3. pk_plus_hill_pd.m
   1. run_dosing_course__two_compartment_PD.m
   2. (DERIVS_2comp_LinearAbs_Linear_Elim__rate_constants_hill.m)
   3. Hill_Effect.m
4. PLOT_sample_figure_codeblock

## Written Texts
1. PH.D. THESIS_gaskins
2. QUAL_gaskins
