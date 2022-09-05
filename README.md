# D.K. Gaskins' Portfolio

## Codes
Note: My coding environment is (MATLAB R2022a with the [curve_fitting_toolbox](https://www.mathworks.com/products/curvefitting.html) and the statistics_toolbox). To run my code, you can download it from GitHub by forking this repository ([How to Fork a Repository](https://docs.github.com/en/get-started/quickstart/fork-a-repo)) and ensuring you have the necessary toolboxes from MATLAB. I am always happy to receive any Bug Findings or Feature Requests, though I'm not in any way expecting them. 

### Main Scripts
Located in deloragaskins/pharm/

These scripts call functions from the Toolbox folder. A list of the dependencies for each script can be found in the "Supporting Functions and Templates Section". 

1. Noncompartmental Analysis is a simple method for fitting key kinetic parameters and for evaluating the exposure of a drug. I used a one compartmental model to generate a concentration vs time profile for IV and PO dose administrations and then added noise. I then fit these curves to fit the following values: 
    * k_elim, c_0,  IV_AUC,  Vd, t_half, CL
    * PO_AUC, C_max, t_max, F

   (code name: _noncompartmental_analysis.m_)
 
2. Simbiology is a commonly used program in PK/PD and QSP circles. To obtain familiarity with the model builder and analyzer, I took data from 'Concepts in Clinical Pharmacokinetics'( ISBN: 978-1-58528-591-4, Page 89 Problem 6-9) and used a two variable model from the PK library to fit the clearances and the compartment volumes. 

   (code name: _two_comp.sbproj_)

3. Administering a loading dose which is higher than the remainder of the doses administered allows the drug concentration to reach steady state sooner. I designated two different dosing regimens, one with a loading dose and a maintenance dose and the other using a constant dose value. The code relies on a supporting function _dosing_choser.m_ which I used to decide the maximum dose, the maintenance dose and the dosing frequency. I obtained the concentration profiles for comparison. 

    (code name: _two_dosing_regimes.m_)
    
4. PK models can be coupled to PD models to determine the physiological/biochemical effect within the body that occurs given a particular dosing regimen. I took a two compartmental PK model and coupled it with a Hill equation to show a decrease in tumor size (n) with a drug treatment. 
   
   (code name: _pk_plus_hill_pd.m_)

### Supporting Functions and Templates
Located in deloragaskins/pharm/Toolbox

* [GENERATING CONCENTRATION VS TIME DATA] A frequent objective in the main scripts is to obtain a concentration vs. time profile for a given model and set of parameters. In each case, a function is needed to represent the differential equations that will be integrated in the model. These types of functions have the string "DERIVS'' as their prefix which is then further specified by model type and followed by a suffix which indicates their parameterization (currently only "rate_constants" although this could be easily expanded to "clearances"). The parameters for the ODEs are passed to these functions as attributes of an object that is an argument specified for the chosen ODE solver. In some cases, the derivatives themselves call a helper function because this expression is used in multiple models (so far there is only one : "Hill_Effect.m"). Additionally, when multiple doses of drug are administered during the time course of interest there is a need to handle the integration in a piece-wise fashion. The files that do this have the string "run_dosing_course" as their prefix which is then further specified by model type. 

* [DETERMINING MAX DOSE and viable DOSE and DOSING INTERVAL] PK can be used to predict viable dosing amount and dosing schedule. The function "dosing_chooser" takes the minimum effective concentration and the maximum tolerated concentration and allows the user to estimate the maximum dose which keeps the concentration under the maximum tolerated concentration. The used can then vary the dose and dosing interval to keep the concentration within the therapeutic window. Once these values are determined, they can be fixed so the function returns them when called.    
  
* [DATA VISUALIZATION] The main scripts often utilize plots to be used at runtime by the user (e.g. to determine constants) or for data display purposes. The code "PLOT_sample_figure_codeblock.m'' can be used as a reference when building customized figures. 

The following list indicates the main files and the functions from the toolbox which support them. 

1. noncompartmental_analysis.m 
   1. DERIVS_1comp_LinearAbs_Linear_Elim__rate_constants.m
2. two_dosing_regimes.m
   1. dosing_chooser.m
   2. run_dosing_course__two_compartment.m
   3. DERIVS_2comp_LinearAbs_Linear_Elim__rate_constants.m 
3. pk_plus_hill_pd.m
   1. run_dosing_course__two_compartment_PD.m
   2. DERIVS_2comp_LinearAbs_Linear_Elim__rate_constants_hill.m
   3. Hill_Effect.m
4. no main file: PLOT_sample_figure_codeblock.m

## THESIS and QUALIFYING EXAM
I uploaded my Ph.D. thesis (PH.D. THESIS_gaskins.pdf) and my qualifying exam (QUAL_gaskins.pdf) here for your perusal. 


