#  Surrogates for Information Dynamics (SID) Toolbox

MATLAB toolbox with implementations of surrogate techniques for assessing the presence of self-dependencies, nonlinearities, and coupling in univariate and bivariate systems.

In order to demonstrate the functionality of the proposed methods, 2 demo scripts are placed within the root directory of the toolbox, /SID, along with a single sample data obtained from (Faes et al., 2011) concerning the RR and respiration series.

In both demo IS estimation and demo MIR estimation, we demonstrate the steps to evaluate the significance of the measure and the presence of nonlinearities within the RR series and respiration series. The script goes through a basic collection of steps which can be condensed as:
• loading and applying z-score normalization on the data,
• defining the parameters needed for the measure estimations,
• applying the measure estimation function with the surr parameter 0,
• applying the measure estimation function with the surr parameter 1 and 2, repeated num surrogates times,
• obtaining the percentile values and comparing them to the originally estimated measure.
