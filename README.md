#  Surrogates for Information Dynamics (SID) Toolbox

MATLAB toolbox with implementations of surrogate techniques for assessing the presence of self-dependencies, nonlinearities, and coupling in univariate and bivariate systems.

The parameters for the main estimation functions surr_ISknn and surr_MIRknn are the following:
* Y - the analyzed multivariate process, organized as a matrix with rows as samples and columns as processes,
* V - the assigned embedding vector,
* jj - index of the column for one investigated process,
* ii - index of the column for the second investigated process (MIR only),
* k - number of neighbors for the KNN estimation,
* metric - distance metric for the KNN estimation, and
* surr - parameter determining which type of surrogate to be created before estimating the measure.

To enable the creation of the proposed surrogates, the functionality is controlled with the surr parameter, where for _surr_ = 1 the function applies the shuffling method on defined columns of the observation matrix, for _surr_ = 2 applies the IAAFT method, and for _surr_ = 0 it does not apply any surrogate method and the measure is estimated from the original data.

In order to demonstrate the functionality of the proposed methods, 2 demo scripts are placed within the root directory of the toolbox, /SID, along with a single sample data obtained from (Faes et al., 2011) concerning the RR and respiration series.

In both demo IS estimation and demo MIR estimation, we demonstrate the steps to evaluate the significance of the measure and the presence of nonlinearities within the RR series and respiration series. The script goes through a basic collection of steps which can be condensed as:
* loading and applying z-score normalization on the data,
* defining the parameters needed for the measure estimations,
* applying the measure estimation function with the surr parameter 0,
* applying the measure estimation function with the surr parameter 1 and 2, repeated num surrogates times,
* obtaining the percentile values and comparing them to the originally estimated measure.

Along with the demo scripts, all the codes used in simulations to test the proposed surrogate framework (Pinto et al., 2024) are presented in the /Simulations folder.

**References:**

* Faes, L., Nollo, G., and Porta, A. (2011). Information domain approach to the investigation of cardiovascular, cardio-pulmonary, and vasculo-pulmonary causal couplings. Frontiers in Physiology 2 NOV. doi:10.3389/fphys.2011.00080

* Pinto, H., Lazic, I., Antonacci, Y., Pernice, R., Gu, D., Barà, C., Faes, L., Rocha, A.P. (2024) Testing Dynamic Correlations and Nonlinearity in Bivariate Time Series through Information Measures and Surrogate Data Analysis. _Submitted to Frontiers in Network Physiology_.
