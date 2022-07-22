Fit normative data (tutorial)
=============================

The following tutorial illustrates how to use the ScanOMetrics package to fit a polynomial model to a normative dataset.

Make sure ScanOMetrics package has been properly [installed](/usage#Installation).

Import the ScanOMetrics_project class::

    from scanometrics.core import ScanOMetrics_project

Instanciate the class with the path to your BIDS dataset folder::

    SOM = ScanOMetrics_project('/home/david/Documents/DATA/test_precomputed_dataset')

Subjects can be loaded using the `load_subjects()` method. The `cov2float` argument is a generic function
to convert variables in `participants.tsv` to float. In this case, we use a lambda function to convert
males 'M' to 0 and females 'F' to 1. The `load_subjects()` function also allows to set a list of covariates
to include, with the `covariates_include` parameter::

    SOM.load_subjects(cov2float={'sex': lambda s: 0 if s=='M' else 1 if s=='F' else s},
                      covariates_include=['age', 'sex'])

The previous piece of code loads participant information from the `participants.tsv` file in your
BIDS directory. Assuming subjects were processed according to a given pipeline, derived metrics can
be loaded using `load_proc_metrics()` method. By default, ScanOMetrics assumes subjects were processed
using the Freesurfer pipeline, which is used to gather all information into a NxM numpy matrix::

    SOM.load_proc_metrics()

Outliers based on standard deviation from samples with matching age can be flagged using `flag_outliers()`::

    SOM.flag_outliers(1.5)

Finally, the normative model (in this case Leave-One-Out Cross-validation of a polynomial model) can be
set with `set_normative_model()` and fitted with `normativeModel.fit()`::

    SOM.set_normative_model('LlocvPolynomial')
    # Sets uncertainty to one as computation of uncertainty not implemented yet.
    import numpy as np
    SOM.normativeModel.fit(SOM.measured_metrics, SOM.metric_names, np.ones((1,SOM.measured_metrics.shape[1])), SOM.outliers, SOM.covariate_values, SOM.covariate_names, flag_opt=1, N_cycl=1, width=0)


