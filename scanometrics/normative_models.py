"""
Module for implementation of developmental models
"""
import abc
import os
import numpy as np
from glob import glob1

import abc
from scanometrics.utils import logging
from scipy.stats.mstats import mquantiles
from scipy.stats import norm
from scipy.stats import f as fstat
from scipy.stats import chi2


class normative_model_template():

    def __init__(self, training_dataset_id, models_folder=None):
        """Constructor for the developmental model. Models are identified by unique combinations of model_name and
        model_dataset used to train it. The model folder is the folder to save/load fitted model parameters.

        :param training_dataset_id: unique ID of training dataset
        :type training_dataset_id: string
        :param models_folder: (optional) folder to save/load model parameters. Defaults to ScanOMetrics install folder.
        :param models_folder: string
        """
        self.training_dataset_id = training_dataset_id  # Code for dataset used for training. Datatype: string
        model_ID = self.model_name + '_' + training_dataset_id  # Model ID suposed to be unique for keeping track, also used to generate model folder in scanometrics installation folder
        self.model_id = model_ID
        if models_folder is None:  # Folder (in scanometrics installation folder by default) used to store models fitted parameters and subsequent loading
            models_folder = os.path.join(os.path.basename(__file__), 'resources', 'normative_models')
        self.model_folder = os.path.join(models_folder, model_ID)
        if not os.path.exists(self.model_folder):
            os.makedirs(self.model_folder)
        # if subjects is None:  # Added at some point but probably not required: users should be able to load a trained model without downloading the training data.
        #     subjects = glob1(os.path.join(bids_folder, 'derivatives', 'scanometrics', 'sub-*'))
        # self.subjects = subjects
        self.parameter_fit = None  # estimated model parameters after fitting on self.X (should be saved/loaded in/from model_folder/parameter_fit.npy)
        self.norm_parameter_fit = None  # estimated model parameters after fitting on self.norm_X (should be saved/loaded in/from model_folder/norm_parameter_fit.npy)
        self.parameter_stat = None  # statistics of parameters trained on self.X (should be saved/loaded in/from model_folder/parameter_stat.npy)
        self.norm_parameter_stat = None  # statistics of parameters trained on self.norm_X (should be saved/loaded in/from model_folder/norm_parameter_stat.npy)
        self.X = None  # measured_metrics used to fit the model (should be small enough to be worth saving, specially as can be used to train several models). Should be loaded from model_folder/<proc_pipeline_ID>_training_X.npy
        self.proc_pipeline_id = None  # ID of pipeline used for generate measured_metrics used to fit the model

    def load_model_parameters(self):
        """
        Loads model parameters.
        """
        if not os.path.exists(os.path.join(self.model_folder, 'parameter_fit.npy')):
            logging.ERROR("File 'parameter_fit.npy' not found in %s. Run 'scanometric_project.normativeModel.fit()' and 'save()' before loading." % (self.model_folder))
        self.parameter_fit = np.load(os.path.join(self.model_folder, 'parameter_fit.npy'))
        self.norm_parameter_fit = np.load(os.path.join(self.model_folder, 'norm_parameter_fit.npy'))
        self.parameter_stat = np.load(os.path.join(self.model_folder, 'parameter_stat.npy'))
        self.norm_parameter_stat = np.load(os.path.join(self.model_folder, 'norm_parameter_stat.npy'))

    def save_model_parameters(self):
        """
        Saves model parameters.
        :return:
        """
        np.save(os.path.join(self.model_folder, 'parameter_fit.npy', self.parameter_fit))
        np.save(os.path.join(self.model_folder, 'norm_parameter_fit.npy', self.norm_parameter_fit))
        np.save(os.path.join(self.model_folder, 'parameter_stat.npy', self.parameter_stat))
        np.save(os.path.join(self.model_folder, 'norm_parameter_stat.npy', self.norm_parameter_stat))

    def load_X(self):
        """
        Load data for analysis and comparison with normative model, or training of normative model.
        """
        file_path = os.path.join(self.model_folder, self.proc_pipeline_id + '_training_X.npy')
        if not os.path.exists(file_path):
            logging.ERROR("Training data file %s not found in %s. Make sure to process training data, load it and save it before reloading." % (self.proc_pipeline_id + '_training_X.npy', self.model_folder))
        self.X = np.load(file_path)

    def save_X(self):
        """
        Save measured_metrics matrix.
        """
        file_path = os.path.join(self.model_folder, self.proc_pipeline_id + '_training_X.npy')
        np.save(file_path, self.X)

    """def fit(self):
        # Should return 'estimated_params' and 'residuals' array for developmental tarjectory fitting."""

class LlocvPolynomial(normative_model_template):
    """
    Normative model based on Leave-One Out Cross-Validation and polynomial fit. 'training_dataset_id' mixes normative
    model name with name of training dataset to keep track of combination used for training.
    """


    def __init__(self, training_dataset_id, models_folder):
        """
        Class constructor. Initialize everythin with empty arrays, matrices and lists.
        :param training_dataset_id: unique ID of dataset being processed.
        :type training_dataset_id: string
        """
        self.model_name = 'LlocvPolynomial'
        super().__init__(training_dataset_id, models_folder)
        self.fit_deg         = None
        self.fit_poly        = None
        self.fit_good        = None
        self.fit_dev         = None
        self.est_dev         = None
        self.fit_sml         = None
        self.fit_lrg         = None
        self.residual        = None
        self.X               = None
        self.metric_names    = None
        self.covariates      = None
        self.covariate_names = None

    def load_model_parameters(self):
        super().load_model_parameters()
        # Some model specific parameters

    def save_model_parameters(self):
        super().save_model_parameters()
        # Some model specific parameters

    def load_X(self):
        super().load_X()
        # Some model specific operation, although can probably be a general function without overriding...

    def save_X(self):
        super().save_X()
        # Some model specific operation, although can probably be a general function without overriding...

    def fit(self, X, metric_names, uncertainty, outliers, covariates, covariate_names, flag_opt, deg_max=None, frac=20, alpha=0.01, N_cycl=10, width=0.5):
        """Computes 'estimated_parameters' and 'residual' matrices for the model, according to 'measured_metrics' array.
        X is the training dataset, given as a numpy array (a copy is created in self.X to be saved/loaded). Outliers in
        the X matrix still have values for computation of residuals.
        TODO: set uncertainty computation"""
        from numpy.polynomial import Polynomial as poly_model

        self.X = X.copy()
        self.metric_names = metric_names
        if 'age' not in covariate_names:
            logging.ERROR("Covariates in normative model requires at least an 'age' variable.")
        else:
            self.covariates = covariates.copy()
            self.covariate_names = covariate_names
        if flag_opt and (deg_max is None):
            deg_max = int(np.floor(len(covariates[:, covariate_names.index('age')])/frac))
        # Convert ages to integer range array
        # Not sure what's the advantage of it compared to continuous/non-uniform values from dataset...
        # Ask Christian about it

        # Supposed to sample uniformly wrt age, but width of zero means all subjects are selected...
        # Setting uniform_selection to ones instead of calling uniform_subsample
        if width == 0:
            uniform_selection = np.ones((1, self.X.shape[0]), dtype='bool')
        else:
            uniform_selection = uniform_subsample_scheme(covariates[:, covariate_names.index('age')], N_cycl, width)

        # Initialize empty arrays
        age      = covariates[:, covariate_names.index('age')].copy()
        age_vec = np.arange(np.floor(age.min()), np.ceil(age.max()) + 1)
        deg_opt  = np.full((self.X.shape[1], N_cycl), np.nan, dtype='int')
        polynom  = np.full((deg_max+1, self.X.shape[1], N_cycl), np.nan)
        fit_good = np.zeros((self.X.shape[1], N_cycl))
        residues = np.full(self.X.shape, np.nan)
        fit_ave  = np.full((len(age_vec), self.X.shape[1]), np.nan)
        fit_dev  = np.full((len(age_vec), self.X.shape[1]), np.nan)
        est_dev  = np.full((len(age_vec), self.X.shape[1]), np.nan)
        sml      = np.full((len(age_vec), self.X.shape[1]), np.nan)
        lrg      = np.full((len(age_vec), self.X.shape[1]), np.nan)

        # Loop through the metrics
        for i in range(self.X.shape[1]):
            polyval_vec = np.zeros((len(age_vec), N_cycl))
            predict_vec = np.zeros((self.X.shape[0], N_cycl))
            fit_selection = (uniform_selection) & (~(np.isnan(self.X[:, i]))) & (~(outliers[:, i]))
            for n in range(N_cycl):
                idx_valid = np.argwhere(~outliers[:, i] & ~np.isnan(self.X[:, i]) & uniform_selection[:, n])[:, 0]  # 0 indexing to get (N,) array TODO: consider changing to array of True/False
                if len(idx_valid) > 2:
                    deg_opt[i, n] = deg_max
                    if flag_opt:
                        SSE = np.zeros(deg_max+1)  # to go from 1 to deg_max+1
                        for d in range(deg_max+1):  # goes from 1 to deg_max+1
                            c = poly_model.fit(age[idx_valid], self.X[idx_valid, i], deg=d)
                            res = self.X[idx_valid, i] - c(age[idx_valid])
                            SSE[d] = res.dot(res.T)
                        for d in range(deg_max):
                            df = len(idx_valid)-(d+1)
                            F = (SSE[d]-SSE[d+1]) / SSE[d+1] * df
                            pF = 1-fstat.cdf(F, 1, df)
                            if pF > alpha/(d+1):
                                deg_opt[i, n] = d
                                break
                    c = poly_model.fit(age[idx_valid], self.X[idx_valid, i], deg=deg_opt[i, n])
                    polynom[:deg_opt[i, n]+1, i, n] = c.coef.copy()
                    polyval_vec[:, n] = c(age_vec)
                    predict_vec[:, n] = c(age)

                    # Chi-squared statistics to test whether the fit is good at all
                    # i.e. residual variance should not be significantly larger than variance of measurement uncertainty
                    T_fit = (len(idx_valid)-1) * np.var(self.X[idx_valid, i]-predict_vec[idx_valid, n]) / uncertainty[0, i]  # <- double check how uncertainty is computed
                    p_fit = 1-chi2.cdf(T_fit, len(idx_valid)-1)
                    if p_fit > alpha:
                        fit_good[i, n] = True
            idx_cycl = np.argwhere(fit_good[i, :])
            meas_is_good = ~np.isnan(self.X[:, i]) & ~outliers[:, i]
            if len(idx_cycl) > 0:
                fit_ave[:, i] = polyval_vec[:, idx_cycl].mean(1)[:, 0]
                fit_dev[:, i] = polyval_vec[:, idx_cycl].std(1, ddof=1)[:, 0]
                for s in range(len(age_vec)):
                    idx_age = np.argwhere((age>=0.9*age_vec[s]) & (age<=1.1*age_vec[s]) & meas_is_good)[:, 0]
                    est_dev[s, i] = np.std(self.X[idx_age, i]-predict_vec[idx_age, idx_cycl].mean(1), ddof=1)
            else:
                fit_ave[:, i] = np.nanmean(polyval_vec, 1)
                if polyval_vec.shape[1] > 1:
                    fit_dev[:, i] = np.nanstd(polyval_vec, 1, ddof=1)
                else:
                    fit_dev[:, i] = np.zeros(polyval_vec.shape[0])
                for s in range(len(age_vec)):
                    idx_age = np.argwhere((age >= 0.9 * age_vec[s]) & (age <= 1.1 * age_vec[s]) & meas_is_good)[:, 0]
                    est_dev[s, i] = np.std(self.X[idx_age, i]-predict_vec[idx_age, :].mean(1), ddof=1)
            sml[:, i] = fit_ave[:, i]-1.96*np.sqrt(fit_dev[:, i]**2+est_dev[:, i]**2)
            lrg[:, i] = fit_ave[:, i]+1.96*np.sqrt(fit_dev[:, i]**2+est_dev[:, i]**2)
            if N_cycl > 0:  # original checks for size of predict_vec, should be equivalent
                residues[:, i] = self.X[:, i] - predict_vec.mean(1)


def uniform_subsample_scheme(samples, n_cyl, width, seed=None):
    """Generates n_cyl subsampling schemes with approximate uniform distribution, by selecting subsamples with probability
    inversely proportional to the density.
    Original implementation from Octave's code. Possible shorter and pythonic implementation to be tested from here:
    https://stackoverflow.com/questions/66476638/downsampling-continuous-variable-to-uniform-distribution

    :param samples: sample from which subsamples should be taken.
    :param n_cyl: number of subsamples sets to generate
    :param width: width of broadening gaussian
    :param seed: seed for np.random.seed() for testing
    :return idx: np.array of size (len(sample),n_cyl), with 0 for excluded samples and 1 for included samples
    """
    if width > 0:
        density = np.zeros(len(samples))
        for s1 in range(len(samples)):
            for s2 in range(len(samples)):
                density[s1] += np.exp(-(samples[s2] - samples[s1]) ** 2 / (2 * width**2)) / (np.sqrt(2 * np.pi) * width)
        prob_sel = 1./density  # sum is approx equal to mean of sample ?
        # prob_sel /= prob_sel.sum()
        idx = (np.tile(prob_sel[:, None], (1, n_cyl)) > np.random.rand(len(samples), n_cyl)).astype('bool')
        # idx = np.random.choice(samples, len(samples), replace=True, p=prob_sel)
    else:
        idx = np.ones((len(samples), n_cyl)).astype('bool')
    return idx
