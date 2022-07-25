"""Core ScanOMetrics classes and methods"""
import multiprocessing
import string
import random
import numpy as np
import os
from glob import glob1
from scanometrics import normative_models, processing
from scanometrics.utils import logging
from csv import DictReader
from scipy.stats import norm
from scipy.stats.mstats import mquantiles


class ScanOMetrics_project:
    """
    Class defining a ScanOMetrics project.
    Subjects should be stored according to BIDS data structure for multiple sessions: https://bids.neuroimaging.io/
    Demographic data should be saved in <bids_database>/participants.tsv
    DICOM data can be converted to bids format using scanometrics.utils.dicom2nifti() function.
    """

    ###############
    # Constructor #
    ###############

    def __init__(self, bids_database, proc_pipeline='freesurfer', subjects_in=None, subjects_out=None, covariates_in=None, covariates_out=None):
        """
        ScanOMetrics_project constructor from bids_database path. BIDS database should at least contain a participants.tsv
        file to load subject names and covariates. Order should match the order in tables in stats2table folder.
        :param bids_database: path to BIDS database
        :type bids_database: string
        :param subject_list: list of subject names to be included in the project (defaults to None to use all sub-* in bids folder)
        :type subject_list: list
        :param metric_proc_pipeline: name of processing pipeline whose 'proc2metric' function should be called
        :type metric_names: string
        """
        # Set directory of bids database with processed subjects
        self.bids_database = bids_database
        # Reset list of subjects (loaded by user calling self.load_subjects())
        self.subject = {}
        # Reset metrics array and metric_names, to be filled when calling self.load_subjects()
        self.measured_metrics = np.array([], dtype='float')
        self.covariate_values = np.array([], dtype='float')
        self.metric_names = []
        self.covariate_names = []
        self.outliers = np.array([], dtype='bool')
        self.set_proc_pipeline(proc_pipeline)
        # Idea would be to refactor everything so that a function proc2metric encodes how to retrieve the data from processing step
        # into a 'measured_metrics' matrix and 'metric_names' list, on which normative models can be trained to generate
        # 'estimated_metrics' and 'metric_residuals'

        # Set and load normative model
        self.normativeModel = None
        # self.set_normative_model() <- should be called by user with appropriate model settings

    ###################
    # SUBJECT LOADING #
    ###################

    def load_subjects(self, subjects_in=None, subjects_out=None, covariates_include=None, covariates_exclude=None,
                      cov2float=None):
        """Load subject data from <self.bids_database>/participants.tsv, with inclusion/exclusion of subjects based on
        their ID, and inclusion/exclusion of covariates based on covariate name. Also keeps track of sessions and adds
        relevant information to self.subject and self.repeated_metrics. Repeated measures can be ignored by manually
        deleting rows after running load_subjects(), to allow user to specify which repeat to keep/exclude at runtime.
        Assumes repeated measures are a separate scan in <self.bids_database>/sub-<ID>/ses-<label>/anat/T1.nii.gz from
        which all metrics are extracted. Session info should be in <self.bids_database>/sub-<ID>/sessions.tsv (one row
        per session, with a mandatory column named 'session_id' and parameter(s) that change(s) between sessions).

        :param subjects_in: list of subjects to include. Defaults to None to load all subjects in participants.tsv
        :param subjects_out: list of subjects to exclude. Defaults to None to load all subjects in subjects_in
        :param covariates_include: list of covariates to include. Defaults to None to load all covariates in participants.tsv
        :param covariates_exclude: list of covariates to exclude. Defaults to None to load all covariates in covariates_include
        :param cov2float: dictionary of functions to convert covariates to float. Functions are associated to covariates
         by matching key (eg {'sex':lambda s : 0 if s=='M' else 1 if s=='F' else s})
        """
        # BS: repeated measures can be a covariate, so should not be tracked in a separate array.
        # Current idea: put everything in covariate_values and covariate_names. Loop a first time to gather all covariates
        # from participants.tsv and sessions.tsv files, remove those in exclude_covariates and those not in include_covariates
        # Then loop again and collect covariate_values in appropriate order, repeating unchanged values accross sessions
        # BS: current implementation drives all subjects to have a sessions.tsv file to save repeated variable value for
        # all subjects, even if there are subjects for which there is only one session.
        # Including sessions in participants.tsv is tempting, despite bids recommendations going against...
        self.subject = {}
        self.covariate_values = []  # Start with a list to be converted to numpy array after looping through tsv files
        # Loop through subjects to check for _sessions.tsv file with potential longitudinal changes on (new) variables
        with open(os.path.join(self.bids_database, 'participants.tsv'), 'r') as f:
            reader = DictReader(f, delimiter='\t')
            self.covariate_names = list(next(reader).keys())
            for row in reader:
                ID = row.pop('participant_id')
                ses_file = os.path.join(self.bids_database, ID, ID + '_sessions.tsv')
                if os.path.exists(ses_file):
                    with open(ses_file, 'r') as f_ses:
                        for k in next(f_ses).strip().split('\t'):
                            if k not in self.covariate_names:
                                self.covariate_names.append(k)
        if 'participant_id' in self.covariate_names:
            self.covariate_names.remove('participant_id')
        if 'session_id' in self.covariate_names:
            self.covariate_names.remove('session_id')
        # Exclude covariates in covariates_exclude or not in covariates_include lists
        for k in self.covariate_names:
            if (covariates_exclude is not None and k in covariates_exclude)\
            or (covariates_include is not None and k not in covariates_include):
                self.covariate_names.remove(k)
        # Loop again, this time filling values for each subject
        with open(os.path.join(self.bids_database, 'participants.tsv'), 'r') as f:
            reader = DictReader(f, delimiter='\t')
            for row in reader:
                ID = row.pop('participant_id')
                # Filter out subject not in subjects_in or in subjects_out lists
                if ((subjects_in is not None) and (ID not in subjects_in)) or ((subjects_out is not None) and (ID in subjects_out)):
                    continue
                # Convert to float
                if cov2float is not None:
                    for k in [k for k in list(cov2float.keys()) if k in row.keys()]:
                        row[k] = cov2float[k](row[k])
                # # Convert everything else to float (might not be needed as list converted to float array below
                # for k in list(row.keys()):
                #     row[k] = float(row[k])
                # Check if sessions.tsv file exists and add values to row
                ses_file = os.path.join(self.bids_database, ID, ID + '_sessions.tsv')
                if os.path.exists(ses_file):
                    with open(ses_file, 'r') as f_ses:
                        ses_reader = DictReader(f_ses, delimiter='\t')
                        for ses_row in ses_reader:  # Loop through repeats, update row dictionary each time
                            ses_id = ses_row.pop('session_id')
                            for k in ses_row.keys():
                                row[k] = ses_row[k]
                            # Convert to float
                            if cov2float is not None:
                                for k in [k for k in list(cov2float.keys()) if k in ses_row.keys()]:
                                    row[k] = cov2float[k](row[k])
                            # Fill missing values with nan
                            for k in [k for k in self.covariate_names if k not in row.keys()]:
                                row[k] = np.nan
                            # Fill in self values according to self.covariate_names order
                            self.covariate_values.append(np.array([row[k] for k in self.covariate_names], dtype='float')[None, :])
                            self.subject[ID+'_'+ses_id] = row
                else:
                    # Fill missing values with nan
                    for k in [k for k in self.covariate_names if k not in row.keys()]:
                        row[k] = np.nan
                    # Fill in self values according to self.covariate_names order
                    self.covariate_values.append(np.array([row[k] for k in self.covariate_names], dtype='float')[None, :])
                    # Keep track of row subject ID
                    self.subject[ID] = row
            self.covariate_values = np.vstack(self.covariate_values)


    ############################
    #   DATA (PRE)PROCESSING   #
    ############################

    def set_proc_pipeline(self, metric_proc_pipeline):
        if hasattr(processing, metric_proc_pipeline):
            self.metric_proc_pipeline = getattr(processing, metric_proc_pipeline)
        else:
            logging.ERROR("Processing pipeline %s not found in scanometrics.processing module. Make sure the module file"
                          " exists, that it has been added to the __init__.py file, and that scanometrics is up-to-date"
                          " with 'pip install -U .'")

    def run_proc_pipeline(self, subject_list=None, n_threads=-1):
        """Runs metric_proc_pipeline(), which generates measured metric values, in parallel using n_threads. Default -1
        allocates all possible threads on the system, as detected by multiprocessing.cpu_count(). Can be used to process
        normative data, or a set of subjects to plot against a trained dataset."""
        if n_threads == -1:
            import multiprocessing
            n_threads = multiprocessing.cpu_count()
        if subject_list is None:
            subject_list = self.subject_list
        # Split subjects and run parallel processing
        for group in [subject_list[i:i+n_threads] for i in range(0,len(subject_list),n_threads)]:
            self.metric_proc_pipeline.run(group)

    def load_proc_metrics(self, subject_list=None):
        """Load metrics computed by processing pipeline. Fills with Nans the values that don't exist for a given subject.
        Can be saved with save_proc_metrics() (eg a 75x1258 matrix requires 760 kB). Covariates is a list of variable
        names in participants.tsv to keep and save in a covariates numpy array.
        TODO: should be used to load both the dev_model data to fit the model, and the SOM data to compare to the model
        :param subject_list:
        :return:
        """
        if subject_list is None:
            subject_list = list(self.subject.keys())
        self.metric_names.append('age')
        self.metric_names, self.measured_metrics = self.metric_proc_pipeline.proc2metric(os.path.join(self.bids_database,'derivatives','scanometrics', self.metric_proc_pipeline.proc_pipeline_name, 'stats2table'), subject_list)

    ###################
    # NORMATIVE MODEL #
    ###################

    def set_normative_model(self, model_name='LlocvPolynomial', dataset_id='ScanOMetrics2022', models_folder=None):
        """
        Sets normative model based on model name (defaults to LLOCV) and a training set ID (defaults to ScanOMetrics2022).
        Intended to initialize normModel dictionary for further fitting or loading already trained model.

        :param model_name:
        :param model_dataset:
        :param models_folder:
        """
        if hasattr(normative_models, model_name):
            self.normativeModel = getattr(normative_models, model_name)(dataset_id, models_folder)
        else:
            logging.ERROR('Model %s not available in scanometrics.normative_models')


    def load_normative_model(self):
        """
        Load model parameters from model's data folder.
        """
        if self.normativeModel is None:
            logging.ERROR("Normative model is not set. Run 'set_normative_model()' before loading.")
        else:
            self.normativeModel.load()


    def train_normative_model(self, bids_folder=None):
        """
        Train model parameters from a bids dataset (defaults to ScanOMetrics_project bids folder).
        """
        if 'model_name' not in self.normativeModel.keys():
            logging.ERROR("Normative model is not set. Run 'set_normative_model()' before training.")
        # Load data

    ######################
    # RUN CORE FUNCTIONS #
    ######################

    def flag_outliers(self, k):
        """
	Label subjects as outlier if morphological value :math:`x\\not\\in[q_{25}-k*IQR;q_{75}+k*IQR]`, where :math:`q_{25}` and :math:`q_{75}`
        are the 25th and 75th percentiles, IQR is the interquartile range, and k sets the threshold of how many IQRs are
        considered for labeling outliers. Subjects are compared to their age matching group [0.9*age,1.1*age]. The
        'metrics' matrix is a NxM matrix with N subjects and M metrics (eg self.measured_metrics or
        self.normativeModel.residuals). The quantiles are computed using numpy.quantile and the 'hazen' method, to
        obtain the same results as the Matlab 'quantile' function.

	:param k: factor of IQRs to be used as threshold for outlier detection.
	:type k: float > 0.
        """
        if k < 0:
            logging.ERROR('Parameter k in flag_outliers() is negative, should be >= 0.')
        if self.measured_metrics.shape == (0,):
            logging.ERROR('SOM project measured metrics is empty. Load metrics before running flag_outliers()')  # Adapt to residue test

        # Recover age vector
        age = self.covariate_values[:, self.covariate_names.index('age')]
        # metrics is shape (N, M)
        # outliers = subjects['outliers'].copy()  # outliers has shape NxM_ (an (N,1) outlier array for each M_ selected metrics)
        #                                         however this means that outlier matrix depends on estimate/residual flag
        #                                         might not need to be copied if initalized here, but consider refactoring
        #                                         if needed
        self.pout = np.zeros(self.measured_metrics.shape[1])
        self.fout = np.zeros(self.measured_metrics.shape[1])
        self.part = np.ones(self.measured_metrics.shape[1])
        self.odds = np.full(self.measured_metrics.shape[1], np.inf)
        # If outliers has shape 0, meaning it's the first time flag_outliers is called, initialize array with X.shape
        if self.outliers.shape == (0,):
            self.outliers = np.zeros(self.measured_metrics.shape).astype('bool')
        for m in range(self.measured_metrics.shape[1]):  # loop through self.measured_metrics m
            thr_lrg = np.zeros(self.measured_metrics.shape[0])
            thr_sml = np.zeros(self.measured_metrics.shape[0])
            for s in range(self.measured_metrics.shape[0]):  # loop through subjects s
                if self.outliers[s, m]:  # skip current subject if already flagged as outlier
                    continue
                age_matches = ((age >= 0.9 * age[s]) & (age <= 1.1 * age[s]))  # First set age_matches within age range
                if age_matches.sum() > 2:  # age_matches always includes current subject, make sure there's at least two more
                    age_matches[s] = False  # remove current subject from matching array
                    q_25, q_75 = np.quantile(self.measured_metrics[age_matches, m], [0.25, 0.75], method='hazen')
                    IQR = q_75 - q_25
                    thr_lrg[s] = q_75 + k * IQR
                    thr_sml[s] = q_25 - k * IQR
                    if self.measured_metrics[s, m] < thr_sml[s] or self.measured_metrics[s, m] > thr_lrg[s]:
                        self.outliers[s, m] = True
            self.fout[m] = self.outliers[:, m].sum() / self.measured_metrics.shape[0]
            if self.fout[m] > 0:
                out_width = np.mean(np.abs(self.measured_metrics[self.outliers[:, m], m]))
                self.pout[m] = norm.cdf(thr_sml.mean(), loc=0, scale=out_width) + 1 - norm.cdf(thr_lrg.mean(), loc=0,
                                                                                          scale=out_width)
                if self.pout[m] > 0:
                    self.part[m] = self.fout[m] / self.pout[m]
                    if self.part[m] > 0:
                        self.odds[m] = (1 - self.part[m]) / self.part[m]
                else:
                    self.part[m] = 1  # if self.pout is below or equal to zero,
            else:
                self.part[m] = 0
                self.odds[m] = np.inf

    def compute_uncertainty(self, frac=0.10):
        """
        Compute uncertainty over metrics in self.measured_metrics. Uncertainty is computed using repeated measures
        when available, otherwise approximates it with a given fraction of the mean across subjects (default frac=0.1).

        :param self: ScanOMmetrics project
        """

