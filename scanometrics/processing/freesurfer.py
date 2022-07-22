"""
Wrapper for freesurfer scripts to be run on <subjects_dir>/<subjid>. Includes pipeline template parameters as
proc_pipeline_name, and methods as run(), and proc2metric().
"""


from shutil import rmtree, which
import subprocess
import numpy as np
import os
from scanometrics.utils import logging
from csv import DictReader
from glob import glob


proc_pipeline_name = 'freesurfer'

def run(subject_dir, subj_id, T1_file=None):
    """Pipeline template method overrided here. Calls freesurfer processing scripts. Single subject processing to allow
    both single or multi-subject processing through loop of external calls. If T1_file is specified, checks that a .mgz
    file is present in <subject_dir>/<subj_id>/mri/orig, and stops with an error otherwise.

    :param subject_dir: path to freesurfer subject directory (eg <bids_database>/derivatives/freesurfer).
    :type subject_dir: string
    :param subj_id: participant code for subject to be analysed.
    :type subj_id: string
    """
    if T1_file is None:
        if len(glob(os.path.join(subjects_dir, subj_id, 'mri', 'orig', '*.mgz'))) == 0:
            logging.ERROR("Freesurfer's recon-all is about to get called without any .mgz file in %s. Make sure to"
                          "convert Freesurfer input file to .mgz and place it there before running recon-all." % (
                          os.path.join(subjects_dir, subj_id, 'mri', 'orig')
                          ))
    run_recon_all(subject_dir, subj_id, T1_file)
    run_stats2tables(subject_dir, subj_id)


def run_recon_all(subjects_dir, subj_id, T1_file=None):
    """Wrapper for freesurfer recon-all. Takes bids T1_file as input and overwrites previous recon-all output. If T1 is
    not specified, takes recon-all assumption that at least one .mgz T1 file is in <subj_id>/mri/orig folder. Intended
    to use in bids freesurfer derivatives folder. Assumes recon-all can be called (i.e. $FREESURFER_HOME is set
    and $FREESURFER_HOME/SetUpFreeSurfer.sh has been sourced).

    :param subjects_dir: path to freesurfer subject directory (eg <bids_database>/derivatives/freesurfer).
    :type subjects_dir: string
    :param subj_id: participant code of subject to run recon-all on.
    :type subj_id: string
    :param T1_file: path to T1 nifti file to be used as input for recon-all.
    :type T1_file: string
    """
    os_env = os.environ()
    if which('recon-all') is None:
        logging.ERROR("""'recon-all command not found. Please set $FREESURFER_HOME to your Freesurfer installation
        directory and run 'source $FREESURFER_HOME/SetUpFreeSurfer.sh'.""")
    cmd = ['recon-all',
           '-subjid', subj_id]
    if T1_file is not None and os.path.exists(os.path.join(subjects_dir, subj_id)):
        rmtree(os.path.join(subjects_dir, subj_id))
        cmd.append += ['-i', T1_file]
    cmd += ['-all',
            '>', os.path.join(subjects_dir, subj_id, 'recon-all.log'),
            '2>&1']
    os_env['SUBJECTS_DIR'] = subjects_dir
    p = subprocess.run(cmd, env=os_env)
    if p.returncode:
        logging.ERROR('recon-all failed with code %d (see %s for details).' % (p.returncode, os.path.join(subjects_dir, subj_id, 'recon-all.log')))


def run_stats2tables(subjects_dir, subj_ids, seg_metrics=['volume'],
                     parc35_metrics=['volume', 'area', 'thickness', 'thicknessstd', 'meancurv', 'gauscurv', 'foldind', 'curvind', 'pctmean'],
                     parc75_metrics=['volume', 'area', 'thickness', 'thicknessstd', 'meancurv', 'gauscurv', 'foldind', 'curvind', 'pctmean'],
                     parcSTD_metrics=['mean', 'std', 'snr']):
    """Wrapper for asegstats2table and aparcstats2table. Assumes subjects were processed with recon-all. Subjects_dir
    should be the bids main folder, containing the derivatives folder."""
    if which('asegstats2table') is None:
        logging.ERROR("""asegstats2table command not found. Please set $FREESURFER_HOME to your Freesurfer installation
        directory and run 'source $FREESURFER_HOME/SetUpFreeSurfer.sh'.""")
    os_env = os.environ()
    os_env['SUBJECTS_DIR'] = subjects_dir
    for metric in seg_metrics:
        p = subprocess.run(['asegstats2table',
                            '--subjects', ' '.join(subj_ids),
                            '--meas', metric,
                            '--all-segs',
                            '--tablefile', os.path.join(subjects_dir, 'derivatives', 'scanometrics', 'stats2table', 'freesurfer', 'aseg_stats_%s.txt' % metric)
                            ], env=os_env)
        if p.returncode:
            logging.ERROR('asegstats2table failed with code %d.' % p.returncode)
    for hemi in ['rh', 'lh']:
        for metric in parc35_metrics:
            if metric not in ['pctmean', 'lgimean']:
                p = subprocess.run(["aparcstats2table",
                                "--subjects", ' '.join(subj_ids),
                                "--parc", "aparc"
                                "--hemi", hemi,
                                "--meas", metric,
                                "--parcs-from-file", os.path.join(os.path.dirname(__file__), 'resources', 'DesikanKilliany_ROIs_select.txt'),
                                "--tablefile", os.path.join(subjects_dir, 'derivatives', 'scanometrics', 'stats2table', 'freesurfer', '%s.aparc_stats_%s.txt' % (hemi, metric))], env=os_env)
                if p.returncode:
                    logging.ERROR("aparcstats2table failed with code %d." % p.returncode)
        for metric in parc75_metrics:
            if metric not in ['pctmean', 'lgimean']:
                p = subprocess.run(['aparcstats2table',
                    '--subjects', ' '.join(subj_ids),
                    '--parc', 'aparc.a2009s'
                    '--hemi', hemi,
                    '--meas', metric,
                    '--parcs-from-file', os.path.join(os.path.dirname(__file__), 'resources', 'Destrieux_ROIs_select.txt'),
                    '--tablefile', os.path.join(subjects_dir, 'derivatives', 'scanometrics', 'stats2table', 'freesurfer', '%s.aparca2009s_stats_%s.txt' % (hemi, metric))], env=os_env)
                if p.returncode:
                    logging.ERROR("aparcstats2table failed with code %d." % p.returncode)
        # GM-WM contrast needs asegstats (table format a little different)
        for metric in parcSTD_metrics:
            p = subprocess.run(["asegstats2table",
                "--inputs", ' '.join([os.path.join(subjects_dir, subj_id, 'stats', '%s.w-g.pct.stats' % hemi) for subj_id in subj_ids]),
                "--meas", metric,
                "--tablefile", os.path.join(subjects_dir, 'derivatives', 'scanometrics', 'stats2table', 'freesurfer', '%s.aparc_stats_pct%s.txt' % (hemi, metric))], env=os_env)
            p = subprocess.run(["asegstats2table",
                "--inputs", ' '.join([os.path.join(subjects_dir, subj_id, 'stats', '%s.w-g.pct.a2009s.stats' % hemi) for subj_id in subj_ids]),
                "--meas", metric,
                "--tablefile", os.path.join(subjects_dir, 'derivatives', 'scanometrics', 'stats2table', 'freesurfer', '%s.aparca2009s_stats_pct%s.txt' % (hemi, metric))], env=os_env)
        p = subprocess.run(["asegstats2table",
            "--inputs", ' '.join([os.path.join(subjects_dir, subj_id, 'stats', '%s.pial_lgi.stats' % hemi) for subj_id in subj_ids]),
            "--meas", "mean",
            "--tablefile", os.path.join(subjects_dir, 'derivatives', 'scanometrics', 'stats2table', 'freesurfer', '%s.aparc_stats_lgimean.txt' % hemi)], env=os_env)
        p = subprocess.run(["asegstats2table",
            "--inputs", ' '.join([os.path.join(subjects_dir, subj_id, 'stats', '%s.pial_lgi.a2009s.stats' % hemi) for subj_id in subj_ids]),
            "--meas", "mean",
            "--tablefile", os.path.join(subjects_dir, 'derivatives', 'scanometrics', 'stats2table', 'freesurfer', '%s.aparca2009s_stats_lgimean.txt' % hemi)], env=os_env)


def proc2metric(stats2table_folder, subjects, seg_metrics=['volume'],
                     parc35_metrics=['volume', 'area', 'thickness', 'thicknessstd', 'meancurv', 'gauscurv', 'foldind', 'curvind', 'pctmean'],
                     parc75_metrics=['volume', 'area', 'thickness', 'thicknessstd', 'meancurv', 'gauscurv', 'foldind', 'curvind', 'pctmean'],
                     parcSTD_metrics=['mean', 'std', 'snr']):
    """
    Wip: idea is to have this function defined for each preprocessing pipeline, that somehow gets all the variables
    from the preprocessing in a list of metrics, probably a table saved as text file that can then be read by
    load_subjects(). Should be the same for all subjects, so consider setting values to Nan when not available for
    certain subjects. Would make sense to put it the processing module. A priori writing a line per metric seems to be
    the most appropriate to do within a loop, to avoid having to load all variables before writing in a single line.
    Function should check for duplicates though, so actually loading everything before hand would also be practical.
    Should not take too much memory. Computes symmetric index on the side if inside a hemi loop, and adds results at end
    :param proc:
    :return:
    """
    metric_names = []
    metric_values = []
    for metric in seg_metrics:
        with open(os.path.join(stats2table_folder, 'aseg_stats_%s.txt' % metric), 'r') as f:
            loaded_subjects = []
            tmp = []
            for row in DictReader(f, delimiter='\t'):
                tmp.append(list(row.values())[1:])  # skip participant ID
                loaded_subjects.append(list(row.values())[0])
            metric_names += ['aseg%s_%s' % (metric, m) for m in list(row.keys())[1:]]  # Read metric names, skipping 1st as corresponds to participant ID
            metric_values += [tmp[loaded_subjects.index(ID)] for ID in subjects]  # make sure order is the same as in subjects
    metric_values = np.array(metric_values, dtype='float')
    for hemi in ['rh', 'lh']:
        for metric in parc35_metrics:
            if metric not in ['pctmean', 'lgimean']:
                with open(os.path.join(stats2table_folder, '%s.aparc_stats_%s.txt' % (hemi, metric)), 'r') as f:
                    loaded_subjects = []
                    tmp = []
                    for row in DictReader(f, delimiter='\t'):
                        tmp.append(list(row.values())[1:-2])  # skip participant ID, and 'BrainSegVolNotVent', 'eTIV' values
                        loaded_subjects.append(list(row.values())[0])
                    metric_names += ['aparc_%s' % m for m in list(row.keys())[1:-2]]  # skip <hemi>.aparc.<metric> entry at start, 'BrainSegVolNotVent', 'eTIV' at the end
                    metric_values = np.hstack((metric_values, np.array([tmp[loaded_subjects.index(ID)] for ID in subjects], dtype='float')))
        for metric in parc75_metrics:
            if metric not in ['pctmean', 'lgimean']:
                with open(os.path.join(stats2table_folder, '%s.aparca2009s_stats_%s.txt' % (hemi, metric)), 'r') as f:
                    loaded_subjects = []
                    tmp = []
                    for row in DictReader(f, delimiter='\t'):
                        tmp.append(list(row.values())[1:-2])  # skip participant ID
                        loaded_subjects.append(list(row.values())[0])
                    metric_names += ['aparc.a2009s_%s' % m for m in list(row.keys())[1:-2]]  # skip <hemi>.aparc.a2009s.<metric> entry
                    metric_values = np.hstack((metric_values, np.array([tmp[loaded_subjects.index(ID)] for ID in subjects], dtype='float')))
        # GM-WM contrast needs asegstats (table format a little different, no ID at each line in CH-first for example)
        for metric in parcSTD_metrics:
            with open(os.path.join(stats2table_folder, '%s.aparc_stats_pct%s.txt' % (hemi, metric)), 'r') as f:
                tmp = []
                for row in DictReader(f, delimiter='\t'):
                    tmp.append(list(row.values())[1:])  # skip subject index as doesn't provide any info
                    # loaded_subjects.append(list(row.values())[0])  <- commented as no ID in pct files (used --inputs option in asegstats2table), using same index as other FS output tables
                metric_names += ['%s_pct%s_%s' % (hemi, metric, m) for m in list(row.keys())[1:]]  # skip <hemi>.aparc.<metric> entry
                metric_values = np.hstack((metric_values, np.array([tmp[loaded_subjects.index(ID)] for ID in subjects], dtype='float')))
            with open(os.path.join(stats2table_folder, '%s.aparca2009s_stats_pct%s.txt' % (hemi, metric)), 'r') as f:
                tmp = []
                for row in DictReader(f, delimiter='\t'):
                    tmp.append(list(row.values())[1:])  # skip participant ID
                    # loaded_subjects.append(list(row.values())[0]) <- commented as no ID in pct files
                metric_names += ['%s_pct%s_%s' % (hemi, metric, m) for m in list(row.keys())[1:]]  # skip <hemi>.aparc.<metric> entry
                metric_values = np.hstack((metric_values, np.array([tmp[loaded_subjects.index(ID)] for ID in subjects], dtype='float')))
        # Optional lgi files if processed with matlab script
        lgi_file = os.path.join(stats2table_folder, '%s.aparc_stats_lgimean.txt' % hemi)
        if os.path.exists(lgi_file):
            with open(lgi_file, 'r') as f:
                loaded_subjects = []
                tmp = []
                for row in DictReader(f, delimiter='\t'):
                    tmp.append(list(row.values())[1:])  # skip participant ID
                    loaded_subjects.append(list(row.values())[0])
                metric_names += ['%s_aparc_lgimean_%s' % (hemi, m) for m in list(row.keys())[1:]]  # skip <hemi>.aparc.<metric> entry
                metric_values = np.hstack((metric_values, np.array([tmp[loaded_subjects.index(ID)] for ID in subjects], dtype='float')))
        lgi_file = os.path.join(stats2table_folder, '%s.aparca2009s_stats_lgimean.txt' % hemi)
        if os.path.exists(lgi_file):
            with open(lgi_file, 'r') as f:
                loaded_subjects = []
                tmp = []
                for row in DictReader(f, delimiter='\t'):
                    tmp.append(list(row.values())[1:])  # skip participant ID
                    loaded_subjects.append(list(row.values())[0])
                metric_names += ['%s_aparca2009s_lgimean_%s' % (hemi, m) for m in list(row.keys())[1:]]  # skip <hemi>.aparc.<metric> entry
                metric_values = np.hstack((metric_values, np.array([tmp[loaded_subjects.index(ID)] for ID in subjects], dtype='float')))
    # Compute symmetric index
    sym_metric_names = []
    sym_metric_values = []
    for i, m in enumerate(metric_names):
        if 'lh' in m:  # lh_ label means there's a column with rh_ label that corresponds (look for lh instead of rh because rh matches entorhinal and other structures with rh in it)
            num = (metric_values[:, metric_names.index(m.replace('lh', 'rh'))] - metric_values[:, i])
            denom = (metric_values[:, metric_names.index(m.replace('lh', 'rh'))] + metric_values[:, i])
            sym_metric_values.append(num / denom)  # Sym index is (rh-lh)/(rh+lh)
            sym_metric_names.append(m.replace('lh', 'symmetryIndex'))
        elif 'Left' in m:  # Left label means there's a column with Right label that corresponds
            num = (metric_values[:, metric_names.index(m.replace('Left', 'Right'))] - metric_values[:, i])
            denom = (metric_values[:, metric_names.index(m.replace('Left', 'Right'))] + metric_values[:, i])
            sym_metric_values.append(num / denom)  # Sym index is (rh-lh)/(rh+lh)
            sym_metric_names.append(m.replace('Left', 'symmetryIndex'))
    return metric_names+sym_metric_names, np.hstack((metric_values, np.array(sym_metric_values, dtype='float').T))


def _generate_random_code(length=10, seed=None):
    """Generate random alphanumerical code, starting with a letter for proper sorting in file browsers.

    :param length: specifies length of random string, defaults to 10
    :type length: int
    :param seed: seed for random string generator, defaults to None
    :type seed: int, optional
    :return: random alphanumerical string with given length
    :rtype: string
    """

    random.seed(seed)
    return ''.join(random.choices(string.ascii_uppercase)+random.choices(string.ascii_uppercase+string.digits,k=length-1))


def _read_fs_table(bids_path, subject_code):
    fs_table = {'fs_version': np.loadtxt(os.path.join(bids_path, 'derivatives', 'freesurfer', 'sub-%s' % subject_code,
                                                      'ses-%d', 'freesurfer_version'))}
    if fs_table['version'] == '6.0.0':
        for hemi in ['lh', 'rh']:
            for ROI in ['bankssts', 'caudalanteriorcingulate', 'caudalmiddlefrontal', 'corpuscallosum', 'cuneus',
                        'entorhinal', 'fusiform', 'inferiorparietal', 'inferiortemporal', 'isthmuscingulate',
                        'lateraloccipital', 'lateralorbitofrontal', 'lingual', 'medialorbitofrontal', 'middletemporal',
                        'parahippocampal', 'paracentral', 'parsopercularis', 'parsorbitalis', 'parstriangularis',
                        'pericalcarine', 'postcentral', 'posteriorcingulate', 'precentral', 'precuneus',
                        'rostralanteriorcingulate', 'rostralmiddlefrontal', 'superiorfrontal', 'superiorparietal',
                        'superiortemporal', 'supramarginal', 'frontalpole', 'temporalpole', 'transversetemporal',
                        'insula']:
                fs_table[ROI]
    return fs_table

