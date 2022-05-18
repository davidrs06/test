"""
Wrapper for freesurfer scripts to be run on <subjects_dir>/<subjid>
"""


from shutil import rmtree, which
import subprocess
import os
from scanometrics.utils import logging

def run_recon_all(subjects_dir, subj_id, T1_file):
    """Wrapper for freesurfer recon-all. Takes bids T1_file as input and overwrites previous recon-all output.
    Intended to use in bids freesurfer derivatives folder.

    :param subjects_dir: path to freesurfer subject directory (eg <bids_database>/derivatives/freesurfer).
    :type subjects_dir: string
    :param subj_id: participant code of subject to run recon-all on.
    :type subj_id: string
    :param T1_file: path to T1 nifti file to be used as input for recon-all.
    :type T1_file: string
    """
    os_env = os.environ.copy()
    os_env['SUBJECTS_DIR'] = subjects_dir
    if which('mri_convert') is None:
        logging.ERROR("""'mri_convert command not found. Please set $FREESURFER_HOME to your Freesurfer installation
        directory and run 'source $FREESURFER_HOME/SetUpFreeSurfer.sh'.""")
    if os.path.exists(os.path.join(subjects_dir, subj_id)):
        rmtree(os.path.join(subjects_dir, subj_id))
    os.makedirs(os.path.join(subjects_dir, subj_id, 'mri', 'orig'))
    cmd = ["mri_convert", T1_file, os.path.join(subjects_dir, subj_id, 'mri', 'orig', '001.mgz')]
    p = subprocess.run(cmd, env=os_env)
    if p.returncode:
        logging.ERROR('Command failed : %s' % (' '.join(cmd)))
    os.makedirs(os.path.join(subjects_dir, subj_id, 'scripts'))
    with open(os.path.join(subjects_dir, subj_id, 'scripts', 'recon-all.log'), 'w') as fs_log_file:
        p = subprocess.run(['recon-all',
                            '-subjid', subj_id,
                            '-all'], env=os_env, stdout=fs_log_file, stderr=fs_log_file)
        if p.returncode:
            logging.ERROR('recon-all failed with code %d (see %s for details).' % (p.returncode, os.path.join(subjects_dir, subj_id, 'scripts', 'recon-all.log')))

def run_stats2tables(subjects_dir, subj_id, seg_metrics = ['volume'],
                     parc35_metrics = ['volume', 'area', 'thickness', 'thicknessstd', 'meancurv', 'gauscurv', 'foldind', 'curvind', 'pctmean'],
                     parc75_metrics = ['volume', 'area', 'thickness', 'thicknessstd', 'meancurv', 'gauscurv', 'foldind', 'curvind', 'pctmean'],
                     parcSTD_metrics = ['mean', 'std', 'snr']):
    """Wrapper for asegstats2table and aparcstats2table. Assumes subject was processed with recon-all."""
    if which('asegstats2table') is None:
        logging.ERROR("""asegstats2table command not found. Please set $FREESURFER_HOME to your Freesurfer installation
        directory and run 'source $FREESURFER_HOME/SetUpFreeSurfer.sh'.""")
    os_env = os.environ.copy()
    os_env['SUBJECTS_DIR'] = subjects_dir
    for metric in seg_metrics:
        p = subprocess.run(['asegstats2table',
                            '--subjects', subj_id,
                            '--meas', metric,
                            '--all-segs',
                            '--tablefile', os.path.join(subjects_dir, subj_id, 'stats2table', 'aseg_stats_%s.txt' % metric)
                            ], env=os_env)
        if p.returncode:
            logging.ERROR('asegstats2table failed with code %d.' % p.returncode)
    for hemi in ['rh','lh']:
        for metric in parc35_metrics:
            if metric not in ['pctmean', 'lgimean']:
                p = subprocess.run(["aparcstats2table",
                                "--subjects", subj_id,
                                "--parc", "aparc"
                                "--hemi", hemi,
                                "--meas", metric,
                                "--parcs-from-file", os.path.join(os.path.dirname(__file__), 'resources', 'DesikanKilliany_ROIs_select.txt'),
                                "--tablefile", os.path.join(subjects_dir, subj_id, 'stats2table', '%s.aparc_stats_%s.txt' % (hemi, metric))], env=os_env)
                if p.returncode:
                    logging.ERROR("aparcstats2table failed with code %d." % p.returncode)
        for metric in parc75_metrics:
            if metric not in ['pctmean', 'lgimean']:
                p = subprocess.run(['aparcstats2table',
                    '--subjects', subj_id,
                    '--parc', 'aparc.a2009s'
                    '--hemi', hemi,
                    '--meas', metric,
                    '--parcs-from-file', os.path.join(os.path.dirname(__file__), 'resources', 'Destrieux_ROIs_select.txt'),
                    '--tablefile', os.path.join(subjects_dir, subj_id, 'stats2table', '%s.aparca2009s_stats_%s.txt' % (hemi, metric))], env=os_env)
                if p.returncode:
                    logging.ERROR("aparcstats2table failed with code %d." % p.returncode)
        # GM-WM contrast needs asegstats (table format a little different)
        for metric in parcSTD_metrics:
            p = subprocess.run(["asegstats2table",
                "--inputs", os.path.join(subjects_dir, subj_id, 'stats', '%s.w-g.pct.stats' % hemi),
                "--meas", metric,
                "--tablefile", os.path.join(subjects_dir, subj_id, 'stats2table', '%s.aparc_stats_pct%s.txt' % (hemi, metric))], env=os_env)
            p = subprocess.run(["asegstats2table",
                "--inputs", os.path.join(subjects_dir, subj_id, 'stats', '%s.w-g.pct.a2009s.stats' % hemi),
                "--meas", metric,
                "--tablefile", os.path.join(subjects_dir, subj_id, 'stats2table', '%s.aparca2009s_stats_pct%s.txt' % (hemi, metric))], env=os_env)
        p = subprocess.run(["asegstats2table",
            "--inputs", os.path.join(subjects_dir, subj_id, 'stats', '%s.pial_lgi.stats' % hemi),
            "--meas", "mean",
            "--tablefile", os.path.join(subjects_dir, subj_id, 'stats2table', '%s.aparc_stats_lgimean.txt' % hemi)], env=os_env)
        p = subprocess.run(["asegstats2table",
            "--inputs", os.path.join(subjects_dir, subj_id, 'stats', '%s.pial_lgi.a2009s.stats' % hemi),
            "--meas", "mean",
            "--tablefile", os.path.join(subjects_dir, subj_id, 'stats2table', '%s.aparca2009s_stats_lgimean.txt' % hemi)], env=os_env)
