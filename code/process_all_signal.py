from ProcessingObj import ProcessingObj
import os
from split import split
from get_all_digits_in_string import get_all_digits_in_string
from AutoDict import AutoDict
from collections import OrderedDict
import pickle
import os.path
from list_files_in_dir import list_files_in_dir_sorted
from constants import *


def get_file_param(s):
    """
    Get the file parameters
    :param s:       as the string with digits
    :return:
    """
    nums = get_all_digits_in_string(s)
    q = int(nums[0])
    t_heat = int(nums[1])
    t_total = int(nums[2]) * SECONDS_IN_MIN
    return q, t_heat, int(t_total)


def process_all(df, check_if_file_exists=False, r_nominal=None, theta_o=None, theta_m=None, add_time_after_max=None, sand=True, peat=True):
    """
    List the files, set up the run and then extract all of the data.
    The data is stored as objects.  The files are stored based on date.

    This function is modified to do the sensitivity analysis so that the quantities can be passed in.
    Run this function first to obtain the function outputs.

    r_nominal
    theta_o
    theta_m
    add_time_after_max

    can be passed into the function

    The data is stored in the file specified by df
    """
    if check_if_file_exists:    # DO NOT RUN IF FILE EXISTS
        if os.path.isfile(df):
            return
    dirs = [x[0] for x in os.walk(MAIN_DIR)]
    dirs_cleaned = []
    dirs_peat = []
    dirs_sand = []
    for elem in dirs:
        if CAL_DIR_NAME in elem or BIN_DIR_NAME in elem:
            continue
        if PEAT_NAME in elem and not elem.endswith(PEAT_NAME) and peat is True:
            dirs_peat.append(elem)
        elif SAND_NAME in elem and not elem.endswith(SAND_NAME) and sand is True:
            dirs_sand.append(elem)
        dirs_cleaned.append(elem)
    dd = OrderedDict()
    dd[SAND_NAME] = dirs_sand
    dd[PEAT_NAME] = dirs_peat
    experiments = AutoDict()
    for soil_type, experiment_dir_list in dd.items():
        print('type: ', soil_type)
        for experiment_dir in experiment_dir_list:
            print('dir: ', experiment_dir)
            date = split(['/', '\\'], experiment_dir)[-1]
            print('date: ', date)
            files = list_files_in_dir_sorted(experiment_dir)  # sort the files by name (important to keep order)
            for file in files:
                print('file: ', file)
                q, t_heat, t_total = get_file_param(file)
                sp = split(['-'], file)
                start_string = sp[0] + HYPEN_STR
                additional_text = HYPEN_STR + sp[-1].replace(CSV_STR, '')
                print('[q = ', q, ', t_heat = ', t_heat,
                      ', start_string = ', start_string, ', additional_text = ', additional_text, ']')
                obj = ProcessingObj(experiment_dir + '/', start_string, additional_text, q, t_heat, t_total,
                                    r_nominal=r_nominal, theta_o=theta_o, theta_m=theta_m,
                                    add_time_after_max=add_time_after_max)
                obj.run()
                experiments[soil_type][date][file] = obj    # stored as soil_type, date, file number
    print('Saving out the object...')
    with open(df, 'wb') as fh:
        pickle.dump(experiments, fh)
    print('DONE')


def create_theta_rho_data_vec(df, df_vec, check_if_file_exists=False):
    """
    Create theta and rho data vectors and store these data vectors in a dictionary
    The process_all() function must be run initially.

    Run this function second to obtain the model outputs

    The data is read from the df file and output in the df_vec file
    :return:
    """
    if check_if_file_exists:    # DO NOT RUN IF FILE EXISTS
        if os.path.isfile(df_vec):
            return
    print('Loading pickle')
    dd = pickle.load(open(df, "rb"))
    experiments = AutoDict()
    print('DONE loading')
    theta_nom_vec = []
    rho_nom_vec = []
    theta_sig_vec = []
    rho_sig_vec = []
    alpha_det_vec = []
    alpha_sig_vec = []
    k_det_vec = []
    kdet_sig_vec = []
    kdet_linear_vec = []
    qav_vec = []
    qknown_vec = []
    experiments = AutoDict()
    for soil_type, soils in dd.items():
        for date, dates in soils.items():
            for file, obj in dates.items():
                theta_w_nom, rho_nom, theta_w_sig, rho_sig, \
                t2_trim1, delta_T2_trim1, dT_synth_nom, dT_synth_sig, \
                cut_time, t2_trim1_cut, r_t_heating_cooling, \
                alpha_det, alpha_sig, k_det, kdet_sig, \
                t1_trim_heating1, delta_T1_trim_heating1, dT_synth_sp, \
                delta_T1_trim, t1_trim, \
                tlinear, kdet_linear, qav = obj.get_outputs()
                """
                theta_w_nom, rho_nom        for the nominal method 
                theta_w_sig, rho_sig        for the signal processing method
                """
                print('----------------------------------------------------')
                print('file  = ', file)
                print('soil_type = ', soil_type, 'date = ', date)
                print('theta_w_nom = ', theta_w_nom, 'rho_nom = ', rho_nom)
                print('theta_w_sig = ', theta_w_sig, 'rho_sig = ', rho_sig)
                print('----------------------------------------------------')
                theta_nom_vec.append(theta_w_nom)
                rho_nom_vec.append(rho_nom)
                theta_sig_vec.append(theta_w_sig)
                rho_sig_vec.append(rho_sig)
                alpha_det_vec.append(alpha_det)
                alpha_sig_vec.append(alpha_sig)
                k_det_vec.append(k_det)
                kdet_sig_vec.append(kdet_sig)
                kdet_linear_vec.append(kdet_linear)
                qav_vec.append(qav)
                qknown_vec.append(obj.q)
                # done soil type, store the data
        print('storing for soil type: ', soil_type)
        experiments[soil_type] = (np.asarray(theta_nom_vec),
                                  np.asarray(rho_nom_vec),
                                  np.asarray(theta_sig_vec),
                                  np.asarray(rho_sig_vec),
                                  np.asarray(alpha_det_vec),
                                  np.asarray(alpha_sig_vec),
                                  np.asarray(k_det_vec),
                                  np.asarray(kdet_sig_vec),
                                  np.asarray(kdet_linear_vec),
                                  np.asarray(qav_vec),
                                  np.asarray(qknown_vec))
        theta_nom_vec = []
        rho_nom_vec = []
        theta_sig_vec = []
        rho_sig_vec = []
        alpha_det_vec = []
        alpha_sig_vec = []
        k_det_vec = []
        kdet_sig_vec = []
        kdet_linear_vec = []
        qav_vec = []
        qknown_vec = []
    # save the experiments out
    print('Saving vectors out...')
    pickle.dump(experiments, open(df_vec, 'wb'))
    print('DONE')


def main():
    process_all(DATA_PICKLE_FILE)
    create_theta_rho_data_vec(DATA_PICKLE_FILE, DATA_PICKLE_FILE_VEC)
    from plots_all_theta_rho import plot_main_all_data
    plot_main_all_data(show_plot=False)


if __name__ == '__main__':
    main()

