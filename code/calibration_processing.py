import numpy as np
from data_file_loader import construct_load_file_string, load_data_file
from common_processing import compute_difference_curve
from average_downsample import average_downsample
from butterworth_low import butter_lowpass
from scipy.signal import filtfilt
from inverse_model_nominal import obtain_sp_vars_from_curve, obtain_r_from_curve
from forward_model_nominal import compute_single_probe_forward_model, compute_delta_T_dual_heating_cooling
from get_size import length
from t_student import FindPlateau
from comparisons import compute_rmse, compute_mb, compute_percentage_diff
from sp_model_signal import sp_model_signal_inv, sp_model_nominal_get_b_c_d
from constants import *


def calibration_processing(path, additional_text, q, use_assumed_q, use_step_detector_q, t_cold, t_heat, t_total, fs,
                           k_assumed, use_assumed_k, rho, c, downsample, fs_down, filt,
                           filt_cutoff=None):
    """
    Function to perform calibration processing for the SP and DP methods.
    This function is used when the heat pulse probe is placed in the agar gel for initial calibration.

    :param path:            as the directory where the calibration data files are placed
    :param additional_text:
    :param q:               strength of the heat pulse                      [W m^-1]
    :param use_assumed_q:   True to use the set q that is passed in
    :param use_step_detector_q: Use step-detector q
    :param t_cold:          time at beginning before the heat pulse         [s]
    :param t_heat:          time of the heat pulse                          [s]
    :param t_total:         total time of the experiment                    [s]
    :param fs:              sample rate                                     [Hz]
    :param k_assumed:       assumed thermal conductivity                    [W m^-1 K^-1]
    :param use_assumed_k:   True to use the assumed k in the calculations rather than obtaining k
    :param rho:             assumed density                                 [kg m^-3]
    :param c:               heat capacity                                   [J kg^-1 K^-1]
    :param downsample:      "sp", "dp", "none" or "both" to downsample
    :param fs_down:         Downsampled frequency if downsampling is done   [Hz]
    :param filt:            "sp", "dp", "none" or "both" to filter using a butterworth filter
    :param filt_cutoff:     as the cutoff frequency of the Butterworth filter [Hz]
    :return:
    """
    # Run the initial operations
    start_string = CAL_BEGIN
    I, Rw, T1, T2, Vout, dV, dac_code, delta_T1, delta_T1_trim, delta_T1_trim_heating, delta_T2, \
        delta_T2_trim, num, qav, qprime, step_q_tfirst, step_q_tsecond, t1, t1_trim, t1_trim_heating, t2, t2_trim, \
        ypk, t,  rmse_q_calc, mb_q_calc, pd_q_calc, t2_trim_heating, delta_T2_trim_heating = \
        initial_compute_hpp(start_string, additional_text,
        downsample, filt, filt_cutoff, fs, fs_down, path, q, t_cold, t_heat, t_total, use_assumed_q,
        use_step_detector_q)

    # SINGLE PROBE
    # determine k from the single probe
    t1_trim_heating0 = t1_trim_heating - t1_trim_heating[0]     # time must start near 0
    t1_trim_heating1 = t1_trim_heating0[1:]
    delta_T1_trim_heating1 = delta_T1_trim_heating[1:]
    dt = t1_trim_heating0[1] - t1_trim_heating0[0]
    # determine the k from signal processing
    k_determined, b_determined = sp_model_signal_inv(delta_T1_trim_heating1, t1_trim_heating1, dt, qav)
    # k is used as starting value for curve-fitting to provide inputs for SP forward model
    kd, bd, cd, dd = obtain_sp_vars_from_curve(q, t1_trim_heating1, delta_T1_trim_heating1, k_determined)

    # obtain the synthetic curve only for the heating section of the single probe
    delta_T_single_synth = compute_single_probe_forward_model(q, kd, t1_trim_heating1,
                                                              bd,
                                                              cd,
                                                              dd)
    rmse_single = compute_rmse(delta_T1_trim_heating1, delta_T_single_synth)
    mb_single = compute_mb(delta_T1_trim_heating1, delta_T_single_synth)
    pd_single = compute_percentage_diff(delta_T1_trim_heating1, delta_T_single_synth)

    if use_assumed_k:
        k = k_assumed
    else:
        k = k_determined

    # COMPUTE DUAL PROBE
    rstart = R0_START                   # starting value as design specifications for system
    t2_trim0 = t2_trim - t2_trim[0]     # time must start near zero
    t2_trim1 = t2_trim0[1:]
    delta_T2_trim1 = delta_T2_trim[1:]
    alpha = k / (rho * c)
    # radius that is used for output comparison using the signal processing
    rfound = obtain_r_from_curve(q, k, alpha, t2_trim1, t_heat, delta_T2_trim1, rstart)
    # obtain the synthetic curve for comparison (heating and cooling)
    # NOTE the use of kd as determined by curve-fitting
    alpha_d = k / (rho * c)
    rfound_d = obtain_r_from_curve(q, k, alpha_d, t2_trim1, t_heat, delta_T2_trim1, rstart)
    deltaT_dual_synth = compute_delta_T_dual_heating_cooling(q, k, alpha_d, rfound_d, t2_trim1, t_heat)

    rmse_dual = compute_rmse(delta_T2_trim1, deltaT_dual_synth)
    mb_dual = compute_mb(delta_T2_trim1, deltaT_dual_synth)
    pd_dual = compute_percentage_diff(delta_T2_trim1, deltaT_dual_synth)

    # return the calculated outputs for this function
    original = (num, Vout, dV, I, Rw, T1, T2, dac_code, qprime)  # original data from the file
    """
    OUTPUTS
    
    delta_T1    = change in T1 for single probe
    delta_T2    = change in T2 for dual probe
    t1          = time vector for delta_T1
    t2          = time vector for delta_T2
       
    delta_T1_trim  = trimmed T1 curve for single probe with cold parts of curve removed
    delta_T2_trim  = trimmed T2 curve for dual probe with cold parts of curve removed 
    t1_trim        = trimmed t1 curve for single probe with cold parts of the curve removed
    t2_trim        = trimmed t2 curve for dual probe with cold parts of the curve removed
    
    delta_T1_trim_heating   = heating part of the single probe curve
    t1_trim_heating         = time vector for heating part of the single probe curve
    
    qav                     = average value of q either found or assumed
    
    
    step_q_tfirst           = first time found by the step detector
    step_q_tsecond          = second time found by the step detector
    
    ypk                     = peak sequence found by the step detector
    
    k                       = thermal conductivity (assumed or determined by the single probe)
    
    alpha                   = calculated value of alpha
    
    rfound                  = dual probe found radius (m)
    NOTE THAT rfound NEEDS TO BE CONVERTED TO UNITS OF MM
    
    delta_T_single_synth   = synthetic curve for the single probe (plot with t1_trim)
    deltaT_dual_synth      = synthetic curve for the dual probe (plot with t2_trim)
    
    t                      = Timestep for original sequence (s)
    
    """
    return original, delta_T1, delta_T2, t1, t2, delta_T1_trim, delta_T2_trim, t1_trim, t2_trim, \
        delta_T1_trim_heating, t1_trim_heating, qav, step_q_tfirst, step_q_tsecond, \
        ypk, k, alpha, rfound, delta_T_single_synth, deltaT_dual_synth, t, rmse_single, mb_single, pd_single, \
        rmse_dual, mb_dual, pd_dual, rmse_q_calc, mb_q_calc, pd_q_calc, t2_trim_heating, delta_T2_trim_heating


def initial_compute_hpp(start_string, additional_text,
                        downsample, filt, filt_cutoff, fs, fs_down, path, q, t_cold, t_heat, t_total, use_assumed_q,
                        use_step_detector_q, return_q_full=False):
    """
    Run the initial computations for the HPP.
    This code is shared between the various methods.

    :param path:            as the directory where the calibration files are placed
    :param q:               strength of the heat pulse                      [W m^-1]
    :param use_assumed_q:   True to use the set q that is passed in
    :param use_step_detector_q: Use step-detector q
    :param t_cold:          time at beginning before the heat pulse         [s]
    :param t_heat:          time of the heat pulse                          [s]
    :param t_total:         total time of the experiment                    [s]
    :param fs:              sample rate                                     [Hz]
    :param downsample:      "sp", "dp", "none" or "both" to downsample
    :param fs_down:         Downsampled frequency if downsampling is done for a specific sequence [Hz]
    :param filt:            "sp", "dp", "none" or "both" to filter using a butterworth filter
    :param filt_cutoff:     as the cutoff frequency of the Butterworth filter [Hz]
    :return:
    """
    t_total_min = int(np.ceil(t_total / SECONDS_IN_MIN))
    fn = construct_load_file_string(path, start_string, q, t_heat, t_total_min, additional_text)
    num, Vout, dV, I, Rw, T1, T2, dac_code, qprime = load_data_file(fn)
    n = length(num)  # total number of elements
    dt = 1 / fs
    # filter the input data using a Butterworth filter if required
    T1_in0 = T1
    T2_in0 = T2
    if filt != NO_CURVE:
        b, a = butter_lowpass(filt_cutoff, fs, order=BUTTERWORTH_ORDER)
        if filt == SP_CURVE:
            T1_in0 = filtfilt(b, a, T1_in0)
        elif filt == DP_CURVE:
            T2_in0 = filtfilt(b, a, T2_in0)
        elif filt == BOTH_CURVE:
            T1_in0 = filtfilt(b, a, T1_in0)
            T2_in0 = filtfilt(b, a, T2_in0)
        else:
            raise ValueError("calibration_processing: filt is not a recognized input")
    # downsample the input data as required
    T1_in1 = T1_in0
    T2_in1 = T2_in0
    T1_fs = fs
    T2_fs = fs
    if downsample != NO_CURVE:
        if filt_cutoff is None:
            raise ValueError("calibration_processing: filt_cutoff must be specified")
        if downsample == SP_CURVE:
            T1_in1 = average_downsample(fs, fs_down, T1_in1)
            T1_fs = fs_down
        elif downsample == DP_CURVE:
            T2_in1 = average_downsample(fs, fs_down, T2_in1)
            T2_fs = fs_down
        elif downsample == BOTH_CURVE:
            T1_in1 = average_downsample(fs, fs_down, T1_in1)
            T2_in1 = average_downsample(fs, fs_down, T2_in1)
            T1_fs = fs_down
            T2_fs = fs_down
        else:
            raise ValueError("calibration_processing: downsample is not a recognized input")
    # compute the difference curves
    delta_T1 = compute_difference_curve(T1_in1, T1_fs, t_cold)
    delta_T2 = compute_difference_curve(T2_in1, T2_fs, t_cold)
    # compute the time vectors for each curve
    t1n = int(np.ceil(t_total * T1_fs))
    t2n = int(np.ceil(t_total * T2_fs))
    t1 = np.linspace(0, t_total-dt, t1n)
    t2 = np.linspace(0, t_total-dt, t2n)
    # trim the cold parts of the curves and the associated time vectors
    n1_first = int(np.ceil(t_cold * T1_fs))
    n2_first = int(np.ceil(t_cold * T2_fs))
    delta_T1_trim = delta_T1[n1_first:]
    t1_trim = t1[n1_first:]
    delta_T2_trim = delta_T2[n2_first:]
    t2_trim = t2[n2_first:]

    # trim only the heating part of the curve for T1 (single probe)
    # NOTE that the trim is done after the calculation of the cold part of the curve
    n1_second = int(np.ceil(t_heat * T1_fs))
    delta_T1_trim_heating = delta_T1_trim[0:n1_second]
    t1_trim_heating = t1_trim[0:n1_second]

    # trim only the heating part of the curve for T2 (dual probe)
    # NOTE that the trim is done after the calculation of the cold part of the curve
    n1_second_dp = int(np.ceil(t_heat * T2_fs))
    delta_T2_trim_heating = delta_T2_trim[0:n1_second_dp]
    t2_trim_heating = t2_trim[0:n1_second_dp]

    # timestep for original sequence
    t = np.linspace(0, t_total-dt, n)
    step_q_tfirst = None
    step_q_tsecond = None
    ypk = None
    if use_step_detector_q:
        # find qav only over the step plateau
        ws = FWIN_PLAT_WIN
        p = FWIN_P
        ps = FindPlateau(p, ws)
        ypk, i_first, i_second = ps.find_plateau(qprime)
        step_q_tfirst = i_first * dt
        step_q_tsecond = i_second * dt
        qav = np.average(qprime[i_first:i_second])
        rmse_q_calc = compute_rmse(q, qprime[i_first:i_second])
        mb_q_calc = compute_mb(q, qprime[i_first:i_second])
        pd_q_calc = compute_percentage_diff(q, qprime[i_first:i_second])
        qfull = qprime[i_first:i_second]
    else:
        # find qav as an average over the normal time of heating
        n_cold = int(np.ceil(t_cold * fs))
        n_hot = int(np.ceil(t_heat * fs))
        nn = n_cold + n_hot
        tq_cut = t[n_cold:nn]
        q_cut = qprime[n_cold:nn]
        qav = np.average(q_cut)
        rmse_q_calc = compute_rmse(q, q_cut)
        mb_q_calc = compute_mb(q, q_cut)
        pd_q_calc = compute_percentage_diff(q, q_cut)
        qfull = q_cut
    if use_assumed_q:
        # use only the assumed q
        qq = q
    else:
        # use the averaged q determined
        qq = qav
    # Return the calculations from this function
    if return_q_full:
        return I, Rw, T1, T2, Vout, dV, dac_code, delta_T1, delta_T1_trim, delta_T1_trim_heating, \
               delta_T2, delta_T2_trim, num, qav, qprime, step_q_tfirst, step_q_tsecond, t1, t1_trim, \
               t1_trim_heating, t2, t2_trim, ypk, t, rmse_q_calc, mb_q_calc, pd_q_calc, t2_trim_heating, \
               delta_T2_trim_heating, qfull  # also return qfull
    else:
        return I, Rw, T1, T2, Vout, dV, dac_code, delta_T1, delta_T1_trim, delta_T1_trim_heating, \
               delta_T2, delta_T2_trim, num, qav, qprime, step_q_tfirst, step_q_tsecond, t1, t1_trim, \
               t1_trim_heating, t2, t2_trim, ypk, t, rmse_q_calc, mb_q_calc, pd_q_calc, t2_trim_heating, \
               delta_T2_trim_heating



