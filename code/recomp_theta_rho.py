from obtain_load_r0 import SaveLoadR0
from calibration_processing import initial_compute_hpp
from curve_time_trim import curve_time_trim
from sp_model_signal import sp_model_signal_inv
from inverse_model_nominal import obtain_k_alpha_from_dual_probe
from forward_model_nominal import compute_delta_T_dual_heating_cooling
from get_max import get_max
from dp_model_recomp import dp_model_recomp_heating_cooling, dp_model_recomp_heating_cooling_trim
from get_theta_rho import get_theta_rho
from inverse_model_nominal import obtain_sp_vars_from_curve
from forward_model_nominal import compute_single_probe_forward_model
from sp_model_signal import sp_model_late_time_inv
from constants import *


def recomp_theta_rho_additional(path, start_string, additional_text, q, t_heat, t_total, end_trim_time,
                                trim_sp_begin=0,
                                r_nominal=None, theta_o=None, theta_m=None, add_time_after_max=None):
    """
    Open a file and recompose theta, rho and other inputs
    :param path:                        as the path where to do the processing
    :param start_string:                as the starting string of the file
    :param additional_text:             as the additional text at the end
    :param q:                           as the heat input into the soil
    :param t_heat:                      as the time of heating
    :param t_total:                     as the total time
    :param end_trim_time:               as the time to trim at the beginning of the sequence
    :param trim_sp_begin:               as the number of points to trim at the beginning of the SP sequence
                                        for inverse model
    :param r_nominal:                   as the nominal radius (to override the calibrated radius)
    :param theta_o:                     fraction of organic content in the soil
    :param theta_m:                     fraction of mineral content in the soil
    :param: add_time_after_max:         time to cut the curve after the maximum
    :return:
    """
    downsample = 'dp'
    filt = 'both'
    filt_cutoff = DP_LOWPASS_FILTER_CUT
    fs = FS_SAMPLE
    fs_down = FS_DOWN_DP
    t_cold = T_COLD_BEGIN

    if add_time_after_max is None:
        add_time_after_max = NOM_TIME_ADD_AFTER_MAX

    sr = SaveLoadR0()
    sr.load_r0(CONST_PATH + CAL_R0_FILENAME + SAVEZ_EXT)
    if r_nominal is None:
        r_nominal = sr.get_r0(q) * 1.0e-3

    # set whether we are dealing with sand or peat
    if 'sand' in start_string:
        Cm_set = None
        Co_set = None
        if theta_o is None:
            theta_o = theta_o_sand
        if theta_m is None:
            theta_m = theta_m_sand
    elif 'peat' in start_string:
        Cm_set = Cm_peat
        Co_set = Co_peat
        if theta_o is None:
            theta_o = theta_o_peat
        if theta_m is None:
            theta_m = theta_m_peat
    else:
        raise ValueError('recomp_theta_rho_additional: start_string must contain sand or peat')

    use_assumed_q = False
    use_step_detector_q = True

    I, Rw, T1, T2, Vout, dV, dac_code, delta_T1, delta_T1_trim, delta_T1_trim_heating, \
    delta_T2, delta_T2_trim, num, qav, qprime, step_q_tfirst, step_q_tsecond, t1, t1_trim, \
    t1_trim_heating, t2, t2_trim, ypk, t, rmse_q_calc, mb_q_calc, pd_q_calc, t2_trim_heating, \
    delta_T2_trim_heating = \
        initial_compute_hpp(start_string, additional_text,
        downsample, filt, filt_cutoff, fs, fs_down, path, q, t_cold, t_heat, t_total, use_assumed_q,
        use_step_detector_q, return_q_full=False)

    dt_sp = t1_trim_heating[1] - t1_trim_heating[0]
    dt_dp = t2_trim_heating[1] - t2_trim_heating[0]

    t1_trim_heating1, delta_T1_trim_heating1 = curve_time_trim(t1_trim_heating, delta_T1_trim_heating)
    t2_trim_heating1, delta_T2_trim_heating1 = curve_time_trim(t2_trim_heating, delta_T2_trim_heating)
    t2_trim1, delta_T2_trim1 = curve_time_trim(t2_trim, delta_T2_trim)

    kdet_sig, _bdet = sp_model_signal_inv(delta_T1_trim_heating1[trim_sp_begin:], t1_trim_heating1[trim_sp_begin:],
                                          dt_sp, qav,
                                          output_model_array=False)

    tlinear, kdet_linear = sp_model_late_time_inv(delta_T1_trim_heating1, t1_trim_heating1, dt_sp, qav,
                                                  entire_heating_time=False, return_all=False)

    th = t_heat
    kstart = 5
    Hstart = 80
    k_det, alpha_det = obtain_k_alpha_from_dual_probe(qav, t2_trim1, th, delta_T2_trim1, kstart, Hstart,
                                                      r_nominal, full=True)
    dT_synth_nom = compute_delta_T_dual_heating_cooling(qav, k_det, alpha_det, r_nominal, t2_trim1, th, split=False)

    k_out_sp, b_out_sp, c_out_sp, d_out_sp = obtain_sp_vars_from_curve(qav, t1_trim_heating1,
                                                                       delta_T1_trim_heating1, kdet_sig)
    dT_synth_sp = compute_single_probe_forward_model(qav, k_out_sp, t1_trim_heating1, b_out_sp, c_out_sp, d_out_sp,
                                                     constrain=True)

    idx_peak, md_peak = get_max(dT_synth_nom)
    idxd_found = idx_peak
    idxd_time = dt_dp * idxd_found
    end_cut_num = int(np.floor(fs_down * end_trim_time))

    cut_time = idxd_time + add_time_after_max  # time when the curve needs to be cut (can also be negative)
    cut_num = int(np.floor(fs_down * cut_time))
    if end_cut_num > 0:
        t2_trim1_cut = t2_trim1[cut_num:-end_cut_num]
        delta_T2_trim1_cut = delta_T2_trim1[cut_num:-end_cut_num]
    else:
        t2_trim1_cut = t2_trim1[cut_num:]
        delta_T2_trim1_cut = delta_T2_trim1[cut_num:]

    r_t_heating_cooling, gamma4 = dp_model_recomp_heating_cooling_trim(fs_down, qav, kdet_sig, t2_trim1_cut, th,
                                                                       delta_T2_trim1_cut, r_nominal, get_gamma4=True)
    alpha_vec = (r_t_heating_cooling**2) / (4.0*gamma4)
    alpha_sig = np.mean(np.abs(alpha_vec))

    theta_w_nom, rho_nom = get_theta_rho(k_det, alpha_det, theta_o, theta_m, Cm_set, Co_set)
    theta_w_sig, rho_sig = get_theta_rho(kdet_sig, alpha_sig, theta_o, theta_m, Cm_set, Co_set)

    dT_synth_sig = compute_delta_T_dual_heating_cooling(qav, kdet_sig, alpha_sig, r_nominal, t2_trim1, th, split=False)
    """
    theta_w_nom                 water content from nominal heating and cooling curve-fitting
    rho_nom                     density from nominal heating and cooling curve-fitting
    theta_w_sig                 water content from signal processing 
    rho_w_sig                   density from signal processing
        
    t2_trim1                    time vector for DP experiment
    delta_T2_trim1              temperature difference vector for DP experiment
    dT_synth_nom                synthetic curve for DP using nominal method
    dT_synth_sig                synthetic curve for DP using signal processing
    
    cut_time                    time when sequence is cut
    t2_trim1_cut                cut time sequence for heating and cooling
    r_t_heating_cooling         r(t) for heating and cooling determined over cut time
    
    alpha_det                   alpha from DP curve-fitting
    alpha_sig                   alpha from signal processing
        
    k_det                       thermal conductivity from curve-fitting
    kdet_sig                    thermal conductivity from signal processing
    
    t1_trim_heating1            trimmed SP time vector
    delta_T1_trim_heating1      trimmed SP time heating vector
    dT_synth_sp                 synthetic model from single probe
    
    delta_T1_trim               change in temperature of SP
    t1_trim                     time vector associated with SP     
    
    tlinear                     time at which curve becomes linear
    kdet_linear                 k detected using linear algorithm                   
    """
    return theta_w_nom, rho_nom, theta_w_sig, rho_sig, \
           t2_trim1, delta_T2_trim1, dT_synth_nom, dT_synth_sig, \
           cut_time, t2_trim1_cut, r_t_heating_cooling, \
           alpha_det, alpha_sig, k_det, kdet_sig, \
           t1_trim_heating1, delta_T1_trim_heating1, dT_synth_sp, \
           delta_T1_trim, t1_trim, \
           tlinear, kdet_linear, qav


