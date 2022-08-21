from scipy import interpolate
import numpy as np

## Interpolate survey sensitivity for given arguments s and q
#  @param log_s_list Array with known log(s) (horizontal table component)
#  @param log_q_list Array with known log(q) (vertical table component)
#  @param sens_list Values array [log q, log s]
#  @param s s argument for searched value
#  @param q q argument for searched value
#  @param rs Points used to interpolate in one dimension (i.e. total points number will be rs**2). It has to at least 2 for linear, 4 for cubic and 8 for quintic approximation. 
#  @param ip_kind Kind of an interpolation - 'linear', 'cubic', 'quintic'. Default: 'linear' 
def InterpolateSens(log_s_list, log_q_list, sens_list, s, q, rs=2, ip_kind='linear'):
    min_diff = 1000.0 # something big
    min_s_diff_index = 0
    min_q_diff_index = 0
    for i in range(0, len(log_s_list)):
        diff = np.abs( log_s_list[i] - np.log10(s) )
        if diff < min_diff:
            min_diff = diff
            min_s_diff_index = i

    rs_2 = rs / 2

    s_indexes = 0

    if np.log10(s) < log_s_list[min_s_diff_index] and min_s_diff_index != 0:
        min_s_diff_index -= 1
    
    if min_s_diff_index >= rs_2 - 1 and min_s_diff_index < len(log_s_list) - rs_2:
        s_indexes = np.arange(min_s_diff_index - (rs_2 - 1), min_s_diff_index + rs_2 + 1).astype(int)
    elif min_s_diff_index < rs_2 - 1:
        mov = rs_2 - min_s_diff_index - 1
        s_indexes = np.arange(min_s_diff_index - (rs_2 - 1) + mov, min_s_diff_index + rs_2 + mov + 1).astype(int)
    elif min_s_diff_index >= len(log_s_list) - rs_2:
        mov = rs_2 - (len(log_s_list) - min_s_diff_index - 1)
        s_indexes = np.arange(min_s_diff_index - (rs_2 - 1) - mov, min_s_diff_index + rs_2 - mov + 1).astype(int)

    min_diff = 1000.0
    for i in range(0, len(log_q_list)):
        diff = np.abs( log_q_list[i] - np.log10(q) )
        if diff < min_diff:
            min_diff = diff
            min_q_diff_index = i

    if np.log10(q) < log_s_list[min_q_diff_index] and min_q_diff_index != 0:
        min_q_diff_index -= 1
    
    if min_q_diff_index >= rs_2 - 1 and min_q_diff_index < len(log_q_list) - rs_2:
        q_indexes = np.arange(min_q_diff_index - (rs_2 - 1), min_q_diff_index + rs_2 + 1).astype(int)
    elif min_s_diff_index < rs_2 - 1:
        mov = rs_2 - min_q_diff_index - 1
        q_indexes = np.arange(min_q_diff_index - (rs_2 - 1) + mov, min_q_diff_index + rs_2 + mov + 1).astype(int)
    elif min_s_diff_index >= len(log_s_list) - rs_2:
        mov = rs_2 - (len(log_s_list) - min_s_diff_index - 1)
        q_indexes = np.arange(min_q_diff_index - (rs_2 - 1) - mov, min_q_diff_index + rs_2 - mov + 1).astype(int)         
        
    ss = []
    qq = []
    zz = []

    for si in s_indexes:
        for qi in q_indexes:
            ss.append(log_s_list[si])
            qq.append(log_q_list[qi])
            zz.append(sens_list[qi, si])

    ss = np.array(ss)
    qq = np.array(qq)
    zz = np.array(zz)  

    func = interpolate.interp2d(ss,qq,zz, kind=ip_kind)

    return func(np.log10(s),np.log10(q))
