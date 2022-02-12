import numpy as np
from scipy.signal import find_peaks
from scipy.interpolate import interp1d
from matplotlib.pyplot import *


def interp(data, t, kind='cubic'):
    """
    ('linear', 'nearest', 'zero', 'slinear', 'quadratic', 'cubic',
    'previous', 'next', where 'zero', 'slinear', 'quadratic' and 'cubic')
    """
    assert len(data) == len(t), "lengths are not equal!"
    uniq = dict()
    for k, v in zip(t, data):
        uniq[k] = v
    data = np.array(list(uniq.values()))
    t = np.array(list(uniq.keys()))
    T_fin = max(t)
    # x = np.linspace(0,T_fin, num = len(data1,25), endpoint=True)
    x = t
    f_theta = interp1d(x, data, kind=kind)
    return f_theta


def preprocess(data, t, show_plots=True, main_title="", N_skip=5, percentiles=99, kind='cubic'):
    """
        Returns spline function of filtered signal
    """
    if show_plots:
        fig, axs = subplots(1, 3, figsize=(15, 2))
        title(main_title)
    ## Дропнуть пики и интерполяцией заполнить пропуски
    height_min = np.percentile(-data, percentiles)
    height_max = np.percentile(data, percentiles)
    height_delta = np.percentile(np.abs(data[1:] - data[:-1]), percentiles)
    # print(height_min)
    # print(height_delta)
    ids1 = find_peaks(data, height=height_max, threshold=-height_delta)[0]
    ids2 = find_peaks(-data, height=height_min, threshold=-height_delta)[0]
    ids = np.concatenate([ids1, ids2])
    # print(ids)

    ## Дропнуть пики
    try:
        drop_ids = np.concatenate(
            [range(id - N_skip, id + N_skip) for id in ids if id > N_skip and id < len(data) - N_skip])
    except:
        drop_ids = []
    # print(drop_ids)
    data_clean = np.delete(data, drop_ids)
    t_clean = np.delete(t, drop_ids)

    ## Интерполировать
    Fdata = interp(data_clean, t_clean, kind)
    data_interp = Fdata(t)

    if show_plots:
        axs[0].plot(data)
        axs[0].scatter(ids, np.take(data, ids), c='r')
        axs[0].set_title('Peaks detected')
        axs[0].set_ylim(-height_min, height_max)
        axs[0].grid()

        axs[1].plot(t, data, linewidth=4)
        axs[1].plot(t, data_interp, linewidth=2)
        axs[1].set_title("Interpolated")
        axs[1].set_ylim(min(data_interp), max(data_interp))
        axs[1].grid()

        axs[2].plot(t, data - data_interp)
        axs[2].set_title("Interp. error")
        axs[2].grid()
        show()
    return Fdata
