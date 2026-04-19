import numpy as np
import matplotlib.pyplot as plt
from astropy import constants as const

solar_radius = const.R_sun.value 
au = const.au.value  # Astronomical Unit in meters
k = const.k_B.value  # Boltzmann constant in J/K
c = const.c.value  # Speed of light in m/s
h = const.h.value  # Planck's constant in J*s
sigma = const.sigma_sb.value  # 5.6704e-8 W/m²/K⁴


def plot_data(
    dfs,
    x_col,
    y_col,
    yerr=None,
    xlabel=None,
    ylabel=None,
    title=None,
    labels=None,
    colors=None,
    savepath=None,
    show=True,

    # ---- style controls (preserved) ----
    title_fontsize=16,
    axis_label_fontsize=14,
    tick_fontsize=12,
    legend_fontsize=11,
    figure_dpi=100,
    figure_size=(10, 7),
    grid_alpha=0.2,
    legend_loc='best'
):
    """
    Pure plotting function with full style control.
    No model computation.
    """

    # ---------- normalize input ----------
    if not isinstance(dfs, list):
        dfs = [dfs]
    n = len(dfs)

    if labels is None:
        labels = [f'Data {i+1}' for i in range(n)]
    if len(labels) != n:
        raise ValueError("labels length mismatch")

    if colors is None:
        colors = plt.cm.tab10(np.linspace(0, 1, n))
    if len(colors) != n:
        raise ValueError("colors length mismatch")

    # ---------- error parser ----------
    def parse_err_for_df(df, err_spec):
        if err_spec is None:
            return None
        if isinstance(err_spec, str):
            return df[err_spec].to_numpy(dtype=float)
        if np.isscalar(err_spec):
            return np.full(len(df), float(err_spec))
        arr = np.asarray(err_spec, dtype=float)
        if len(arr) != len(df):
            raise ValueError("Error array length mismatch")
        return arr

    # ---------- figure ----------
    fig, ax = plt.subplots(
        figsize=figure_size,
        dpi=figure_dpi
    )
    fig.patch.set_facecolor('white')

    # ---------- styling ----------
    ax.set_xlabel(xlabel or x_col, fontsize=axis_label_fontsize)
    ax.set_ylabel(ylabel or y_col, fontsize=axis_label_fontsize)

    if title is not None:
        ax.set_title(title, fontsize=title_fontsize, pad=12)

    ax.tick_params(
        axis='both',
        labelsize=tick_fontsize,
        direction='in',
        length=5,
        width=1.1
    )

    for spine in ax.spines.values():
        spine.set_linewidth(1.1)

    ax.grid(alpha=grid_alpha, linestyle=':', linewidth=0.8)

    # ---------- plot data ----------
    for df, lab, col in zip(dfs, labels, colors):
        x = df[x_col].to_numpy(dtype=float)
        y = df[y_col].to_numpy(dtype=float)
        cur_yerr = parse_err_for_df(df, yerr)

        ax.errorbar(
            x, y,
            yerr=cur_yerr,
            fmt='o',
            capsize=3,
            color=col,
            alpha=0.85,
            markersize=5,
            elinewidth=1.1,
            capthick=1.1,
            label=lab
        )

    # ---------- legend ----------
    ax.legend(
        loc=legend_loc,
        fontsize=legend_fontsize,
        frameon=True,
        fancybox=False,
        edgecolor='black'
    )

    ax.set_ylim(0, ax.get_ylim()[1])

    fig.tight_layout()

    # ---------- output ----------
    if savepath:
        fig.savefig(savepath, dpi=figure_dpi, bbox_inches='tight', facecolor='white')

    if show:
        plt.show()

    return fig


def blackbody_plot_with_temps(
    dfs,
    x_col,
    y_col,
    temp_params,  # [min_temp, max_temp, num_temps]
    yerr=None,
    xlabel=None,
    ylabel=None,
    title='Blackbody Curves for Different Temperatures',
    labels=None,
    colors=None,
    savepath=None,
    show=True,
    show_box_zoom=False,
    box=None,
    # ---- Control de fuentes / estilo ----
    title_fontsize=16,
    axis_label_fontsize=14,
    tick_fontsize=12,
    legend_fontsize=11,
    figure_dpi=100,
    figure_size=(10, 7),
    grid_alpha=0.2,
    legend_loc='best'
):
    """
    Plot blackbody curves for different temperatures and overlay experimental data.

    temp_params: [min_temp, max_temp, num_temps]
    box: [(x1, y1), (x2, y2)]
    """

    # ------------------------- normalizar entrada --------------------------
    if not isinstance(dfs, list):
        dfs = [dfs]
    n = len(dfs)

    if labels is None:
        labels = [f'Data {i+1}' for i in range(n)]
    if len(labels) != n:
        raise ValueError("labels debe tener la misma longitud que dfs")

    if colors is None:
        colors = plt.cm.tab10(np.linspace(0, 1, n))
    if len(colors) != n:
        raise ValueError("colors debe tener la misma longitud que dfs")

    # ------------------------- temperatures --------------------------
    min_temp, max_temp, num_temps = temp_params
    temperatures = np.linspace(min_temp, max_temp, num_temps)

    # ------------------------- helpers --------------------------
    def parse_err_for_df(df, err_spec, fallback_len):
        if err_spec is None:
            return None
        if isinstance(err_spec, str):
            return df[err_spec].to_numpy(dtype=float)
        if np.isscalar(err_spec):
            return np.full(fallback_len, float(err_spec), dtype=float)
        arr = np.asarray(err_spec, dtype=float)
        if arr.shape[0] != fallback_len:
            raise ValueError("Error array debe tener la misma longitud que el dataframe")
        return arr



    def blackbody_model(wavelength_nm, T):
        """
        wavelength_nm en nm
        T en K

        devuelve irradiancia espectral en W / m^2 / nm
        """
        lam = wavelength_nm * 1e-9  # nm -> m

        exponent = (h * c) / (lam * k * T)
        exponent = np.clip(exponent, 1e-12, 700)

        B_lambda = (2 * h * c**2) / (lam**5) / (np.exp(exponent) - 1.0)
        # W / (m^2 sr m)

        F_lambda_m = B_lambda 
        # W / (m^2 m)

        F_lambda_nm = F_lambda_m * 1e-9
        # W / (m^2 nm)

        return F_lambda_nm

    if show_box_zoom:
        if box is None or len(box) != 2:
            raise ValueError("box debe ser [(x1, y1), (x2, y2)] cuando show_box_zoom=True")
        (x1, y1), (x2, y2) = box
        x_box_min, x_box_max = sorted([float(x1), float(x2)])
        y_box_min, y_box_max = sorted([float(y1), float(y2)])

    # ------------------------- figura --------------------------
    fig, ax = plt.subplots(
        figsize=figure_size,
        dpi=figure_dpi
    )
    fig.patch.set_facecolor('white')

    if show_box_zoom:
        from matplotlib.patches import Rectangle
        from mpl_toolkits.axes_grid1.inset_locator import inset_axes

    def style_axes(ax, xlabel_text, ylabel_text, title_text, t_fs, lab_fs, tick_fs):
        ax.set_xlabel(xlabel_text, fontsize=lab_fs)
        ax.set_ylabel(ylabel_text, fontsize=lab_fs)
        if title_text is not None:
            ax.set_title(title_text, fontsize=t_fs, pad=12)
        ax.tick_params(axis='both', labelsize=tick_fs, direction='in', length=5, width=1.1)
        for spine in ax.spines.values():
            spine.set_linewidth(1.1)
        ax.grid(alpha=grid_alpha, linestyle=':', linewidth=0.8)

    style_axes(
        ax,
        xlabel or x_col,
        ylabel or y_col,
        title,
        title_fontsize,
        axis_label_fontsize,
        tick_fontsize
    )

    # ------------------------- find x range --------------------------
    all_x = []
    for df in dfs:
        x = df[x_col].to_numpy(dtype=float)
        all_x.extend(x)
    min_x = np.min(all_x)
    max_x = np.max(all_x)
    x_line = np.linspace(min_x, max_x, 1000)

    # ------------------------- plot data --------------------------
    for df, lab, col in zip(dfs, labels, colors):
        x = df[x_col].to_numpy(dtype=float)
        y = df[y_col].to_numpy(dtype=float)
        cur_yerr = parse_err_for_df(df, yerr, len(df))

        ax.errorbar(
            x, y,
            yerr=cur_yerr,
            fmt='o',
            capsize=3,
            color=col,
            alpha=0.85,
            markersize=5,
            elinewidth=1.1,
            capthick=1.1,
            label=lab
        )

    # ------------------------- plot blackbody --------------------------
    for temp in temperatures:
        y_line = blackbody_model(x_line, temp)
        ax.plot(
            x_line, y_line,
            label=f'T = {temp:.0f} K'
        )

    if show_box_zoom:
        ax.add_patch(
            Rectangle(
                (x_box_min, y_box_min),
                x_box_max - x_box_min,
                y_box_max - y_box_min,
                fill=False,
                edgecolor='black',
                linewidth=1.2,
                linestyle='--'
            )
        )

        ax_zoom = inset_axes(ax, width='38%', height='38%', loc='center right')

        for df, col in zip(dfs, colors):
            x = df[x_col].to_numpy(dtype=float)
            y = df[y_col].to_numpy(dtype=float)
            ax_zoom.plot(x, y, 'o', color=col, alpha=0.85, markersize=4)

        for temp in temperatures:
            y_line = blackbody_model(x_line, temp)
            ax_zoom.plot(x_line, y_line)

        ax_zoom.set_xlim(x_box_min, x_box_max)
        ax_zoom.set_ylim(y_box_min, y_box_max)
        ax_zoom.set_title('Zoom', fontsize=max(tick_fontsize - 1, 8), pad=6)
        ax_zoom.tick_params(axis='both', labelsize=max(tick_fontsize - 2, 8), direction='in')
        ax_zoom.grid(alpha=grid_alpha, linestyle=':', linewidth=0.6)
        for spine in ax_zoom.spines.values():
            spine.set_linewidth(1.0)

    ax.legend(
        loc=legend_loc,
        fontsize=legend_fontsize,
        frameon=True,
        fancybox=False,
        edgecolor='black'
    )

    ax.set_ylim(0, ax.get_ylim()[1])

    fig.tight_layout()

    if savepath:
        fig.savefig(savepath, dpi=figure_dpi, bbox_inches='tight', facecolor='white')

    if show:
        plt.show()

    return fig

def blackbody_temperature_comparison_with_errors(
    dfs,
    x_col,
    y_col,
    yerr=None,                    # str (columna) o array o escalar o None
    xlabel=None,
    ylabel=None,
    title='Blackbody Fit with Residuals',
    labels=None,
    colors=None,
    savepath=None,
    show=True,
    temperature = 5850,

    # ---- Control de fuentes / estilo ----
    title_fontsize=16,
    axis_label_fontsize=14,
    tick_fontsize=12,
    legend_fontsize=11,
    residual_title_fontsize=None,
    residual_label_fontsize=None,
    residual_tick_fontsize=None,
    residual_legend_fontsize=None,

    figure_dpi=100,
    figure_size=(10, 7),
    grid_alpha=0.2,
    legend_loc='best'
):
    """
    Ajuste de cuerpo negro para uno o varios datasets con errores en y
    + gráfico de residuales normalizados.

    El modelo usa:
        F_lambda = B_lambda(T) * Omega_sun

    donde:
        - wavelength en nm
        - irradiance en W / m^2 / nm
        - Omega_sun se calcula a partir del diámetro angular solar:
          0°32'11.41''
    """

    # ------------------------- defaults residuales -------------------------
    if residual_title_fontsize is None:
        residual_title_fontsize = axis_label_fontsize + 1
    if residual_label_fontsize is None:
        residual_label_fontsize = axis_label_fontsize
    if residual_tick_fontsize is None:
        residual_tick_fontsize = tick_fontsize
    if residual_legend_fontsize is None:
        residual_legend_fontsize = max(8, legend_fontsize - 1)

    # ------------------------- normalizar entrada --------------------------
    if not isinstance(dfs, list):
        dfs = [dfs]
    n = len(dfs)

    if labels is None:
        labels = [f'Data {i+1}' for i in range(n)]
    if len(labels) != n:
        raise ValueError("labels debe tener la misma longitud que dfs")

    if colors is None:
        colors = plt.cm.tab10(np.linspace(0, 1, n))
    if len(colors) != n:
        raise ValueError("colors debe tener la misma longitud que dfs")

    # ------------------------- helpers --------------------------
    def parse_err_for_df(df, err_spec, fallback_len):
        if err_spec is None:
            return None
        if isinstance(err_spec, str):
            return df[err_spec].to_numpy(dtype=float)
        if np.isscalar(err_spec):
            return np.full(fallback_len, float(err_spec), dtype=float)
        arr = np.asarray(err_spec, dtype=float)
        if arr.shape[0] != fallback_len:
            raise ValueError("Error array debe tener la misma longitud que el dataframe")
        return arr

    def calculate_r_squared(y_true, y_pred):
        ss_res = np.sum((y_true - y_pred) ** 2)
        ss_tot = np.sum((y_true - np.mean(y_true)) ** 2)
        return 1 - ss_res / ss_tot if ss_tot != 0 else 0.0



    def blackbody_model(wavelength_nm, T):
        """
        wavelength_nm en nm
        T en K

        devuelve irradiancia espectral en W / m^2 / nm
        """
        lam = wavelength_nm * 1e-9  # nm -> m

        exponent = (h * c) / (lam * k * T)
        exponent = np.clip(exponent, 1e-12, 700)

        B_lambda = (2 * h * c**2) / (lam**5) / (np.exp(exponent) - 1.0)
        # W / (m^2 sr m)

        F_lambda_m = B_lambda 
        # W / (m^2 m)

        F_lambda_nm = F_lambda_m * 1e-9
        # W / (m^2 nm)

        return F_lambda_nm

    # ------------------------- figura (2 filas) --------------------------
    fig, (ax_main, ax_res) = plt.subplots(
        2, 1,
        figsize=figure_size,
        dpi=figure_dpi,
        gridspec_kw={'height_ratios': [2, 1]}
    )
    fig.patch.set_facecolor('white')

    def style_axes(ax, xlabel_text, ylabel_text, title_text, t_fs, lab_fs, tick_fs):
        ax.set_xlabel(xlabel_text, fontsize=lab_fs)
        ax.set_ylabel(ylabel_text, fontsize=lab_fs)
        if title_text is not None:
            ax.set_title(title_text, fontsize=t_fs, pad=12)
        ax.tick_params(axis='both', labelsize=tick_fs, direction='in', length=5, width=1.1)
        for spine in ax.spines.values():
            spine.set_linewidth(1.1)
        ax.grid(alpha=grid_alpha, linestyle=':', linewidth=0.8)

    style_axes(
        ax_main,
        xlabel or x_col,
        ylabel or y_col,
        title,
        title_fontsize,
        axis_label_fontsize,
        tick_fontsize
    )

    # ------------------------- loop datasets --------------------------
    all_results = []
    all_residuals = []

    for df, lab, col in zip(dfs, labels, colors):
        x = df[x_col].to_numpy(dtype=float)
        y = df[y_col].to_numpy(dtype=float)
        cur_yerr = parse_err_for_df(df, yerr, len(df))

        # filtrar válidos
        mask = np.isfinite(x) & np.isfinite(y)
        if cur_yerr is not None:
            mask &= np.isfinite(cur_yerr) & (cur_yerr > 0)

        x = x[mask]
        y = y[mask]
        if cur_yerr is not None:
            cur_yerr = cur_yerr[mask]

        if len(x) < 3:
            raise ValueError(f"Dataset '{lab}': se requieren al menos 3 puntos válidos")

        # ordenar para graficar bien
        order = np.argsort(x)
        x = x[order]
        y = y[order]
        if cur_yerr is not None:
            cur_yerr = cur_yerr[order]

        # se toma solo la temperatura y ya

        temperature = temperature
        temperature_err = 1

        y_pred = blackbody_model(x, temperature)
        r2 = calculate_r_squared(y, y_pred)

        # residuales
        if cur_yerr is not None:
            residuals = (y - y_pred) / cur_yerr
        else:
            residuals = y - y_pred

        all_residuals.append((x, residuals, lab, col))

        # --- main plot ---
        ax_main.errorbar(
            x, y,
            yerr=cur_yerr,
            fmt='o',
            capsize=3,
            color=col,
            alpha=0.85,
            markersize=5,
            elinewidth=1.1,
            capthick=1.1,
            label=lab
        )

        x_line = np.linspace(np.min(x), np.max(x), 800)
        y_line = blackbody_model(x_line, temperature)
        ax_main.plot(
            x_line, y_line,
            '--',
            color='blue',
            linewidth=1.8,
            label='Plank curve at T=' + str(temperature) + ' K'
        )

        all_results.append({
            'dataset': lab,
            'temperature': temperature,
            'temperature_err': temperature_err,
            'r_squared': r2
        })

        print(f"\n--- Results for {lab} ---")
        print(f"Estimated temperature: {temperature:.6g} ± {temperature_err:.6g} K")
        print(f"R²:                    {r2:.6g}")

    ax_main.legend(
        loc=legend_loc,
        fontsize=legend_fontsize,
        frameon=True,
        fancybox=False,
        edgecolor='black'
    )

    ax_main.set_ylim(0, ax_main.get_ylim()[1])

    # ------------------------- residuales panel --------------------------
    style_axes(
        ax_res,
        xlabel or x_col,
        'Residuals',
        'Normalized Residuals' if yerr is not None else 'Residuals',
        residual_title_fontsize,
        residual_label_fontsize,
        residual_tick_fontsize
    )

    for x_data, res, lab, col in all_residuals:
        ax_res.scatter(
            x_data, res,
            color=col,
            alpha=0.5,
            s=30,
            edgecolors='none',
            label=lab
        )

    ax_res.axhline(0, color='black', linestyle='--', linewidth=1.2, alpha=0.8)

    if len(all_residuals) > 0:
        max_abs = max(np.max(np.abs(r)) for _, r, _, _ in all_residuals if len(r) > 0)
        if np.isfinite(max_abs) and max_abs > 0:
            ax_res.set_ylim(-1.1 * max_abs, 1.1 * max_abs)

    ax_res.legend(
        loc=legend_loc,
        fontsize=residual_legend_fontsize,
        frameon=True,
        fancybox=False,
        edgecolor='black'
    )

    fig.tight_layout()

    if savepath:
        fig.savefig(savepath, dpi=figure_dpi, bbox_inches='tight', facecolor='white')

    if show:
        plt.show()

    return all_results, fig


def blackbody_fit_with_errors(
    dfs,
    x_col,
    y_col,
    yerr=None,                    # str (columna) o array o escalar o None
    xlabel=None,
    ylabel=None,
    title='Blackbody Fit with Residuals',
    labels=None,
    colors=None,
    savepath=None,
    show=True,

    # ---- Control de fuentes / estilo ----
    title_fontsize=16,
    axis_label_fontsize=14,
    tick_fontsize=12,
    legend_fontsize=11,
    residual_title_fontsize=None,
    residual_label_fontsize=None,
    residual_tick_fontsize=None,
    residual_legend_fontsize=None,

    figure_dpi=100,
    figure_size=(10, 7),
    grid_alpha=0.2,
    legend_loc='best'
):
    """
    Ajuste de cuerpo negro para uno o varios datasets con errores en y
    + gráfico de residuales normalizados.

    El modelo usa:
        F_lambda = B_lambda(T) * Omega_sun

    donde:
        - wavelength en nm
        - irradiance en W / m^2 / nm
        - Omega_sun se calcula a partir del diámetro angular solar:
          0°32'11.41''
    """

    # ------------------------- defaults residuales -------------------------
    if residual_title_fontsize is None:
        residual_title_fontsize = axis_label_fontsize + 1
    if residual_label_fontsize is None:
        residual_label_fontsize = axis_label_fontsize
    if residual_tick_fontsize is None:
        residual_tick_fontsize = tick_fontsize
    if residual_legend_fontsize is None:
        residual_legend_fontsize = max(8, legend_fontsize - 1)

    # ------------------------- normalizar entrada --------------------------
    if not isinstance(dfs, list):
        dfs = [dfs]
    n = len(dfs)

    if labels is None:
        labels = [f'Data {i+1}' for i in range(n)]
    if len(labels) != n:
        raise ValueError("labels debe tener la misma longitud que dfs")

    if colors is None:
        colors = plt.cm.tab10(np.linspace(0, 1, n))
    if len(colors) != n:
        raise ValueError("colors debe tener la misma longitud que dfs")

    # ------------------------- helpers --------------------------
    def parse_err_for_df(df, err_spec, fallback_len):
        if err_spec is None:
            return None
        if isinstance(err_spec, str):
            return df[err_spec].to_numpy(dtype=float)
        if np.isscalar(err_spec):
            return np.full(fallback_len, float(err_spec), dtype=float)
        arr = np.asarray(err_spec, dtype=float)
        if arr.shape[0] != fallback_len:
            raise ValueError("Error array debe tener la misma longitud que el dataframe")
        return arr

    def calculate_r_squared(y_true, y_pred):
        ss_res = np.sum((y_true - y_pred) ** 2)
        ss_tot = np.sum((y_true - np.mean(y_true)) ** 2)
        return 1 - ss_res / ss_tot if ss_tot != 0 else 0.0

    # ------------------------- ángulo sólido solar --------------------------
    # diámetro angular: 0°32'11.41''
    diameter_arcsec = 32 * 60 + 11.41
    alpha_arcsec = diameter_arcsec / 2.0
    alpha_rad = alpha_arcsec / 206265.0
    Omega_sun = 2 * np.pi * (1 - np.cos(alpha_rad))

    def blackbody_model(wavelength_nm, T):
        """
        wavelength_nm en nm
        T en K

        devuelve irradiancia espectral en W / m^2 / nm
        """
        lam = wavelength_nm * 1e-9  # nm -> m

        exponent = (h * c) / (lam * k * T)
        exponent = np.clip(exponent, 1e-12, 700)

        B_lambda = (2 * h * c**2) / (lam**5) / (np.exp(exponent) - 1.0)
        # W / (m^2 sr m)

        F_lambda_m = B_lambda * Omega_sun
        # W / (m^2 m)

        F_lambda_nm = F_lambda_m * 1e-9
        # W / (m^2 nm)

        return F_lambda_nm

    # ------------------------- figura (2 filas) --------------------------
    fig, (ax_main, ax_res) = plt.subplots(
        2, 1,
        figsize=figure_size,
        dpi=figure_dpi,
        gridspec_kw={'height_ratios': [2, 1]}
    )
    fig.patch.set_facecolor('white')

    def style_axes(ax, xlabel_text, ylabel_text, title_text, t_fs, lab_fs, tick_fs):
        ax.set_xlabel(xlabel_text, fontsize=lab_fs)
        ax.set_ylabel(ylabel_text, fontsize=lab_fs)
        if title_text is not None:
            ax.set_title(title_text, fontsize=t_fs, pad=12)
        ax.tick_params(axis='both', labelsize=tick_fs, direction='in', length=5, width=1.1)
        for spine in ax.spines.values():
            spine.set_linewidth(1.1)
        ax.grid(alpha=grid_alpha, linestyle=':', linewidth=0.8)

    style_axes(
        ax_main,
        xlabel or x_col,
        ylabel or y_col,
        title,
        title_fontsize,
        axis_label_fontsize,
        tick_fontsize
    )

    # ------------------------- loop datasets --------------------------
    all_results = []
    all_residuals = []

    for df, lab, col in zip(dfs, labels, colors):
        x = df[x_col].to_numpy(dtype=float)
        y = df[y_col].to_numpy(dtype=float)
        cur_yerr = parse_err_for_df(df, yerr, len(df))

        # filtrar válidos
        mask = np.isfinite(x) & np.isfinite(y)
        if cur_yerr is not None:
            mask &= np.isfinite(cur_yerr) & (cur_yerr > 0)

        x = x[mask]
        y = y[mask]
        if cur_yerr is not None:
            cur_yerr = cur_yerr[mask]

        if len(x) < 3:
            raise ValueError(f"Dataset '{lab}': se requieren al menos 3 puntos válidos")

        # ordenar para graficar bien
        order = np.argsort(x)
        x = x[order]
        y = y[order]
        if cur_yerr is not None:
            cur_yerr = cur_yerr[order]

        # ajuste
        popt, pcov = curve_fit(
            blackbody_model,
            x, y,
            p0=[5778.0],
            sigma=cur_yerr,
            absolute_sigma=(cur_yerr is not None),
            bounds=([3000.0], [8000.0]),
            maxfev=20000
        )

        temperature = popt[0]
        temperature_err = np.sqrt(np.diag(pcov))[0]

        y_pred = blackbody_model(x, temperature)
        r2 = calculate_r_squared(y, y_pred)

        # residuales
        if cur_yerr is not None:
            residuals = (y - y_pred) / cur_yerr
        else:
            residuals = y - y_pred

        all_residuals.append((x, residuals, lab, col))

        # --- main plot ---
        ax_main.errorbar(
            x, y,
            yerr=cur_yerr,
            fmt='o',
            capsize=3,
            color=col,
            alpha=0.85,
            markersize=5,
            elinewidth=1.1,
            capthick=1.1,
            label=lab
        )

        x_line = np.linspace(np.min(x), np.max(x), 800)
        y_line = blackbody_model(x_line, temperature)
        ax_main.plot(
            x_line, y_line,
            '--',
            color=col,
            linewidth=1.8,
            label='Black body fit'
        )

        all_results.append({
            'dataset': lab,
            'temperature': temperature,
            'temperature_err': temperature_err,
            'r_squared': r2
        })

        print(f"\n--- Results for {lab} ---")
        print(f"Estimated temperature: {temperature:.6g} ± {temperature_err:.6g} K")
        print(f"R²:                    {r2:.6g}")

    ax_main.legend(
        loc=legend_loc,
        fontsize=legend_fontsize,
        frameon=True,
        fancybox=False,
        edgecolor='black'
    )

    ax_main.set_ylim(0, ax_main.get_ylim()[1])

    # ------------------------- residuales panel --------------------------
    style_axes(
        ax_res,
        xlabel or x_col,
        'Residuals',
        'Normalized Residuals' if yerr is not None else 'Residuals',
        residual_title_fontsize,
        residual_label_fontsize,
        residual_tick_fontsize
    )

    for x_data, res, lab, col in all_residuals:
        ax_res.plot(
            x_data, res,
            color=col,
            alpha=0.7,
            #s=30,
            edgecolors='none',
            label=lab
        )

    ax_res.axhline(0, color='black', linestyle='--', linewidth=1.2, alpha=0.8)

    if len(all_residuals) > 0:
        max_abs = max(np.max(np.abs(r)) for _, r, _, _ in all_residuals if len(r) > 0)
        if np.isfinite(max_abs) and max_abs > 0:
            ax_res.set_ylim(-1.1 * max_abs, 1.1 * max_abs)

    ax_res.legend(
        loc=legend_loc,
        fontsize=residual_legend_fontsize,
        frameon=True,
        fancybox=False,
        edgecolor='black'
    )

    fig.tight_layout()

    if savepath:
        fig.savefig(savepath, dpi=figure_dpi, bbox_inches='tight', facecolor='white')

    if show:
        plt.show()

    return all_results, fig


