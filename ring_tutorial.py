"""
ring_tutorial.py
================
Учебный модуль к методичке "Кольцевые резонаторы: от физики к геометрии".

Конфигурация: all-pass ring (одно кольцо + один прямой волновод).
Платформа: SOI strip 220 x 500 nm, TE-мода, рабочая длина волны ~1550 нм.

Все длины — в микрометрах (мкм), длины волн — в мкм, потери — в 1/мкм.

Запуск:  python ring_tutorial.py   -> построит все 5 рисунков в папку images/
"""

import os
import numpy as np
import matplotlib
matplotlib.use("Agg")            # рисуем в файл, без окна
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, FancyArrowPatch, Ellipse

# ----------------------------------------------------------------------
# Константы платформы (SOI strip 220 x 500 нм, TE, ~1550 нм)
# Это типовые "учебные" значения. В реальном дизайне их берут из mode solver.
# ----------------------------------------------------------------------
N_EFF      = 2.45      # эффективный индекс моды
N_G        = 4.20      # групповой индекс (учитывает дисперсию)
ALPHA_PROP = 4.6e-5    # потери распространения, 1/мкм  (~2 дБ/см)
KAPPA0     = 0.60      # макс. связь при нулевом зазоре gap=0
GAMMA      = 7.0       # скорость затухания связи с зазором, 1/мкм
LAMBDA0    = 1.55      # рабочая длина волны, мкм

# Линейная дисперсия n_eff, согласованная с n_g:
#   n_g = n_eff - lambda * dn_eff/dlambda  =>  dn_eff/dlambda = (n_eff - n_g)/lambda0
DNEFF_DLAMBDA = (N_EFF - N_G) / LAMBDA0   # 1/мкм


# ----------------------------------------------------------------------
# Уровень 1: геометрия -> базовые величины
# ----------------------------------------------------------------------
def round_trip_length(R):
    """Длина одного оборота света по кольцу: L = 2*pi*R."""
    return 2.0 * np.pi * R


def n_eff_of_lambda(lam):
    """Эффективный индекс с учётом (линейной) дисперсии."""
    return N_EFF + DNEFF_DLAMBDA * (lam - LAMBDA0)


def amplitude_loss(R, alpha=ALPHA_PROP):
    """Множитель A = exp(-alpha*L/2): доля АМПЛИТУДЫ, дожившая до конца оборота."""
    return np.exp(-alpha * round_trip_length(R) / 2.0)


def kappa(gap, kappa0=KAPPA0, gamma=GAMMA):
    """Cross-coupling: какая доля амплитуды перепрыгивает в кольцо. kappa = k0*exp(-gamma*gap)."""
    return kappa0 * np.exp(-gamma * gap)


def self_coupling(gap):
    """Self-coupling t: что прошло мимо. Из закона сохранения t^2 + kappa^2 = 1."""
    return np.sqrt(1.0 - kappa(gap) ** 2)


# ----------------------------------------------------------------------
# Уровень 2: спектр пропускания
# ----------------------------------------------------------------------
def transmission(lam, R, gap):
    """
    Пропускание all-pass кольца T(lambda) = |b_out / b_in|^2.
        T = (t^2 - 2 t A cos(phi) + A^2) / (1 - 2 t A cos(phi) + (t A)^2)
    где phi = 2*pi*n_eff*L/lambda — фаза за один оборот.
    """
    L = round_trip_length(R)
    A = amplitude_loss(R)
    t = self_coupling(gap)
    phi = 2.0 * np.pi * n_eff_of_lambda(lam) * L / lam
    num = t ** 2 - 2.0 * t * A * np.cos(phi) + A ** 2
    den = 1.0 - 2.0 * t * A * np.cos(phi) + (t * A) ** 2
    return num / den


# ----------------------------------------------------------------------
# Уровень 3: метрики (FOMs)
# ----------------------------------------------------------------------
def FSR(R, lam=LAMBDA0):
    """Free Spectral Range — расстояние между резонансами: FSR = lambda^2 / (n_g * L)."""
    return lam ** 2 / (N_G * round_trip_length(R))


def Q_intrinsic(R, lam=LAMBDA0, alpha=ALPHA_PROP):
    """Внутренняя добротность (только потери в кольце)."""
    return 2.0 * np.pi * N_G / (lam * alpha)


def Q_coupling(R, gap, lam=LAMBDA0):
    """Добротность, связанная с утечкой через coupler."""
    L = round_trip_length(R)
    return np.pi * N_G * L / (lam * kappa(gap) ** 2)


def Q_loaded(R, gap, lam=LAMBDA0):
    """Полная (наблюдаемая) добротность: 1/Q = 1/Q_int + 1/Q_coupling."""
    return 1.0 / (1.0 / Q_intrinsic(R, lam) + 1.0 / Q_coupling(R, gap, lam))


def extinction_ratio_dB(R, gap):
    """Глубина провала в дБ: ER = -10*log10(T_min / T_max)."""
    A = amplitude_loss(R)
    t = self_coupling(gap)
    T_min = ((t - A) / (1.0 - t * A)) ** 2
    T_max = ((t + A) / (1.0 + t * A)) ** 2
    return -10.0 * np.log10(T_min / T_max)


def critical_gap(R, kappa0=KAPPA0, gamma=GAMMA):
    """Зазор критической связи (t = A): максимальный ER, провал до нуля."""
    k_crit = np.sqrt(1.0 - amplitude_loss(R) ** 2)
    return -np.log(k_crit / kappa0) / gamma


# ----------------------------------------------------------------------
# Уровень 5: обратная задача (inverse design)
# ----------------------------------------------------------------------
def design_for_FSR_Q(FSR_target, Q_target, lam=LAMBDA0):
    """
    Подобрать радиус R и зазор gap под заданные FSR и нагруженную добротность Q.
    Возвращает словарь с параметрами и пометкой о выполнимости.
    """
    R = lam ** 2 / (2.0 * np.pi * N_G * FSR_target)     # шаг 1: радиус из FSR
    Qi = Q_intrinsic(R, lam)                            # шаг 2: потолок добротности
    if Q_target >= Qi:
        return {"feasible": False, "R": R, "Q_intrinsic": Qi,
                "reason": "Q_target >= Q_intrinsic: потери в кольце не дадут такую добротность"}
    Qc = Qi * Q_target / (Qi - Q_target)                # шаг 3: нужная Q_coupling
    L = round_trip_length(R)
    k = np.sqrt(np.pi * N_G * L / (lam * Qc))           # -> kappa
    gap = -np.log(k / KAPPA0) / GAMMA                   # шаг 4: kappa -> gap
    return {"feasible": gap >= 0.05, "R": R, "gap": gap, "kappa": k,
            "Q_intrinsic": Qi, "Q_coupling": Qc,
            "fab_ok": gap >= 0.05}  # технологический предел ~50 нм


# ======================================================================
#                         Р И С У Н К И
# ======================================================================
plt.rcParams.update({"font.size": 12, "figure.dpi": 130,
                     "axes.grid": True, "grid.alpha": 0.3})

COL_BUS  = "#2b6cb0"
COL_RING = "#c53030"
COL_OK   = "#2f855a"


def fig1_schematic(path="images/fig1_schematic.png"):
    """Рисунок к уровню 1: устройство all-pass кольцевого резонатора."""
    fig, ax = plt.subplots(figsize=(7.4, 4.8))
    ax.set_aspect("equal"); ax.axis("off")

    R = 1.0
    cy = 0.65                       # центр кольца
    gap = 0.34                      # видимый зазор
    ring_bottom = cy - R
    bus_y = ring_bottom - gap       # волновод ниже кольца с явным зазором

    # кольцо
    ax.add_patch(Circle((0, cy), R, fill=False, lw=9, color=COL_RING))
    # прямой волновод (bus)
    ax.plot([-2.2, 2.2], [bus_y, bus_y], lw=9, color=COL_BUS,
            solid_capstyle="round")

    # входной и выходной пучки
    ax.annotate("", xy=(-2.18, bus_y), xytext=(-2.85, bus_y),
                arrowprops=dict(arrowstyle="-|>", color="black", lw=1.6))
    ax.annotate("", xy=(2.85, bus_y), xytext=(2.18, bus_y),
                arrowprops=dict(arrowstyle="-|>", color="black", lw=1.6))
    ax.text(-2.9, bus_y + 0.16, r"$b_{in}$", ha="right", fontsize=15)
    ax.text(2.9, bus_y + 0.16, r"$b_{out}$", ha="left", fontsize=15)

    # область связи (coupler)
    coupler_y = (bus_y + ring_bottom) / 2.0
    ax.add_patch(Ellipse((0, coupler_y), 1.1, 0.62, fill=False,
                 lw=1.6, ls="--", color="gray"))
    ax.text(0.0, coupler_y - 0.42, r"coupler:  $t,\ \kappa$",
            ha="center", fontsize=13, color="gray")

    # радиус
    ax.annotate("", xy=(R*np.cos(np.pi/4), cy + R*np.sin(np.pi/4)), xytext=(0, cy),
                arrowprops=dict(arrowstyle="-|>", color="black", lw=1.4))
    ax.text(0.30, cy + 0.42, r"$R$", fontsize=15)

    # зазор: явная вертикальная двойная стрелка справа
    ax.annotate("", xy=(1.35, ring_bottom), xytext=(1.35, bus_y),
                arrowprops=dict(arrowstyle="<->", color="black", lw=1.4))
    ax.text(1.45, coupler_y, r"$g$", fontsize=15, va="center")

    # стрелка "свет бежит по кольцу"
    ax.annotate("", xy=(-R*np.cos(0.25), cy + R*np.sin(0.25)),
                xytext=(-R*np.cos(0.75), cy + R*np.sin(0.75)),
                arrowprops=dict(arrowstyle="-|>", color=COL_RING, lw=2,
                                connectionstyle="arc3,rad=0.3"))

    ax.set_xlim(-3.3, 3.3); ax.set_ylim(-0.9, 1.95)
    ax.set_title("Уровень 1. All-pass кольцевой резонатор", fontsize=13)
    fig.tight_layout(); fig.savefig(path, bbox_inches="tight"); plt.close(fig)


def fig2_transmission(path="images/fig2_transmission.png"):
    """Рисунок к уровню 2: спектр пропускания, резонансы и FSR."""
    R, gap = 10.0, 0.34
    lam = np.linspace(1.530, 1.570, 60000)
    T = transmission(lam, R, gap)

    fig, ax = plt.subplots(figsize=(7.6, 4.4))
    ax.plot(lam * 1000, T, lw=1.2, color=COL_BUS)
    ax.set_xlabel("Длина волны, нм")
    ax.set_ylabel(r"Пропускание $T$")
    ax.set_ylim(-0.03, 1.05)

    # найдём положения двух соседних резонансов и покажем FSR
    idx = np.where((T[1:-1] < T[:-2]) & (T[1:-1] < T[2:]))[0] + 1
    idx = idx[T[idx] < 0.5]
    if len(idx) >= 2:
        l1, l2 = lam[idx[1]] * 1000, lam[idx[2]] * 1000
        ax.annotate("", xy=(l2, 1.0), xytext=(l1, 1.0),
                    arrowprops=dict(arrowstyle="<->", color=COL_RING, lw=1.5))
        ax.text((l1 + l2) / 2, 1.02, f"FSR ≈ {l2 - l1:.1f} нм",
                ha="center", color=COL_RING, fontsize=12)

    ax.set_title(f"Уровень 2. Спектр пропускания (R={R:.0f} мкм, g={gap*1000:.0f} нм)",
                 fontsize=12)
    fig.tight_layout(); fig.savefig(path, bbox_inches="tight"); plt.close(fig)


def nearest_resonance(R, gap, target=1.553):
    """Найти длину волны резонанса (минимум T) рядом с target."""
    lam = np.linspace(target - 0.006, target + 0.006, 400000)
    return lam[np.argmin(transmission(lam, R, gap))]


def fig3_metrics(path="images/fig3_metrics.png"):
    """Рисунок к уровню 3: Q и ER на одном резонансе + три режима связи."""
    R = 10.0
    # положение резонанса не зависит от связи -> найдём один раз
    lam_res = nearest_resonance(R, 0.30)
    fig, (axL, axR) = plt.subplots(1, 2, figsize=(11.5, 4.4))

    # --- левая панель: один резонанс, FWHM (Q) и ER ---
    g = 0.33                              # чуть мимо критики -> глубокий, но КОНЕЧНЫЙ ER
    lam = np.linspace(lam_res - 6e-5, lam_res + 6e-5, 60000)   # окно ±0.06 нм
    T = transmission(lam, R, g)
    axL.plot((lam - lam_res) * 1000, T, lw=1.6, color=COL_BUS)
    Tmin = T.min()
    half = (1.0 + Tmin) / 2.0
    above = (lam[T <= half] - lam_res) * 1000
    if len(above) > 1:
        axL.hlines(half, above[0], above[-1], color=COL_RING, lw=1.6)
        axL.annotate(r"FWHM $\Delta\lambda \Rightarrow Q=\lambda/\Delta\lambda$",
                     (above[-1], half), (above[-1] + 0.004, half + 0.06),
                     fontsize=11, color=COL_RING, va="center")
    axL.annotate("", xy=(0, 1.0), xytext=(0, Tmin),
                 arrowprops=dict(arrowstyle="<->", color=COL_OK, lw=1.5))
    axL.text(0.004, 0.5, f"ER ≈ {extinction_ratio_dB(R, g):.0f} дБ",
             color=COL_OK, fontsize=12)
    axL.set_xlabel(f"Отстройка от {lam_res*1000:.1f} нм, нм")
    axL.set_ylabel(r"$T$"); axL.set_ylim(-0.03, 1.08)
    axL.set_title(f"Один резонанс: Q$_L$ ≈ {Q_loaded(R, g):,.0f}".replace(",", " "),
                  fontsize=11)

    # --- правая панель: три режима связи (центрированы на резонансе) ---
    gc = critical_gap(R)
    regimes = [(0.20, "overcoupled (t<A)",  "#dd6b20"),
               (gc,   "critical (t=A)",     COL_OK),
               (0.50, "undercoupled (t>A)", "#805ad5")]
    lam2 = np.linspace(lam_res - 2e-4, lam_res + 2e-4, 120000)   # окно ±0.2 нм
    for g_, label, c in regimes:
        axR.plot((lam2 - lam_res) * 1000, transmission(lam2, R, g_), lw=1.6,
                 color=c, label=f"{label}, g={g_*1000:.0f} нм")
    axR.set_xlabel(f"Отстройка от {lam_res*1000:.1f} нм, нм")
    axR.set_ylabel(r"$T$"); axR.set_ylim(-0.03, 1.08)
    axR.legend(fontsize=9, loc="lower left")
    axR.set_title("Три режима связи", fontsize=11)

    fig.suptitle("Уровень 3. Что мы измеряем: Q, ER и режим связи", fontsize=13)
    fig.tight_layout(); fig.savefig(path, bbox_inches="tight"); plt.close(fig)


def fig4_geometry(path="images/fig4_geometry.png"):
    """Рисунок к уровню 4: геометрия -> физика. alpha(R) и kappa(gap)."""
    fig, (axL, axR) = plt.subplots(1, 2, figsize=(11.5, 4.4))

    # alpha(R): потери распространения + изгибные потери (иллюстративно)
    R = np.linspace(2.0, 20.0, 400)
    R_crit, R0, a0 = 4.0, 1.0, 5e-3          # иллюстративные параметры изгиба
    alpha_bend = a0 * np.exp(-(R - R_crit) / R0)
    alpha_total = ALPHA_PROP + alpha_bend
    axL.semilogy(R, alpha_total * 1e4 * 10 / np.log(10), lw=2, color=COL_RING,
                 label=r"$\alpha(R)$ полные")  # перевод в дБ/см для наглядности
    axL.axhline(ALPHA_PROP * 1e4 * 10 / np.log(10), ls="--", color="gray",
                label=r"$\alpha_{prop}\approx2$ дБ/см")
    axL.set_xlabel("Радиус R, мкм"); axL.set_ylabel("Потери, дБ/см")
    axL.legend(fontsize=10); axL.set_title("Малый радиус → изгибные потери", fontsize=11)

    # kappa(gap)
    g = np.linspace(0.05, 0.6, 400)
    axR.plot(g * 1000, kappa(g), lw=2, color=COL_BUS, label=r"$\kappa(g)=\kappa_0 e^{-\gamma g}$")
    gc = critical_gap(10.0)
    axR.axvline(gc * 1000, ls="--", color=COL_OK,
                label=f"крит. связь @R=10: g≈{gc*1000:.0f} нм")
    axR.set_xlabel("Зазор g, нм"); axR.set_ylabel(r"$\kappa$")
    axR.legend(fontsize=10); axR.set_title("Больше зазор → слабее связь", fontsize=11)

    fig.suptitle("Уровень 4. От геометрии (R, g) к физике (α, κ)", fontsize=13)
    fig.tight_layout(); fig.savefig(path, bbox_inches="tight"); plt.close(fig)


def fig5_designmap(path="images/fig5_designmap.png"):
    """Рисунок к уровню 5: карта дизайна — линии FSR и Q в плоскости (R, gap)."""
    R = np.linspace(5.0, 25.0, 220)
    g = np.linspace(0.10, 0.45, 220)
    RR, GG = np.meshgrid(R, g)

    FSR_map = FSR(RR) * 1000.0                       # нм
    Q_map = Q_loaded(RR, GG)

    fig, ax = plt.subplots(figsize=(7.8, 5.2))
    # Q как фон
    pc = ax.pcolormesh(RR, GG * 1000, np.log10(Q_map), shading="auto", cmap="viridis")
    cb = fig.colorbar(pc, ax=ax); cb.set_label(r"$\log_{10} Q_{loaded}$")

    # линии постоянного FSR (зависит только от R)
    cs1 = ax.contour(RR, GG * 1000, FSR_map, levels=[4, 6, 8, 10, 12],
                     colors="white", linewidths=1.2)
    ax.clabel(cs1, fmt="FSR=%.0f нм", fontsize=9)

    # линия критической связи
    g_crit = np.array([critical_gap(r) for r in R]) * 1000
    ax.plot(R, g_crit, color=COL_RING, lw=2.5, label="критическая связь (max ER)")

    # пример точки дизайна
    d = design_for_FSR_Q(FSR_target=0.006, Q_target=2.0e4)
    ax.plot(d["R"], d["gap"] * 1000, "*", ms=18, color="yellow",
            markeredgecolor="black",
            label=f"пример: FSR=6нм, Q=20k → R={d['R']:.1f}мкм, g={d['gap']*1000:.0f}нм")

    ax.set_xlabel("Радиус R, мкм"); ax.set_ylabel("Зазор g, нм")
    ax.legend(fontsize=9, loc="upper right")
    ax.set_title("Уровень 5. Карта дизайна: выбираем R и g под FSR и Q", fontsize=12)
    fig.tight_layout(); fig.savefig(path, bbox_inches="tight"); plt.close(fig)


if __name__ == "__main__":
    # работаем в папке самого скрипта и создаём images/, если её нет
    os.chdir(os.path.dirname(os.path.abspath(__file__)))
    os.makedirs("images", exist_ok=True)
    fig1_schematic()
    fig2_transmission()
    fig3_metrics()
    fig4_geometry()
    fig5_designmap()
    print("OK: 5 рисунков сохранены в images/")
    # короткая сводка чисел для проверки
    print(f"FSR(R=10)      = {FSR(10)*1000:.2f} нм")
    print(f"Q_intrinsic    = {Q_intrinsic(10):.0f}")
    print(f"crit gap(R=10) = {critical_gap(10)*1000:.0f} нм")
    print(f"ER @crit       = {extinction_ratio_dB(10, critical_gap(10)):.1f} дБ")
    print("inverse:", design_for_FSR_Q(0.006, 2.0e4))
