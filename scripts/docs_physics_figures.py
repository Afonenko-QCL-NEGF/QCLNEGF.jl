#!/usr/bin/env python3
"""Render the scientific atlas from Julia-exported CSV; no solver duplication.

Run examples/10_physics_atlas.jl first. Analytic teaching illustrations are
explicitly labelled in every corresponding figure and in the atlas chapter.
Requires Python 3, numpy and matplotlib; no HDF5 or server environment.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
COLORS = ["#22577a", "#c35c2c", "#427b58", "#92649a", "#c0942e"]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data", type=Path, default=ROOT / "docs/src/assets/physics/data")
    parser.add_argument("--output", type=Path, default=ROOT / "docs/src/assets/physics")
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10,
        "axes.spines.top": False, "axes.spines.right": False,
        "axes.prop_cycle": matplotlib.cycler(color=COLORS), "axes.titlepad": 10,
        "figure.dpi": 110, "savefig.dpi": 180, "svg.fonttype": "none", "svg.hashsalt": "qcl-negf-physics-atlas"})

    def load(name):
        return np.genfromtxt(args.data / name, delimiter=",", names=True,
            dtype=None, encoding="utf-8")

    def save(fig, name):
        fig.savefig(args.output / (name + ".svg"), bbox_inches="tight", metadata={"Date": None})
        plt.close(fig)

    structure = load("structure.csv")
    z, ec = structure["z_nm"], structure["Ec_meV"]
    layer_data = load("layers.csv")
    parameters = {row["name"]:float(row["value"]) for row in load("parameters.csv")}
    period = parameters["period"]
    def layers(ax):
        for layer in layer_data:
            if layer["material"]=="AlGaAs":
                ax.axvspan(layer["left_nm"],layer["right_nm"],color="#22577a",alpha=.08,lw=0)
            if str(layer["doped"]).lower()=="true":
                ax.axvspan(layer["left_nm"],layer["right_nm"],color="#c0942e",alpha=.13,lw=0)
        ax.set_xlim(0,period)
        ax.set_xlabel(r"Координата $z$, нм")

    fig, axes = plt.subplots(3,1,figsize=(10,8),sharex=True,layout="constrained")
    axes[0].step(z,ec,where="mid",label="Край зоны без поля")
    axes[0].step(z,ec+structure["bias_meV"],where="mid",label="С заданным полем")
    axes[0].set(ylabel="Энергия, мэВ", title="Реальная геометрия reference design: один период 29,61 нм")
    axes[0].legend(loc="upper right",fontsize=9)
    axes[1].step(z,structure["mass_ratio"],where="mid",color=COLORS[2])
    axes[1].set_ylabel(r"Эффективная масса $m_z/m_0$")
    axes[2].step(z,structure["donors_m3"]/1e23,where="mid",color=COLORS[4])
    axes[2].set_ylabel(r"Доноры, $10^{23}$ м$^{-3}$")
    for ax in axes: layers(ax)
    save(fig,"reference2019_structure")

    metrics, waves = load("basis_metrics.csv"), load("basis_envelopes.csv")
    methods = [("energy","Энергетические состояния: без поля"),
               ("pzp_tails",r"Локализация $PzP$: без поля"),
               ("bloch_wannier","Функции Ванье: периодическая граница"),
               ("wannier_stark","Состояния Ванье–Штарка: заданное поле")]
    fig, axes = plt.subplots(4,1,figsize=(11,12),sharex=True,layout="constrained")
    for ax,(method,title) in zip(axes,methods):
        for number,state in enumerate(range(11,16)):
            row = waves[(waves["method"]==method)&(waves["state"]==state)]
            probability=row["probability_per_nm"]
            y=number+0.82*probability/max(probability)
            ax.plot(row["z_nm"],y,color=COLORS[number],lw=1.35)
            ax.fill_between(row["z_nm"],number,y,color=COLORS[number],alpha=.13)
        for pos in np.arange(-2,4)*period: ax.axvline(pos,color=".75",lw=.7)
        ax.axvspan(0,period,color="#c0942e",alpha=.1)
        ax.set(title=title,yticks=range(5),yticklabels=range(11,16),ylabel="Номер состояния")
    axes[-1].set(xlabel=r"Координата $z$, нм",xlim=(-2*period,3*period))
    fig.suptitle("Плотности вероятности в окне пяти периодов\nВысота каждой кривой нормирована отдельно; сдвиг по вертикали не является энергией",fontsize=12)
    save(fig,"basis_localization")

    operators=load("basis_operators.csv")
    fig, axes=plt.subplots(2,3,figsize=(13,8),layout="constrained")
    for i,method in enumerate(["energy","pzp_tails"]):
        row=operators[operators["method"]==method]
        arrays=[np.hypot(row["H_real_meV"],row["H_imag_meV"]),
                np.hypot(row["Z_real_nm"],row["Z_imag_nm"]),row["speed_nm_ps"]]
        for ax,a,title in zip(axes[i],arrays,[r"$|H_{ab}|$, мэВ",r"$|z_{ab}|$, нм",r"$|v_{ab}|$, нм/пс"]):
            values=a.reshape(25,25)
            im=ax.imshow(values,origin="lower",cmap="viridis",extent=[.5,25.5,.5,25.5])
            ax.set(title=title,xlabel="Индекс b",ylabel="Индекс a")
            fig.colorbar(im,ax=ax,shrink=.8)
        axes[i,0].text(-.28,.5,"Энергетический базис" if i==0 else "Локализованный базис",
            transform=axes[i,0].transAxes,rotation=90,ha="center",va="center")
    fig.suptitle("Один оператор в двух базисах: диагональность не сохраняется при локализации")
    save(fig,"basis_operators")

    row=operators[operators["method"]=="pzp_tails"]
    h=(row["H_real_meV"]+1j*row["H_imag_meV"]).reshape(25,25)
    position=(row["Z_real_nm"]+1j*row["Z_imag_nm"]).reshape(25,25)
    velocity=1j*(h@position-position@h)/parameters["reduced_planck_constant"]
    phase=np.linspace(-np.pi,np.pi,201)
    coherence=np.linspace(0,.5,101)
    rho_ba=coherence[:,None]*np.exp(-1j*phase[None,:])
    average_velocity=2*np.real(velocity[10,11]*rho_ba)
    fig,ax=plt.subplots(figsize=(9,5),layout="constrained")
    limit=np.max(np.abs(average_velocity))
    im=ax.pcolormesh(phase/np.pi,coherence,average_velocity,shading="auto",cmap="RdBu_r",
        vmin=-limit,vmax=limit,rasterized=True)
    ax.set(xlabel="Фаза ρ₁₂ / π",ylabel="Модуль когерентности |ρ₁₂|",
        title="Одинаковые заселения ρ₁₁ = ρ₂₂ = 1/2; различный перенос")
    fig.colorbar(im,ax=ax,label="Средняя скорость, нм/пс")
    fig.suptitle("Два реальных состояния PzP и заданная матрица плотности; это не ток SCBA")
    save(fig,"coherence_current")

    fig,axes=plt.subplots(1,2,figsize=(11,4.5),layout="constrained")
    for method,title in methods[:3]:
        row=metrics[(metrics["method"]==method)&(metrics["state"]>=11)&(metrics["state"]<=15)]
        axes[0].plot(np.arange(1,6),row["spread_nm"],"o-",label=title.split(":")[0])
        axes[1].plot(np.arange(1,6),row["outside_central_weight"],"o-")
    axes[0].set(ylabel="Среднеквадратичная ширина, нм",xlabel="Состояние в выбранной группе")
    axes[1].set(ylabel="Вес за пределами центрального периода",xlabel="Состояние в выбранной группе",ylim=(0,1))
    axes[0].legend(fontsize=8)
    fig.suptitle("Ширина и хвосты: характеристики представления, а не ошибка транспорта")
    save(fig,"localization_metrics")

    spectrum=load("spectral_population.csv")
    E=np.unique(spectrum["energy_meV"])
    zz=np.unique(spectrum["z_nm"])
    fig,axes=plt.subplots(1,2,figsize=(12,5),layout="constrained",sharey=True)
    maximum=spectrum["ldos_per_meV_nm"].max()
    for ax,key,title in zip(axes,["ldos_per_meV_nm","occupied_per_meV_nm"],
            [r"Доступные состояния: $A(z,E)/(2\pi)$",r"Занятые состояния: $f(E)A(z,E)/(2\pi)$"]):
        im=ax.pcolormesh(zz,E,spectrum[key].reshape(len(E),len(zz)),shading="auto",
            cmap="magma",norm=LogNorm(vmin=maximum*1e-4,vmax=maximum),rasterized=True)
        ax.step(z,ec,where="mid",color="#90d9e1",lw=1,alpha=.7)
        ax.set(title=title,xlabel="Координата z, нм",ylim=(E[0],E[-1]))
    axes[0].set_ylabel("Энергия, мэВ")
    fig.colorbar(im,ax=axes,label=r"Плотность, (мэВ·нм)$^{-1}$",shrink=.8)
    fig.suptitle("Конечное окно без рассеяния: η = 2 мэВ; заданные μ = 60 мэВ и T = 200 K")
    save(fig,"spectral_population")

    hartree=load("hartree.csv")
    fig,axes=plt.subplots(3,1,figsize=(10,8),sharex=True,layout="constrained")
    axes[0].step(z,hartree["donors_m3"]/1e23,where="mid",label="Ионизованные доноры")
    axes[0].plot(z,hartree["electrons_m3"]/1e23,label="Пробная электронная плотность")
    axes[0].set_ylabel(r"Плотность, $10^{23}$ м$^{-3}$")
    axes[0].legend(fontsize=9)
    axes[1].plot(z,hartree["hartree_meV"],color=COLORS[2])
    axes[1].axhline(0,color=".6",lw=.7)
    axes[1].set_ylabel(r"Энергия Хартри $U_H$, мэВ")
    axes[2].step(z,ec+structure["bias_meV"],where="mid",label="Без Хартри")
    axes[2].step(z,hartree["total_band_meV"],where="mid",label="Один отклик Пуассона")
    axes[2].set_ylabel("Край зоны, мэВ")
    axes[2].legend(fontsize=9)
    for ax in axes: layers(ax)
    fig.suptitle("Периодическое уравнение Пуассона: нейтральный пробный заряд, среднее U_H = 0")
    save(fig,"hartree_response")

    ff=load("form_factors.csv")
    fig,axes=plt.subplots(1,3,figsize=(13,4.2),layout="constrained")
    axes[0].plot(ff["q_per_nm"],ff["diagonal_F2"],label=r"$|F_{aa}|^2$")
    axes[0].plot(ff["q_per_nm"],ff["transition_F2"],label=r"$|F_{ab}|^2$, a ≠ b")
    axes[0].set(title="Реальные функции PzP",ylabel="Квадрат форм-фактора",xlabel=r"$q_z$, нм$^{-1}$")
    axes[0].legend(fontsize=9)
    for key,label in [("normalized_screening",r"$q_s^2/(q^2+q_s^2)$"),("normalized_gaussian_IFR",r"$\exp(-q^2\Lambda^2/4)$")]:
        axes[1].plot(ff["q_per_nm"],ff[key],label=label)
    axes[1].set(title="Отдельные факторы, не полные ядра",ylabel="Нормированный фактор",xlabel=r"$q$, нм$^{-1}$")
    axes[1].legend(fontsize=9)
    temp=np.linspace(20,350,250)
    bose=1/np.expm1(parameters["phonon_energy"]/(parameters["boltzmann_constant"]*temp))
    axes[2].plot(temp,bose,label=r"Поглощение: $N_{LO}$")
    axes[2].plot(temp,bose+1,label=r"Испускание: $N_{LO}+1$")
    axes[2].set(title="Бозе-факторы: ℏω_LO = 36,7 мэВ",xlabel="Температура фононов, K",ylabel="Статистический множитель")
    axes[2].legend(fontsize=9)
    save(fig,"scattering_factors")

    disp=load("dispersion.csv")
    fig,ax=plt.subplots(figsize=(9,5),layout="constrained")
    for state in range(1,6):
        row=disp[disp["state"]==state]
        ax.plot(row["k_per_nm"],row["energy_meV"],label=f"Подзона {state}")
    ax.set(title="Поперечное движение: BDD-оператор одной ячейки с переменной массой",
        xlabel=r"Модуль поперечного волнового вектора, нм$^{-1}$",ylabel="Энергия, мэВ")
    ax.legend()
    save(fig,"band_dispersion")

    optical=load("optical_parameters.csv")
    transition=float(optical["transition_meV"])
    gamma=float(optical["halfwidth_meV"])
    photon=np.linspace(max(0,transition-20),transition+20,400)
    lineshape=gamma**2/((photon-transition)**2+gamma**2)
    fig,axes=plt.subplots(1,2,figsize=(11,4.5),layout="constrained")
    for inversion in [-.3,0,.3]:
        axes[0].plot(photon,inversion*lineshape,label=f"n₂ − n₁ = {inversion:+.1f}")
    axes[0].axhline(0,color=".6",lw=.7)
    axes[0].set(title="Двухуровневая иллюстрация",xlabel="Энергия фотона, мэВ",ylabel="Усиление, условные единицы")
    axes[0].legend(fontsize=9)
    axes[1].hlines([0,transition],.25,.75,color=[COLORS[0],COLORS[1]],lw=2)
    axes[1].annotate("",xy=(.5,0),xytext=(.5,transition),arrowprops={"arrowstyle":"->","color":COLORS[2],"lw":2})
    axes[1].text(.8,transition/2,f"ΔE = {transition:.2f} мэВ\n|z₁₂| = {float(optical['dipole_nm']):.2f} нм",va="center")
    axes[1].set(xlim=(0,1.5),xticks=[],ylabel="Энергия относительно уровня 1, мэВ",title="Параметры из реального H одной ячейки")
    fig.suptitle("Заданные заселения и ширина линии; это не расчёт усиления работающего reference design")
    save(fig,"optical_response")

    fig,axes=plt.subplots(1,2,figsize=(11,4.5),layout="constrained")
    energy=np.linspace(-5,5,501)
    sigma=1/(energy+1j*.7)
    axes[0].plot(energy,sigma.real,label=r"$\mathrm{Re}\,\Sigma^R$")
    axes[0].plot(energy,-2*sigma.imag,label=r"$\Gamma=-2\mathrm{Im}\,\Sigma^R$")
    axes[0].axhline(0,color=".6",lw=.7)
    axes[0].set(xlabel="Энергия относительно резонанса, условные единицы",ylabel="Собственная энергия",title="Причинный полюс: сдвиг и уширение")
    axes[0].legend()
    for alpha in [.2,.6,1.0]:
        value=0j
        residuals=[]
        for _ in range(220):
            candidate=.8**2/(.2+.06j-value)
            residuals.append(abs(candidate-value)/.8)
            value=(1-alpha)*value+alpha*candidate
        axes[1].semilogy(range(1,221),np.maximum(residuals,1e-16),label=f"α = {alpha}")
    axes[1].axhline(1e-8,color=".5",ls=":",lw=1)
    axes[1].set(xlabel="Номер итерации",ylabel="Невязка неподвижной точки",title="Скалярная модель SCBA: влияние смешивания")
    axes[1].legend()
    fig.suptitle("Аналитические учебные модели; параметры безразмерны")
    save(fig,"scba_fixed_point")
    print(f"Rendered 11 figures (SVG) in {args.output}")


if __name__ == "__main__":
    main()
