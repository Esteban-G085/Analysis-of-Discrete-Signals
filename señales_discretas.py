"""Señales discretas: muestreo, Nyquist-Shannon, aliasing y FFT.

Dos figuras de 4x2, una fila por señal: a la izquierda la curva continua con las
muestras superpuestas, a la derecha el espectro. La primera muestrea a 100 Hz
(cumple Nyquist-Shannon) y la segunda a 25 Hz (aliasing intencional). Al mover el
cursor por cualquier panel se leen t, x(t), A, F y, en frecuencia, |X(f)|.
"""

from math import pi, radians

import matplotlib.pyplot as plt
import numpy as np

T = 2.0          # duración de la observación
T_PANEL = 1.0    # los paneles de tiempo muestran solo el primer segundo
t_cont = np.linspace(0, T, 20000)

FS_NYQUIST = 100.0   # cumple fs >= 2 · f_max
FS_ALIASING = 25.0   # no cumple: su límite de Nyquist es 12.5 Hz


def senoidal(t, A, F, fase=0.0):
    return A * np.sin(2 * pi * F * t + radians(fase))


def f_alias(f, fs):
    """Frecuencia a la que se pliega una componente f al muestrearla a fs."""
    return abs(f - fs * round(f / fs))


def muestrear(x_c, fs):
    """(t_m, x_m, f, |X(f)|) de la curva continua x_c muestreada a fs."""
    N = int(fs * T)
    t_m = np.arange(N) / fs
    x_m = np.interp(t_m, t_cont, x_c)
    X = np.fft.fft(x_m)
    f = np.fft.fftfreq(N, d=1 / fs)[: N // 2]
    return t_m, x_m, f, 2.0 / N * np.abs(X[: N // 2])


# "componentes" son las frecuencias y amplitudes que la FFT debe encontrar; con
# fs baja cada una se pliega a f_alias y es lo que se ve en el espectro.
SENALES = [
    {
        "nombre": "S1: Senoidal simple",
        "color": "blue",
        "x": senoidal(t_cont, 1.5, 10.0),
        "componentes": [(10.0, 1.5)],
    },
    {
        "nombre": "S2: Multitono",
        "color": "green",
        "x": senoidal(t_cont, 1.0, 5.0) + 0.8 * senoidal(t_cont, 1.0, 20.0, 90.0),
        "componentes": [(5.0, 1.0), (20.0, 0.8)],
    },
    {
        "nombre": "S3: Modulada en AM",
        "color": "orange",
        "x": (1.0 + 0.5 * senoidal(t_cont, 1.0, 2.0, 90.0)) * senoidal(t_cont, 1.0, 15.0),
        "componentes": [(13.0, 0.25), (15.0, 1.0), (17.0, 0.25)],
    },
    {
        "nombre": "S4: Decaimiento exponencial",
        "color": "purple",
        "x": 2.0 * np.exp(-1.2 * t_cont) * senoidal(t_cont, 1.0, 30.0),
        "componentes": [(30.0, 2.0)],
    },
]

F_MAX = max(f for s in SENALES for f, _ in s["componentes"])


def frecuencias_observadas(senal, fs):
    """"10 Hz" o "20 Hz -> 5 Hz": qué frecuencias aparecerán en el espectro."""
    partes = []
    for f, _ in senal["componentes"]:
        fa = f_alias(f, fs)
        partes.append(f"{f:g} Hz" if np.isclose(fa, f) else f"{f:g} Hz -> {fa:g} Hz")
    return ", ".join(partes)


def resumen():
    print(f"\nf_max del conjunto = {F_MAX:g} Hz  ->  Nyquist exige fs >= {2 * F_MAX:g} Hz")
    for s in SENALES:
        f_max = max(f for f, _ in s["componentes"])
        print(f'\n  {s["nombre"]}   f_max = {f_max:g} Hz')
        for nombre, fs in (("cumple", FS_NYQUIST), ("ALIASING", FS_ALIASING)):
            print(f"      fs = {fs:>5.1f} Hz   Nyquist = {fs / 2:>5.1f} Hz   {nombre:8}"
                  f"  ->  {frecuencias_observadas(s, fs)}")
    print("\nCursor sobre cualquier panel: t, x(t), A, F | f, |X(f)|, F\n")


def dibujar(fs, titulo, color_muestras):
    fig, ejes = plt.subplots(4, 2, figsize=(14, 10))
    fig.suptitle(titulo, fontsize=14, weight="bold")
    lecturas = {}   # eje -> (dominio, señal, eje_x, eje_y, frecuencias observadas)
    marcas = {}     # eje -> (punto, etiqueta), se crean al primer hover

    for i, s in enumerate(SENALES):
        t_m, x_m, f, mag = muestrear(s["x"], fs)
        observadas = frecuencias_observadas(s, fs)
        lecturas[ejes[i, 0]] = ("tiempo", s, t_m, x_m, observadas)
        lecturas[ejes[i, 1]] = ("frecuencia", s, f, mag, observadas)

        ejes[i, 0].plot(t_cont, s["x"], color=s["color"], alpha=0.35, label="Continua")
        ejes[i, 0].plot(t_m, x_m, color=color_muestras, marker=".", ls="-", label=f"Muestras (fs={fs:g} Hz)")
        ejes[i, 0].set(title=f'{s["nombre"]} — tiempo', ylabel="Amplitud", xlim=(0, T_PANEL))
        ejes[i, 0].grid(True, alpha=0.3)
        ejes[i, 0].legend(fontsize=8, loc="upper right")

        espectro = ejes[i, 1].stem(f, mag)
        for artista in (espectro.markerline, espectro.stemlines):
            artista.set_color(color_muestras)
        ejes[i, 1].axvline(fs / 2, color="black", ls="--", lw=1, alpha=0.6)
        ejes[i, 1].annotate(f"Nyquist {fs / 2:g} Hz", xy=(fs / 2, 1.0),
                            xycoords=("data", "axes fraction"), xytext=(-4, -10),
                            textcoords="offset points", fontsize=8, ha="right", color="black")
        ejes[i, 1].set(title=f'{s["nombre"]} — espectro', ylabel="|X(f)|", xlim=(0, fs / 2))
        ejes[i, 1].grid(True, alpha=0.3)

    ejes[3, 0].set_xlabel("Tiempo t (s)")
    ejes[3, 1].set_xlabel("Frecuencia f (Hz)")
    fig.tight_layout()

    def consultar(evento):
        datos = lecturas.get(evento.inaxes)
        if datos is None:
            return
        dominio, s, xs, ys, observadas = datos
        j = int(np.abs(xs - evento.xdata).argmin())   # el punto más cercano al cursor
        if evento.inaxes not in marcas:
            marcas[evento.inaxes] = (
                evento.inaxes.plot([], [], "o", ms=7, mec="black", mew=1.2)[0],
                evento.inaxes.annotate("", xy=(0, 0), xytext=(14, 12), textcoords="offset points",
                                       fontsize=8, family="monospace", color="white",
                                       bbox=dict(boxstyle="round,pad=0.35", fc="#111111",
                                                 ec="white", alpha=0.9, lw=0.5)),
            )
        punto, etiqueta = marcas[evento.inaxes]
        punto.set_data([xs[j]], [ys[j]])
        punto.set_color(s["color"])
        # la caja se abre hacia el lado del panel con más espacio libre
        a_la_izquierda = xs[j] < sum(evento.inaxes.get_xlim()) / 2
        etiqueta.set_position((xs[j], ys[j]))
        etiqueta.xytext = (14 if a_la_izquierda else -14, 12)
        etiqueta.set_ha("left" if a_la_izquierda else "right")
        etiqueta.set_text(
            f'{s["nombre"]}\nt = {xs[j]:.3f} s\nx(t) = {ys[j]:.3f}\n'
            f"A (señal) = {np.abs(s['x']).max():.2f}\nF = {observadas}"
            if dominio == "tiempo" else
            f'{s["nombre"]}\nf = {xs[j]:.2f} Hz\n|X(f)| = {ys[j]:.3f}\nF = {observadas}'
        )
        evento.inaxes.figure.canvas.draw_idle()

    def salir(evento):
        if evento.inaxes in marcas:
            marcas[evento.inaxes][1].set_visible(False)
            evento.inaxes.figure.canvas.draw_idle()

    fig.canvas.mpl_connect("motion_notify_event", consultar)
    fig.canvas.mpl_connect("axes_leave_event", salir)
    return fig


def main():
    resumen()
    dibujar(FS_NYQUIST, f"Muestreo que cumple Nyquist-Shannon (fs = {FS_NYQUIST:g} Hz)", "black")
    dibujar(FS_ALIASING, f"Aliasing: fs = {FS_ALIASING:g} Hz < 2 · f_max = {2 * F_MAX:g} Hz", "red")
    plt.show()


if __name__ == "__main__":
    main()
