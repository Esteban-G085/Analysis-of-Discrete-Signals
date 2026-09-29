"""Autocomprobación del muestreo, el aliasing y la FFT. Ejecutar: python test_senales_discretas.py"""

import numpy as np

from señales_discretas import (
    F_MAX,
    FS_ALIASING,
    FS_NYQUIST,
    SENALES,
    f_alias,
    frecuencias_observadas,
    muestrear,
)


def pico(f, mag, objetivo):
    return mag[int(np.abs(f - objetivo).argmin())]


def test_nyquist_declarado():
    assert F_MAX == 30.0
    assert FS_NYQUIST >= 2 * F_MAX
    assert FS_ALIASING < 2 * F_MAX


def test_tamano_de_las_muestras():
    for fs in (FS_NYQUIST, FS_ALIASING):
        N = int(fs * 2.0)
        t_m, x_m, f, mag = muestrear(SENALES[0]["x"], fs)
        assert len(t_m) == len(x_m) == N
        assert np.isclose(t_m[1] - t_m[0], 1 / fs)
        assert len(f) == len(mag) == N // 2


def test_f_alias():
    assert f_alias(30.0, 25.0) == 5.0    # 30 Hz se pliega a 5 Hz
    assert f_alias(20.0, 25.0) == 5.0    # 20 Hz se pliega justo encima de 5 Hz
    assert f_alias(17.0, 25.0) == 8.0
    assert f_alias(10.0, 25.0) == 10.0   # bajo el límite de 12.5 Hz: no hay aliasing
    assert f_alias(30.0, FS_NYQUIST) == 30.0


def test_escala_de_la_fft():
    # una senoidal pura de amplitud 1.5 da exactamente 1.5: la escala 2/N es correcta
    _, _, f, mag = muestrear(SENALES[0]["x"], FS_NYQUIST)
    assert np.isclose(pico(f, mag, 10.0), 1.5)


def test_picos_sin_aliasing():
    for s in SENALES:
        _, _, f, mag = muestrear(s["x"], FS_NYQUIST)
        for f_esperada, A in s["componentes"]:
            # 2/N reproduce A en señales periódicas; S4 decae y reparte su energía en
            # bins vecinos, de modo que su pico queda en 0.76 de A (sección 2.4)
            assert pico(f, mag, f_esperada) >= 0.3 * A, f'{s["nombre"]} @ {f_esperada} Hz'


def test_picos_con_aliasing():
    for s in SENALES:
        _, _, f, mag = muestrear(s["x"], FS_ALIASING)
        # el pico del espectro cae en la frecuencia plegada, no en la original
        dominante = max(s["componentes"], key=lambda c: c[1])[0]
        assert f[mag.argmax()] == f_alias(dominante, FS_ALIASING), s["nombre"]


def test_s1_no_sufren_aliasing():
    # 10 Hz está bajo el límite de 12.5 Hz: a 25 Hz se sigue leyendo como 10 Hz
    t_m, x_m, _, _ = muestrear(SENALES[0]["x"], FS_ALIASING)
    assert np.allclose(x_m, 1.5 * np.sin(2 * np.pi * 10.0 * t_m))
    assert "10 Hz" in frecuencias_observadas(SENALES[0], FS_ALIASING)
    assert "30 Hz -> 5 Hz" in frecuencias_observadas(SENALES[3], FS_ALIASING)


if __name__ == "__main__":
    for nombre, prueba in sorted(globals().items()):
        if nombre.startswith("test_"):
            prueba()
            print("ok", nombre)
    print("todo correcto")
