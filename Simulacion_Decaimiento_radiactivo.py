# This is a sample Python script.

# ------------------------------------------------------------------------
# LIBRERÍAS
# ------------------------------------------------------------------------
import numpy as np  # cálculo numérico vectorizado (arrays, exp, log, etc.)
import matplotlib.pyplot as plt  # generación de gráficas científicas
from scipy.integrate import cumulative_trapezoid  # integración numérica (regla del trapecio) para acumular energía

# ------------------------------------------------------------------------
# 1. CONSTANTES FÍSICAS Y FACTORES DE CONVERSIÓN
# ------------------------------------------------------------------------
AVOGADRO = 6.02214076e23  # [núcleos/mol] número de Avogadro: convierte masa -> número de átomos
MEV_A_JOULE = 1.602176634e-13  # [J/MeV] factor de conversión de MeV (unidad nuclear típica) a Joules (unidad SI de energía)
SEGUNDOS_POR_ANIO = 360 * 24 * 3600  # [s] año juliano, usado para convertir vidas medias dadas en años a segundos
SEGUNDOS_POR_DIA = 86400  # [s] un día en segundos, usado para el I-131 (vida media dada en días)

# ------------------------------------------------------------------------
# 2. DATOS NUCLEARES DE CADA ISÓTOPO
# ------------------------------------------------------------------------
# Vida media (T½): tiempo físico medido experimentalmente en el que decae la
# mitad de una muestra. Masa atómica: para convertir gramos -> núcleos.
# Energía por decaimiento (valor Q aproximado, en MeV): energía térmica que
# se considera depositada en cada desintegración (simplificación educativa;
# en la realidad una fracción se pierde en neutrinos no detectables).
isotopos = {
    "Uranio-235": {
        "vida_media_s": 7.038e8 * SEGUNDOS_POR_ANIO,  # T½ ≈ 703.8 millones de años (decaimiento alfa)
        "masa_atomica": 235.0439,  # [g/mol]
        "energia_por_decaimiento_MeV": 4.678,  # [MeV] energía total de la desintegración alfa (Q)
        "color": "#1f77b4",  # color para graficar (azul)
    },
    "Cesio-137": {
        "vida_media_s": 30.17 * SEGUNDOS_POR_ANIO,  # T½ ≈ 30.17 años (decaimiento beta menos)
        "masa_atomica": 136.907,  # [g/mol]
        "energia_por_decaimiento_MeV": 1.176,  # [MeV] Q total (beta + desexcitación gamma del Ba-137m)
        "color": "#ff7f0e",  # color para graficar (naranja)
    },
    "Yodo-131": {
        "vida_media_s": 8.02 * SEGUNDOS_POR_DIA,  # T½ ≈ 8.02 días (decaimiento beta menos, uso médico/reactores)
        "masa_atomica": 130.906,  # [g/mol]
        "energia_por_decaimiento_MeV": 0.971,  # [MeV] Q total (beta + gamma característico de 364 keV)
        "color": "#2ca02c",  # color para graficar (verde)
    },
}

masa_inicial_g = 10.0  # [g] masa inicial arbitraria e idéntica para los 3 isótopos, para poder compararlos en igualdad de condiciones

# ------------------------------------------------------------------------
# 3. MALLA TEMPORAL (ESCALA LOGARÍTMICA)
# ------------------------------------------------------------------------
# Como las vidas medias difieren en varios órdenes de magnitud (~8 días vs
# ~700 millones de años), se usa una malla de tiempo logarítmica: esto
# permite observar en una sola gráfica tanto el decaimiento casi instantáneo
# del I-131 como el decaimiento imperceptible del U-235 durante el mismo
# intervalo de observación humano (200 años).
t_min_anios = 1e-3  # [años] ≈ 8.77 horas, tiempo inicial (no puede ser 0 en escala log)
t_max_anios = 100  # [años] horizonte de observación (una vida humana varias veces)
n_puntos = 4000  # número de puntos de muestreo temporal (resolución de la simulación)

# np.logspace genera puntos equiespaciados en escala log10 entre dos exponentes
t_anios = np.logspace(np.log10(t_min_anios), np.log10(t_max_anios), n_puntos)  # [años]
t_segundos = t_anios * SEGUNDOS_POR_ANIO  # conversión a segundos (unidad base de las constantes de decaimiento)

# ------------------------------------------------------------------------
# 4. SIMULACIÓN FÍSICA: DECAIMIENTO Y ENERGÍA, ISÓTOPO POR ISÓTOPO
# ------------------------------------------------------------------------
resultados = {}  # diccionario donde guardamos las curvas calculadas de cada isótopo

for nombre, datos in isotopos.items():

    # --- Paso 4.1: constante de decaimiento λ a partir de la vida media ---
    # λ representa la probabilidad de decaimiento por núcleo y por segundo.
    vida_media_s = datos["vida_media_s"]
    lambda_decaimiento = np.log(2) / vida_media_s  # [1/s], derivado de N(T½) = N0/2

    # --- Paso 4.2: número inicial de núcleos N0 a partir de la masa ---
    # N0 = (masa / masa_atomica) da los moles; multiplicado por Avogadro da átomos.
    N0 = (masa_inicial_g / datos["masa_atomica"]) * AVOGADRO  # [núcleos]

    # --- Paso 4.3: energía liberada por cada decaimiento, en Joules ---
    E_por_decaimiento_J = datos["energia_por_decaimiento_MeV"] * MEV_A_JOULE  # [J]

    # --- Paso 4.4: número de núcleos restantes N(t), ley exponencial ---
    # Esta es la solución analítica exacta de dN/dt = -λN.
    N_t = N0 * np.exp(-lambda_decaimiento * t_segundos)  # [núcleos] en cada instante t

    # --- Paso 4.5: actividad A(t), es decir, desintegraciones por segundo ---
    # A(t) = λ·N(t); físicamente es lo que mediría un contador Geiger (en Bq).
    A_t = lambda_decaimiento * N_t  # [Bq] = [desintegraciones/s]

    # --- Paso 4.6: potencia térmica instantánea liberada P(t) ---
    # Cada desintegración deposita E_por_decaimiento_J julios de energía.
    P_t = A_t * E_por_decaimiento_J  # [W] = [J/s]

    # --- Paso 4.7: ENERGÍA ACUMULADA mediante INTEGRACIÓN NUMÉRICA (SciPy) ---
    # cumulative_trapezoid aproxima el área bajo la curva P(t) usando la regla
    # del trapecio en cada subintervalo de la malla temporal, acumulando el
    # resultado. initial=0.0 fija E(t_min) = 0 como condición de referencia.
    E_acumulada_J = cumulative_trapezoid(P_t, t_segundos, initial=0.0)  # [J]

    # --- Paso 4.8: validación contra la solución analítica exacta ---
    # Integrando la ley exponencial entre t_min (primer punto de la malla) y
    # un tiempo t se obtiene una fórmula cerrada:
    #   E(t) = E_decaimiento * (N(t_min) - N(t))
    # OJO: usamos N(t_min) y no N0 como referencia, porque la malla temporal
    # (logarítmica, por necesidad) arranca en t_min > 0, no en t = 0. Para
    # isótopos de vida media muy corta (p. ej. I-131, T½ ≈ 8 días) ya ha
    # decaído una fracción no despreciable entre t=0 y t_min ≈ 8.8 horas, así
    # que comparar contra N0 introduciría un desfase artificial y no un error
    # real de la integración numérica.
    E_analitica_J = E_por_decaimiento_J * (N_t[0] - N_t)  # [J], referenciada a t_min

    # Evitamos división por cero / cancelación catastrófica cuando la energía
    # liberada todavía es despreciable (típico en U-235 al inicio de la malla).
    energia_maxima_posible = E_por_decaimiento_J * N0  # [J] si decayera el 100% de la muestra (desde t=0)
    mascara_valida = E_analitica_J > (energia_maxima_posible * 1e-9)
    if np.any(mascara_valida):
        error_relativo = np.max(
            np.abs(E_acumulada_J[mascara_valida] - E_analitica_J[mascara_valida])
            / E_analitica_J[mascara_valida]
        )
    else:
        error_relativo = np.nan  # decaimiento tan lento que no hay energía significativa que comparar

    # Guardamos todo lo necesario para graficar y reportar resultados
    resultados[nombre] = {
        "lambda": lambda_decaimiento,
        "vida_media_s": vida_media_s,
        "N0": N0,
        "N_t": N_t,
        "E_acumulada_J": E_acumulada_J,
        "energia_maxima_posible_J": energia_maxima_posible,
        "error_relativo_integracion": error_relativo,
        "color": datos["color"],
    }

# ------------------------------------------------------------------------
# 5. REPORTE NUMÉRICO EN CONSOLA
# ------------------------------------------------------------------------
print("=" * 78)
print(f"RESUMEN DE SIMULACIÓN (masa inicial de cada isótopo: {masa_inicial_g} g)")
print("=" * 78)
for nombre, r in resultados.items():
    porcentaje_decaido = (1 - r["N_t"][-1] / r["N0"]) * 100  # % de núcleos ya desintegrados al final de la simulación
    energia_liberada_final = r["E_acumulada_J"][-1]
    print(f"\n{nombre}:")
    print(f"  Vida media               : {r['vida_media_s'] / SEGUNDOS_POR_ANIO:.4e} años")
    print(f"  Núcleos iniciales (N0)   : {r['N0']:.4e}")
    print(f"  % decaído en {t_max_anios} años   : {porcentaje_decaido:.6f} %")
    print(f"  Energía acumulada final  : {energia_liberada_final:.4e} J")
    print(f"  Energía máxima posible   : {r['energia_maxima_posible_J']:.4e} J (si decayera el 100%)")
    if not np.isnan(r["error_relativo_integracion"]):
        print(f"  Error rel. integ. numérica vs analítica: {r['error_relativo_integracion']:.2e}")
    else:
        print("  Error rel. integ. numérica vs analítica: N/A (decaimiento insignificante en este rango)")

# ------------------------------------------------------------------------
# 6. VISUALIZACIÓN: DOS GRÁFICAS PROFESIONALES EN UNA SOLA IMAGEN
# ------------------------------------------------------------------------
plt.rcParams.update({
    "font.size": 11,
    "axes.titlesize": 13,
    "axes.titleweight": "bold",
    "axes.labelsize": 11,
    "figure.facecolor": "white",
    "axes.facecolor": "white",
})

fig, (ax_decay, ax_energy) = plt.subplots(1, 2, figsize=(15, 6.5))

# ---- Subgráfica izquierda: curvas de decaimiento N(t)/N0 en escala log-log ----
# Escala logarítmica en ambos ejes: permite comparar isótopos cuyas vidas
# medias difieren en más de 8 órdenes de magnitud en una sola gráfica legible.
for nombre, r in resultados.items():
    fraccion_restante = r["N_t"] / r["N0"]
    # np.clip evita valores exactamente 0.0 (underflow numérico de e^-x para x muy grande)
    # que no pueden representarse en escala logarítmica (log(0) = -inf).
    fraccion_restante_clip = np.clip(fraccion_restante, 1e-300, None)
    vida_media_anios = r["vida_media_s"] / SEGUNDOS_POR_ANIO
    ax_decay.loglog(
        t_anios, fraccion_restante_clip,
        color=r["color"], linewidth=2.2,
        label=f"{nombre}  (T½ = {vida_media_anios:.3g} años)",
    )
    # línea vertical punteada marcando la vida media de cada isótopo, como referencia visual
    ax_decay.axvline(vida_media_anios, color=r["color"], linestyle=":", alpha=0.5, linewidth=1.3)

ax_decay.set_xlabel("Tiempo transcurrido (años, escala logarítmica)")
ax_decay.set_ylabel("Fracción de núcleos restantes  N(t)/N₀  (escala log)")
ax_decay.set_title("Decaimiento radiactivo: N(t) = N₀·e^(-λt)")
ax_decay.grid(True, which="both", linestyle=":", alpha=0.5)
ax_decay.legend(loc="lower left", fontsize=9, framealpha=0.9)
ax_decay.set_ylim(1e-8, 2)  # límite inferior razonable para visualizar sin ruido numérico extremo

# ---- Subgráfica derecha: energía térmica acumulada liberada ----
for nombre, r in resultados.items():
    energia_clip = np.clip(r["E_acumulada_J"], 1e-6, None)  # evita ceros exactos en escala log
    ax_energy.loglog(
        t_anios, energia_clip,
        color=r["color"], linewidth=2.2, label=nombre,
    )
    # línea horizontal discontinua: energía máxima posible si decayera el 100% de la muestra
    ax_energy.axhline(r["energia_maxima_posible_J"], color=r["color"], linestyle="--", alpha=0.35, linewidth=1.3)

ax_energy.set_xlabel("Tiempo transcurrido (años, escala logarítmica)")
ax_energy.set_ylabel("Energía térmica acumulada liberada, E(t)  [J]  (escala log)")
ax_energy.set_title("Energía acumulada (integración numérica: SciPy trapecio)")
ax_energy.grid(True, which="both", linestyle=":", alpha=0.5)
ax_energy.legend(loc="lower right", fontsize=9, framealpha=0.9)

fig.suptitle(
    "Simulación de Decaimiento Radiactivo: U-235, Cs-137 y I-131 (muestra inicial de 10 g cada uno)",
    fontsize=14, fontweight="bold", y=1.02,
)
plt.tight_layout()

# ------------------------------------------------------------------------
# 7. GUARDADO DE LA FIGURA
# ------------------------------------------------------------------------
ruta_salida = "decaimiento_radiactivo.png"
plt.savefig(ruta_salida, dpi=160, bbox_inches="tight")
print(f"\nGráfica guardada en: {ruta_salida}")