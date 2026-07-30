# Simulación de Decaimiento Radiactivo

Este proyecto simula el decaimiento radiactivo de tres isótopos (Uranio-235, Cesio-137 y Yodo-131) comparando métodos de integración numérica contra soluciones analíticas.

## Resultados
- **Gráfica:** `decaimiento_radiactivo.png`
- **Masa inicial:** 10.0 g por isótopo.

## Cómo ejecutar
```bash
python mi_script.py
"""
==============================================================================
 SIMULACIÓN DE DECAIMIENTO RADIACTIVO Y LIBERACIÓN DE ENERGÍA TÉRMICA
 Isótopos: Uranio-235 (U-235), Cesio-137 (Cs-137), Yodo-131 (I-131)
==============================================================================
 
FÍSICA DE FONDO (resumen):
---------------------------
Todo núcleo radiactivo decae siguiendo una ley de probabilidad constante por
unidad de tiempo, llamada "constante de decaimiento" (lambda, λ). Esto da
lugar a la célebre ecuación diferencial:
 
        dN/dt = -λ · N(t)
 
cuya solución analítica es la conocida "ley exponencial de decaimiento":
 
        N(t) = N0 · e^(-λt)
 
donde N0 es el número inicial de núcleos. La constante λ se relaciona con la
vida media (T½, el tiempo en que la mitad de los núcleos decae) mediante:
 
        λ = ln(2) / T½
 
Cada vez que un núcleo decae, libera una energía característica (el "valor Q"
de la desintegración), que en la práctica se reparte entre partículas alfa,
beta, gamma y (en el caso beta) neutrinos que escapan sin depositar energía
térmica localmente. Para esta simulación usamos la energía total aproximada
depositada por decaimiento (simplificación pedagógica estándar).
 
La "actividad" A(t) [Bq, desintegraciones/segundo] es:
 
        A(t) = λ · N(t)
 
y la "potencia térmica instantánea" liberada es:
 
        P(t) = A(t) · E_decaimiento   [Watts = Joules/segundo]
 
La energía térmica ACUMULADA hasta un tiempo t es el área bajo la curva de
potencia, es decir, la integral:
 
        E(t) = ∫₀ᵗ P(t') dt'
 
En este script, dicha integral se calcula NUMÉRICAMENTE mediante la regla del
trapecio (scipy.integrate.cumulative_trapezoid), y se valida comparándola
contra la solución analítica exacta E(t) = E_decaimiento · (N0 - N(t)).
==============================================================================
"""
