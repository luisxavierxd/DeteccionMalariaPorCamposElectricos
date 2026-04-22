import csv
import numpy as np

from campoElectrico import X_ROJA, X_AZUL, LARGO_ROJA, LARGO_AZUL

# ── Constantes físicas ────────────────────────────────────────────────────────
K_FIS = 9e9
Q_FIS = 1e-5
NQ_SIM = 150  # cargas por placa; menos que NQ_VIS para reducir coste computacional

# ── Propiedades de las células ────────────────────────────────────────────────
Q_SANA      = 0.8e-6   # carga efectiva célula sana (C)
Q_INFECTADA = 1.4e-6   # carga efectiva célula infectada (C)
RUIDO       = 0.3e-6   # dispersión gaussiana entre células de la misma categoría

# ── Parámetros del integrador ─────────────────────────────────────────────────
N_PARTICULAS = 30
DT           = 0.02
MASA         = 1
MAX_STEPS    = 2000
V_INI        = 10


def _posiciones_placa():
    """Devuelve las coordenadas (x, y) de los puntos discretos de cada placa."""
    yp = np.linspace(-LARGO_ROJA / 2, LARGO_ROJA / 2, NQ_SIM)
    xp = np.full(NQ_SIM, X_ROJA)
    yn = np.linspace(-LARGO_AZUL / 2, LARGO_AZUL / 2, NQ_SIM)
    xn = np.full(NQ_SIM, X_AZUL)
    return xp, yp, xn, yn


def simular(n_particulas=N_PARTICULAS, semilla=42):
    """Integra las trayectorias de n_particulas células con Euler explícito.

    Retorna lista de dicts con claves: xs, ys, tipo, carga, color.
    """
    xp_s, yp_s, xn_s, yn_s = _posiciones_placa()
    np.random.seed(semilla)
    trayectorias = []

    for _ in range(n_particulas):
        x = np.random.uniform(-0.2, 0.2)
        y = 2.5
        vx, vy = 0.0, -V_INI

        if np.random.rand() < 0.5:
            q_base, color, tipo = Q_SANA,      'lime', 'Sana'
        else:
            q_base, color, tipo = Q_INFECTADA, 'red',  'Infectada'

        q_part = max(q_base + np.random.normal(0, RUIDO), 1e-9)
        xs, ys = [x], [y]

        for _ in range(MAX_STEPS):
            # Fuerza de la placa positiva
            rxp = x - xp_s;  ryp = y - yp_s
            rp  = np.sqrt(rxp**2 + ryp**2);  mask_p = rp > 0.1
            Fx = np.sum(K_FIS * q_part *  Q_FIS * rxp[mask_p] / rp[mask_p]**3)
            Fy = np.sum(K_FIS * q_part *  Q_FIS * ryp[mask_p] / rp[mask_p]**3)

            # Fuerza de la placa negativa
            rxn = x - xn_s;  ryn = y - yn_s
            rn  = np.sqrt(rxn**2 + ryn**2);  mask_n = rn > 0.1
            Fx += np.sum(K_FIS * q_part * -Q_FIS * rxn[mask_n] / rn[mask_n]**3)
            Fy += np.sum(K_FIS * q_part * -Q_FIS * ryn[mask_n] / rn[mask_n]**3)

            # Euler explícito: v → r
            vx += (Fx / MASA) * DT
            vy += (Fy / MASA) * DT
            x  += vx * DT
            y  += vy * DT
            xs.append(x);  ys.append(y)

            if y < -2.5 or abs(x) > 3:
                break

        trayectorias.append({
            'xs':    xs,
            'ys':    ys,
            'tipo':  tipo,
            'carga': q_part,
            'color': color,
        })

    return trayectorias


def exportar_csv(trayectorias, ruta):
    """Guarda todas las trayectorias en un CSV.

    Columnas: particula_id, paso, x, y, tipo, carga
    """
    with open(ruta, 'w', newline='', encoding='utf-8') as f:
        writer = csv.writer(f)
        writer.writerow(['particula_id', 'paso', 'x', 'y', 'tipo', 'carga'])
        for pid, tray in enumerate(trayectorias):
            for paso, (xi, yi) in enumerate(zip(tray['xs'], tray['ys'])):
                writer.writerow([pid, paso, xi, yi, tray['tipo'], tray['carga']])
    print(f"CSV exportado: {ruta}")
