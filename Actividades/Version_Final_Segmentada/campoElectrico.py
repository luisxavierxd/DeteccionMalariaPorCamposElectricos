import numpy as np

# ── Geometría de las placas ───────────────────────────────────────────────────
DX_PLACA   = -1.0
LARGO_ROJA =  4
LARGO_AZUL =  1
ANCHO      = -0.2

X_ROJA = DX_PLACA
X_AZUL = DX_PLACA + 2

NQ_VIS = 500  # cargas por placa para el campo de fondo (visualización)


def calcular_campo(resolucion=1000, nq=NQ_VIS):
    """Calcula Ex, Ey y V en una malla cuadrada [-3,3]×[-3,3].

    Usa k=1, q=1 (normalizados) para que la visualización de potencial
    sea adimensional y comparable entre ejecuciones.

    Retorna: xx, yy, X, Y, Ex, Ey, V
    """
    yp = np.linspace(-LARGO_ROJA / 2, LARGO_ROJA / 2, nq)
    xp = np.full(nq, X_ROJA)
    yn = np.linspace(-LARGO_AZUL / 2, LARGO_AZUL / 2, nq)
    xn = np.full(nq, X_AZUL)

    xx = np.linspace(-3, 3, resolucion)
    yy = np.linspace(-3, 3, resolucion)
    X, Y = np.meshgrid(xx, yy)

    Ex = np.zeros_like(X)
    Ey = np.zeros_like(Y)
    V  = np.zeros_like(X)

    for i in range(nq):
        # Placa positiva — campo
        rx = X - xp[i];  ry = Y - yp[i]
        rp = np.sqrt(rx**2 + ry**2);  rp[rp < 0.15] = np.nan
        Ex += rx / rp**3
        Ey += ry / rp**3

        # Placa negativa — campo
        rx = X - xn[i];  ry = Y - yn[i]
        rn = np.sqrt(rx**2 + ry**2);  rn[rn < 0.15] = np.nan
        Ex -= rx / rn**3
        Ey -= ry / rn**3

        # Potencial (clamp en lugar de NaN para evitar artefactos en contourf)
        rp2 = np.sqrt((X - xp[i])**2 + (Y - yp[i])**2);  rp2[rp2 < 0.15] = 0.15
        rn2 = np.sqrt((X - xn[i])**2 + (Y - yn[i])**2);  rn2[rn2 < 0.15] = 0.15
        V += 1.0 / rp2 - 1.0 / rn2

    return xx, yy, X, Y, Ex, Ey, V
