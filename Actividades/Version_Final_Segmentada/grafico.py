import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

from campoElectrico import X_ROJA, X_AZUL, LARGO_ROJA, LARGO_AZUL, ANCHO


def graficar(xx, yy, X, Y, Ex, Ey, V, trayectorias, ruta_salida):
    """Genera la figura con campo de fondo y trayectorias superpuestas.

    Campo de fondo: potencial (contourf 'bone') + equipotenciales (contour) +
    líneas de campo (streamplot verde).
    Trayectorias: verde para células sanas, rojo para infectadas.
    """
    fig, ax = plt.subplots(figsize=(9, 7))

    V_plot = np.clip(V, -500, 500)
    cf = ax.contourf(X, Y, V_plot, levels=60, cmap='bone')
    plt.colorbar(cf, ax=ax, label='Potencial eléctrico V')
    ax.contour(X, Y, V_plot, levels=10, colors='white', linewidths=0.5)

    with np.errstate(invalid='ignore'):
        Ex_s = np.where(np.isfinite(Ex), Ex, 0.0)
        Ey_s = np.where(np.isfinite(Ey), Ey, 0.0)

    # Enmascarar el interior de las placas para el streamplot
    ancho_abs = abs(ANCHO)
    mask_roja = ((X >= X_ROJA - ancho_abs) & (X <= X_ROJA + ancho_abs) &
                 (Y >= -LARGO_ROJA / 2)     & (Y <= LARGO_ROJA / 2))
    mask_azul = ((X >= X_AZUL - ancho_abs) & (X <= X_AZUL + ancho_abs) &
                 (Y >= -LARGO_AZUL / 2)     & (Y <= LARGO_AZUL / 2))
    Ex_s[mask_roja | mask_azul] = np.nan
    Ey_s[mask_roja | mask_azul] = np.nan

    ax.streamplot(xx, yy, Ex_s, Ey_s,
                  color='#00cc44', linewidth=1.0, density=1.5, arrowsize=1.2)

    ax.add_patch(Rectangle((X_ROJA - ancho_abs, -LARGO_ROJA / 2),
                            2 * ancho_abs, LARGO_ROJA,
                            color='red', alpha=0.7, label='Placa + (roja)'))
    ax.add_patch(Rectangle((X_AZUL - ancho_abs, -LARGO_AZUL / 2),
                            2 * ancho_abs, LARGO_AZUL,
                            color='blue', alpha=0.7, label='Placa − (azul)'))

    vistos = set()
    for tray in trayectorias:
        label = tray['tipo'] if tray['tipo'] not in vistos else '_nolegend_'
        vistos.add(tray['tipo'])
        ax.plot(tray['xs'], tray['ys'],
                color=tray['color'], linewidth=1.2, alpha=0.9,
                label=label, zorder=5)

    ax.set_xlim(-3, 3);  ax.set_ylim(-3, 3)
    ax.set_aspect('equal')
    ax.set_xlabel(r'$x$ (m)');  ax.set_ylabel(r'$y$ (m)')
    ax.set_title('Trayectorias de células (sanas vs infectadas)', fontsize=13)
    ax.legend(loc='upper right')
    ax.grid(True, linestyle='--', alpha=0.2)

    plt.tight_layout()
    plt.savefig(ruta_salida, dpi=150, bbox_inches='tight')
    print(f"Imagen guardada: {ruta_salida}")
