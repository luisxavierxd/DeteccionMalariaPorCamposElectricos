import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

# Ruta de la carpeta donde está el script
carpeta = os.path.dirname(__file__)

# Parámetros físicos
k = 9e9
q = 1
d = 0.5

# Posiciones de las cargas
q1_pos = np.array([-d, 0])
q2_pos = np.array([ d, 0])

# Malla de puntos
x, y = np.meshgrid(np.arange(-3, 3.1, 0.1), np.arange(-3, 3.1, 0.1))

# Campo debido a carga positiva
rx1 = x - q1_pos[0]; ry1 = y - q1_pos[1]
r1  = np.sqrt(rx1**2 + ry1**2)
Ex1 = k * q * rx1 / r1**3
Ey1 = k * q * ry1 / r1**3

# Campo debido a carga negativa
rx2 = x - q2_pos[0]; ry2 = y - q2_pos[1]
r2  = np.sqrt(rx2**2 + ry2**2)
Ex2 = -k * q * rx2 / r2**3
Ey2 = -k * q * ry2 / r2**3

# Campo total
Ex = Ex1 + Ex2
Ey = Ey1 + Ey2

# Normalizar para quiver (solo dirección, como MATLAB auto-escala)
E_mag = np.sqrt(Ex**2 + Ey**2)
Ex_n = Ex / E_mag
Ey_n = Ey / E_mag

# ── Figura ────────────────────────────────────────────────────────────────
fig, ax = plt.subplots(figsize=(6, 6))

ax.quiver(x, y, Ex_n, Ey_n,
          color='#0072BD', linewidth=0.8, scale=40, width=0.003)

ax.plot(*q1_pos, 'o', markersize=10, color='red',   label=r'$+q$')
ax.plot(*q2_pos, 'o', markersize=10, color='blue',  label=r'$-q$')

ax.set_aspect('equal')
ax.set_xlim(-1, 1)
ax.set_ylim(-1, 1)
ax.grid(True, linestyle='--', alpha=0.3)
ax.tick_params(labelsize=12)
ax.xaxis.set_major_formatter(ticker.FormatStrFormatter('%.1f'))
ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%.1f'))

ax.set_title(r'Campo eléctrico: dipolo de cargas puntuales',
             fontsize=15)
ax.set_xlabel(r'$x$ (m)', fontsize=14)
ax.set_ylabel(r'$y$ (m)', fontsize=14)
ax.legend(fontsize=11)

plt.tight_layout()
plt.savefig(os.path.join(carpeta, 'campo_electrico_dipolo.png'), dpi=150, bbox_inches='tight')
plt.show()