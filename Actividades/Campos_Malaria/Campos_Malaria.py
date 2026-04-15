import os
import numpy as np
import matplotlib                      # Debe importarse antes de pyplot.
matplotlib.use("Agg")                  # Backend sin GUI (escritura a archivo).
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from matplotlib.patches import Rectangle

carpeta = os.path.dirname(os.path.abspath(__file__))

# ── Parámetros físicos ────────────────────────────────────────────────────
k  = 1
q  = 1
Nq = 500

# ── Geometría ─────────────────────────────────────────────────────────────
v           = 0.5
dx_placa    = -1.0
largo_roja  = 4
largo_azul  = 1
ancho       = -0.2

x_roja = dx_placa
x_azul = dx_placa + 2

# ── Posiciones de las cargas ──────────────────────────────────────────────
yp = np.linspace(-largo_roja/2, largo_roja/2, Nq)
xp = np.full(Nq, x_roja);  zp = np.zeros(Nq)

yn = np.linspace(-largo_azul/2, largo_azul/2, Nq)
xn = np.full(Nq, x_azul);  zn = np.zeros(Nq)

# ── Malla 2D ─────────────────────────────────────────────────────────────
xx = np.linspace(-3, 3, 1000)
yy = np.linspace(-3, 3, 1000)
X, Y = np.meshgrid(xx, yy)
Z    = np.zeros_like(X)

Ex = np.zeros_like(X)
Ey = np.zeros_like(Y)
Ez = np.zeros_like(Z)
V  = np.zeros_like(X)

for k_i in range(Nq):
    # ── Carga positiva (+q) ──
    rx = X - xp[k_i];  ry = Y - yp[k_i];  rz = Z - zp[k_i]
    rp = np.sqrt(rx**2 + ry**2 + rz**2)
    rp[rp < 0.15] = np.nan
    Ex += k * q * rx / rp**3
    Ey += k * q * ry / rp**3
    Ez += k * q * rz / rp**3

    # ── Carga negativa (-q) ──
    rx = X - xn[k_i];  ry = Y - yn[k_i];  rz = Z - zn[k_i]
    rn = np.sqrt(rx**2 + ry**2 + rz**2)
    rn[rn < 0.15] = np.nan
    Ex += k * (-q) * rx / rn**3
    Ey += k * (-q) * ry / rn**3
    Ez += k * (-q) * rz / rn**3

    # ── Potencial eléctrico V = k*q/rp - k*q/rn ──────────────────────────
    rp2 = np.sqrt((X - xp[k_i])**2 + (Y - yp[k_i])**2)
    rp2[rp2 < 0.15] = 0.15  # evitar singularidades/arreglo cuadrado blanco
    rn2 = np.sqrt((X - xn[k_i])**2 + (Y - yn[k_i])**2)
    rn2[rn2 < 0.15] = 0.15  # evitar singularidades
    V += k * q / rp2 - k * q / rn2

# ── Figura 2D — Potencial + equipotenciales + streamlines ─────────────
fig2, ax2 = plt.subplots(figsize=(9, 7))

V_plot = np.clip(V, -500, 500)

# pcolor equivalente → contourf con colormap 'bone'
cf = ax2.contourf(X, Y, V_plot, levels=60, cmap='bone')
plt.colorbar(cf, ax=ax2, label='Potencial eléctrico V')

# Líneas equipotenciales blancas (contour(...,10,'-w'))
ax2.contour(X, Y, V_plot, levels=10, colors='white', linewidths=0.5)

# streamplot — máscara para que las líneas no atraviesen las placas
with np.errstate(invalid='ignore'):
    Ex_s = np.where(np.isfinite(Ex), Ex, 0.0)
    Ey_s = np.where(np.isfinite(Ey), Ey, 0.0)

mask_roja = (
    (X >= x_roja - abs(ancho)) & (X <= x_roja + abs(ancho)) &
    (Y >= -largo_roja/2)       & (Y <= largo_roja/2)
)
mask_azul = (
    (X >= x_azul - abs(ancho)) & (X <= x_azul + abs(ancho)) &
    (Y >= -largo_azul/2)       & (Y <= largo_azul/2)
)
Ex_s[mask_roja | mask_azul] = np.nan
Ey_s[mask_roja | mask_azul] = np.nan

ax2.streamplot(xx, yy, Ex_s, Ey_s,
               color='#00cc44', linewidth=1.0, density=1.5, arrowsize=1.2)

# Placas
ax2.add_patch(Rectangle((x_roja-abs(ancho), -largo_roja/2), 2*abs(ancho), largo_roja,
                         color='red', alpha=0.7, label='Placa + (roja)'))
ax2.add_patch(Rectangle((x_azul-abs(ancho), -largo_azul/2), 2*abs(ancho), largo_azul,
                         color='blue', alpha=0.7, label='Placa − (azul)'))

ax2.set_xlim(-3,3); ax2.set_ylim(-3,3)
ax2.set_xlabel(r'$x$ (m)'); ax2.set_ylabel(r'$y$ (m)')
ax2.set_title('Potencial eléctrico y líneas equipotenciales', fontsize=13)
ax2.legend(loc='upper right')
ax2.set_aspect('equal')
ax2.grid(True, linestyle='--', alpha=0.2)
plt.tight_layout()
plt.savefig(os.path.join(carpeta, 'campo_placas_malaria.png'), dpi=150, bbox_inches='tight')

plt.show()
print("Listo: campo_placas_malaria.png")