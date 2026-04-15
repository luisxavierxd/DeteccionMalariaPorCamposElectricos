import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

carpeta = os.path.dirname(os.path.abspath(__file__))

# ── Geometría (compartida) ────────────────────────────────────────────────
dx_placa   = -1.0
largo_roja = 4
largo_azul = 1
ancho      = -0.2

x_roja = dx_placa
x_azul = dx_placa + 2

Nq = 500
yp = np.linspace(-largo_roja/2, largo_roja/2, Nq)
xp = np.full(Nq, x_roja)
yn = np.linspace(-largo_azul/2, largo_azul/2, Nq)
xn = np.full(Nq, x_azul)

# ── Campo normalizado (k=1, q=1) para visualización ───────────────────────
k_vis, q_vis = 1, 1

xx = np.linspace(-3, 3, 1000)
yy = np.linspace(-3, 3, 1000)
X, Y = np.meshgrid(xx, yy)

Ex = np.zeros_like(X)
Ey = np.zeros_like(Y)
V  = np.zeros_like(X)

for i in range(Nq):
    rx = X - xp[i]; ry = Y - yp[i]
    rp = np.sqrt(rx**2 + ry**2); rp[rp < 0.15] = np.nan
    Ex += k_vis * q_vis * rx / rp**3
    Ey += k_vis * q_vis * ry / rp**3

    rx = X - xn[i]; ry = Y - yn[i]
    rn = np.sqrt(rx**2 + ry**2); rn[rn < 0.15] = np.nan
    Ex += k_vis * (-q_vis) * rx / rn**3
    Ey += k_vis * (-q_vis) * ry / rn**3

    rp2 = np.sqrt((X-xp[i])**2 + (Y-yp[i])**2); rp2[rp2 < 0.15] = 0.15
    rn2 = np.sqrt((X-xn[i])**2 + (Y-yn[i])**2); rn2[rn2 < 0.15] = 0.15
    V += k_vis*q_vis/rp2 - k_vis*q_vis/rn2

# ── Simulación de partículas (k=9e9, q=1e-5, unidades físicas) ───────────
k_fis, q_fis = 9e9, 1e-5
Nq_sim = 150   # menos puntos para la dinámica (más rápido, suficiente precisión)

yp_s = np.linspace(-largo_roja/2, largo_roja/2, Nq_sim)
xp_s = np.full(Nq_sim, x_roja)
yn_s = np.linspace(-largo_azul/2, largo_azul/2, Nq_sim)
xn_s = np.full(Nq_sim, x_azul)

N_particulas = 30
dt = 0.02
m  = 1
max_steps = 2000
v_ini = 10

q_sana      = 0.8e-6
q_infectada = 1.4e-6
ruido       = 0.3e-6

np.random.seed(42)
trayectorias = []

for i in range(N_particulas):
    x = np.random.uniform(-0.2, 0.2)
    y = 2.5
    vx, vy = 0.0, -v_ini

    if np.random.rand() < 0.5:
        q_base = q_sana;      color = 'lime';  label = 'Sana'
    else:
        q_base = q_infectada; color = 'red';   label = 'Infectada'

    q_part = max(q_base + np.random.normal(0, ruido), 1e-9)
    xs, ys = [x], [y]

    for _ in range(max_steps):
        rxp = x - xp_s; ryp = y - yp_s
        rp  = np.sqrt(rxp**2 + ryp**2); mask_p = rp > 0.1
        Fx = np.sum(k_fis * q_part * q_fis  * rxp[mask_p] / rp[mask_p]**3)
        Fy = np.sum(k_fis * q_part * q_fis  * ryp[mask_p] / rp[mask_p]**3)

        rxn = x - xn_s; ryn = y - yn_s
        rn  = np.sqrt(rxn**2 + ryn**2); mask_n = rn > 0.1
        Fx += np.sum(k_fis * q_part * (-q_fis) * rxn[mask_n] / rn[mask_n]**3)
        Fy += np.sum(k_fis * q_part * (-q_fis) * ryn[mask_n] / rn[mask_n]**3)

        vx += (Fx/m)*dt; vy += (Fy/m)*dt
        x  += vx*dt;     y  += vy*dt
        xs.append(x); ys.append(y)

        if y < -2.5 or abs(x) > 3:
            break

    trayectorias.append((xs, ys, color, label))

# ── Visualización idéntica a referencia + trayectorias encima ─────────────
fig, ax = plt.subplots(figsize=(9, 7))

V_plot = np.clip(V, -500, 500)
cf = ax.contourf(X, Y, V_plot, levels=60, cmap='bone')
plt.colorbar(cf, ax=ax, label='Potencial eléctrico V')

ax.contour(X, Y, V_plot, levels=10, colors='white', linewidths=0.5)

with np.errstate(invalid='ignore'):
    Ex_s = np.where(np.isfinite(Ex), Ex, 0.0)
    Ey_s = np.where(np.isfinite(Ey), Ey, 0.0)

mask_roja = ((X >= x_roja-abs(ancho)) & (X <= x_roja+abs(ancho)) &
             (Y >= -largo_roja/2)      & (Y <= largo_roja/2))
mask_azul = ((X >= x_azul-abs(ancho)) & (X <= x_azul+abs(ancho)) &
             (Y >= -largo_azul/2)      & (Y <= largo_azul/2))
Ex_s[mask_roja | mask_azul] = np.nan
Ey_s[mask_roja | mask_azul] = np.nan

ax.streamplot(xx, yy, Ex_s, Ey_s,
              color='#00cc44', linewidth=1.0, density=1.5, arrowsize=1.2)

ax.add_patch(Rectangle((x_roja-abs(ancho), -largo_roja/2),
                        2*abs(ancho), largo_roja,
                        color='red', alpha=0.7, label='Placa + (roja)'))
ax.add_patch(Rectangle((x_azul-abs(ancho), -largo_azul/2),
                        2*abs(ancho), largo_azul,
                        color='blue', alpha=0.7, label='Placa − (azul)'))

_vistos = set()
for xs, ys, color, label in trayectorias:
    lbl = label if label not in _vistos else '_nolegend_'
    _vistos.add(label)
    ax.plot(xs, ys, color=color, linewidth=1.2, alpha=0.9, label=lbl, zorder=5)

ax.set_xlim(-3, 3); ax.set_ylim(-3, 3)
ax.set_aspect('equal')
ax.set_xlabel(r'$x$ (m)'); ax.set_ylabel(r'$y$ (m)')
ax.set_title('Trayectorias de células (sanas vs infectadas)', fontsize=13)
ax.legend(loc='upper right')
ax.grid(True, linestyle='--', alpha=0.2)

plt.tight_layout()
plt.savefig(os.path.join(carpeta, 'trayectorias_malaria.png'), dpi=150, bbox_inches='tight')
print("Listo: trayectorias_malaria.png")