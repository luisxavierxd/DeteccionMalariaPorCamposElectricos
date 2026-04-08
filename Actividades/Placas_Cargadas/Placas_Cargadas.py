import os
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

# Ruta de la carpeta donde está el script
carpeta = os.path.dirname(__file__)

# ── Parámetros físicos ────────────────────────────────────────────────────
k  = 9e9          # constante de Coulomb (k=1 para simplificar como la prof sugiere)
k  = 1            # simplificación sugerida en las pistas
q  = 1            # magnitud de carga (signo se aplica abajo)
Nq = 3            # número de cargas por placa

# ── Geometría de las placas ───────────────────────────────────────────────
v        = 0.5    # semiancho de las placas en y,z
dx_placa = -1.0   # desplazamiento global en x → placa positiva en x=-1, negativa en x=1
dy_placa =  0.0
dz_placa =  0.0

largo_roja  = 3   # largo placa positiva (roja)
largo_azul  = 3   # largo placa negativa (azul)
ancho       = -0.2
alto        =  1

x_roja = dx_placa          # x=-1
x_azul = dx_placa + 2      # x=+1

# ── Posiciones de las cargas (linspace a lo largo del largo de la placa) ──
yp = np.linspace(-(1 - 0)*largo_roja/2, (1 - 0)*largo_roja/2, Nq)  # cargas +q
xp = np.full(Nq, x_roja)
zp = np.zeros(Nq)

yn = np.linspace(-(1 - 0)*largo_azul/2, (1 - 0)*largo_azul/2, Nq)  # cargas -q
xn = np.full(Nq, x_azul)
zn = np.zeros(Nq)

# ── Malla 2D (plano z=0, donde están las cargas) ──────────────────────────
xx = np.linspace(-3, 3, 30)
yy = np.linspace(-3, 3, 30)
X, Y = np.meshgrid(xx, yy)
Z    = np.zeros_like(X)

Ex = np.zeros_like(X)
Ey = np.zeros_like(Y)
Ez = np.zeros_like(Z)

for k_i in range(Nq):
    # Carga positiva en placa roja
    rx = X - xp[k_i];  ry = Y - yp[k_i];  rz = Z - zp[k_i]
    r  = np.sqrt(rx**2 + ry**2 + rz**2)
    r[r < 0.15] = np.nan          # evitar singularidades
    Ex += k * q * rx / r**3
    Ey += k * q * ry / r**3
    Ez += k * q * rz / r**3

    # Carga negativa en placa azul
    rx = X - xn[k_i];  ry = Y - yn[k_i];  rz = Z - zn[k_i]
    r  = np.sqrt(rx**2 + ry**2 + rz**2)
    r[r < 0.15] = np.nan
    Ex += k * (-q) * rx / r**3
    Ey += k * (-q) * ry / r**3
    Ez += k * (-q) * rz / r**3

# Normalizar para visualización (solo dirección)
E_mag = np.sqrt(Ex**2 + Ey**2 + Ez**2)
with np.errstate(invalid='ignore'):
    Ex_n = np.where(np.isfinite(E_mag) & (E_mag > 0), Ex/E_mag, 0)
    Ey_n = np.where(np.isfinite(E_mag) & (E_mag > 0), Ey/E_mag, 0)
    Ez_n = np.where(np.isfinite(E_mag) & (E_mag > 0), Ez/E_mag, 0)

# ── Función para construir vértices de un prisma rectangular ──────────────
def prisma_vertices(xc, yc, zc, ancho, largo, v):
    """
    Devuelve (vertices, caras) para un prisma rectangular centrado en (xc,yc,zc).
    Equivalente a crear_prisma de MATLAB.
    """
    a = abs(ancho)
    l = largo / 2
    h = v
    verts = np.array([
        [xc - a, yc - l, zc - h],  # v0
        [xc - a, yc + l, zc - h],  # v1
        [xc + a, yc + l, zc - h],  # v2
        [xc + a, yc - l, zc - h],  # v3
        [xc - a, yc - l, zc + h],  # v4
        [xc - a, yc + l, zc + h],  # v5
        [xc + a, yc + l, zc + h],  # v6
        [xc + a, yc - l, zc + h],  # v7
    ])
    # 6 caras del cubo (índices de vértices)
    faces_idx = [
        [0, 1, 2, 3],  # bottom
        [4, 5, 6, 7],  # top
        [0, 1, 5, 4],  # front
        [2, 3, 7, 6],  # back
        [0, 3, 7, 4],  # left
        [1, 2, 6, 5],  # right
    ]
    faces = [[verts[i] for i in f] for f in faces_idx]
    return faces

# ── Figura 3D ─────────────────────────────────────────────────────────────
fig = plt.figure(figsize=(9, 7))
ax  = fig.add_subplot(111, projection='3d')

# Placa roja (+)
faces_roja = prisma_vertices(x_roja, 0, 0, ancho, largo_roja, v)
poly_roja  = Poly3DCollection(faces_roja, alpha=0.5, facecolor='red', edgecolor='darkred', linewidth=0.4)
ax.add_collection3d(poly_roja)

# Placa azul (-)
faces_azul = prisma_vertices(x_azul, 0, 0, ancho, largo_azul, v)
poly_azul  = Poly3DCollection(faces_azul, alpha=0.5, facecolor='blue', edgecolor='navy', linewidth=0.4)
ax.add_collection3d(poly_azul)

# Posiciones de las cargas (marcadores)
ax.scatter(xp, yp, zp, color='yellow', s=80, zorder=5, depthshade=False)
ax.scatter(xn, yn, zn, color='yellow', s=80, zorder=5, depthshade=False)

# Campo eléctrico (quiver3)
scale = 0.35
ax.quiver(X, Y, Z,
          Ex_n * scale, Ey_n * scale, Ez_n * scale,
          color='#00cc44', linewidth=0.7, arrow_length_ratio=0.3)

# ── Formato ───────────────────────────────────────────────────────────────
ax.set_xlim(-3, 3)
ax.set_ylim(-3, 3)
ax.set_zlim(-3, 3)
ax.view_init(elev=30, azim=30)
ax.set_xlabel(r'$x$ (m)', fontsize=12)
ax.set_ylabel(r'$y$ (m)', fontsize=12)
ax.set_zlabel(r'$z$ (m)', fontsize=12)
ax.set_title('Campo eléctrico: arreglo de cargas en dos placas',
             fontsize=13)
ax.grid(True, linestyle='--', alpha=0.3)

plt.tight_layout()
plt.savefig(os.path.join(carpeta, 'campo_placas_paralelas.png'), dpi=150, bbox_inches='tight')
plt.show()
print("Listo: campo_placas_paralelas.png")