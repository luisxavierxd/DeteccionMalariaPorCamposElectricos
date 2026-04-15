import os
import numpy as np
import matplotlib                      # Debe importarse antes de pyplot.
matplotlib.use("Agg")                  # Backend sin GUI (escritura a archivo).
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from matplotlib.patches import Rectangle
from matplotlib.animation import FuncAnimation

carpeta = os.path.dirname(os.path.abspath(__file__))

# ── Parámetros físicos ────────────────────────────────────────────────────
k  = 8.99e9   # Constante de Coulomb (N·m²/C²) - valor realista
q  = 1.6e-19  # Carga elemental (C)
Nq = 3

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
xx = np.linspace(-3, 3, 200)
yy = np.linspace(-3, 3, 200)
X, Y = np.meshgrid(xx, yy)
Z    = np.zeros_like(X)

Ex = np.zeros_like(X)
Ey = np.zeros_like(Y)
Ez = np.zeros_like(Z)
V  = np.zeros_like(X)

print("Calculando campo eléctrico y potencial...")
for k_i in range(Nq):
    # ── Carga positiva (+q) ──
    rx = X - xp[k_i];  ry = Y - yp[k_i];  rz = Z - zp[k_i]
    rp = np.sqrt(rx**2 + ry**2 + rz**2)
    rp[rp < 0.05] = np.nan
    Ex += k * q * rx / rp**3
    Ey += k * q * ry / rp**3
    Ez += k * q * rz / rp**3

    # ── Carga negativa (-q) ──
    rx = X - xn[k_i];  ry = Y - yn[k_i];  rz = Z - zn[k_i]
    rn = np.sqrt(rx**2 + ry**2 + rz**2)
    rn[rn < 0.05] = np.nan
    Ex += k * (-q) * rx / rn**3
    Ey += k * (-q) * ry / rn**3
    Ez += k * (-q) * rz / rn**3

    # ── Potencial eléctrico V = k*q/rp - k*q/rn ──────────────────────────
    rp2 = np.sqrt((X - xp[k_i])**2 + (Y - yp[k_i])**2)
    rp2[rp2 < 0.05] = 0.05  # evitar singularidades
    rn2 = np.sqrt((X - xn[k_i])**2 + (Y - yn[k_i])**2)
    rn2[rn2 < 0.05] = 0.05  # evitar singularidades
    V += k * q / rp2 - k * q / rn2

# ── Gradiente del campo eléctrico (para fuerza dielectroforética) ─────────
# Calculamos |E|² y su gradiente
E_magnitude_sq = Ex**2 + Ey**2 + Ez**2
grad_E_magnitude_sq_y, grad_E_magnitude_sq_x = np.gradient(E_magnitude_sq, yy, xx)

# ============================================================================
# SIMULACIÓN DE CÉLULAS: SANAS VS INFECTADAS (basado en diapositivas)
# ============================================================================

# ── Parámetros de simulación de trayectorias ──────────────────────────────
n_celulas = 50           # Número de células a simular
dt = 1e-5                # Paso de tiempo (s)
b = 0.4e-9               # Coeficiente de arrastre (ajustado para escala)
ke = k                  # Constante de Coulomb para la simulación

# Posiciones iniciales de las células (caen desde arriba)
x_inicial = 0.0
y_inicial = 2.5

# Almacenar resultados
resultados = []  # (clase, x_final, tiempo_vuelo)

# Función para simular el dipolo (célula con separación de cargas dx)
def simular_celula(clase, dx, qE, masa, dt, max_steps=5000):
    """
    Simula una célula como un dipolo con separación dx entre cargas.
    Calcula la fuerza neta usando Coulomb + arrastre.
    
    clase: 0 = sana, 1 = infectada
    dx: separación entre cargas (momento dipolar)
    qE: carga efectiva
    masa: masa de la célula
    """
    xe = x_inicial
    ye = y_inicial
    Vx = 0.0
    Vy = 0.0
    
    trayectoria_x = [xe]
    trayectoria_y = [ye]
    
    t = 0
    
    while ye > -2.5 and ye < 3 and xe > -2.5 and xe < 2.5 and t < max_steps * dt:
        Fx = 0
        Fy = 0
        
        # Simulación del dipolo: cargas en xe ± dx
        # Posición de las dos cargas del dipolo
        x_neg = xe - dx   # carga negativa del dipolo
        x_pos = xe + dx   # carga positiva del dipolo
        
        # Interacción con la placa negativa (cargas -q en xn, yn)
        for k_i in range(Nq):
            # Para la carga negativa del dipolo (en x_neg)
            rx_neg = x_neg - xn[k_i]
            ry_neg = ye - yn[k_i]
            r_neg = np.sqrt(rx_neg**2 + ry_neg**2)
            if r_neg > 0.01:
                # Fuerza sobre carga negativa del dipolo (atraída por placa negativa?)
                # qE es la carga efectiva de la célula
                Fx += -ke * qE * (-q) * rx_neg / r_neg**3
                Fy += -ke * qE * (-q) * ry_neg / r_neg**3
            
            # Para la carga positiva del dipolo (en x_pos)
            rx_pos = x_pos - xn[k_i]
            ry_pos = ye - yn[k_i]
            r_pos = np.sqrt(rx_pos**2 + ry_pos**2)
            if r_pos > 0.01:
                Fx += -ke * qE * q * rx_pos / r_pos**3
                Fy += -ke * qE * q * ry_pos / r_pos**3
        
        # Fuerza de arrastre (proporcional a la velocidad)
        F_arrastre_x = -b * Vx
        F_arrastre_y = -b * Vy
        
        # Segunda ley de Newton: F_total = m * a
        ax = (Fx + F_arrastre_x) / masa
        ay = (Fy + F_arrastre_y) / masa
        
        # Euler
        Vx += ax * dt
        Vy += ay * dt
        xe += Vx * dt
        ye += Vy * dt
        
        trayectoria_x.append(xe)
        trayectoria_y.append(ye)
        t += dt
    
    return np.array(trayectoria_x), np.array(trayectoria_y), t

print("\n" + "="*60)
print("SIMULACIÓN DE CÉLULAS SANAS vs INFECTADAS")
print("="*60)
print("Basado en diapositivas:")
print("  - Célula sana (dx pequeño): momento dipolar pequeño → poca desviación")
print("  - Célula infectada (dx grande): momento dipolar grande → mayor desviación")
print()

# Simular células sanas e infectadas
trayectorias_sanas = []
trayectorias_infectadas = []

for i in range(n_celulas):
    # Determinar si es sana o infectada (50% probabilidad)
    if np.random.rand() > 0.5:
        # CÉLULA SANA (dx pequeño)
        clase = 0
        dx = 0.05 + np.random.rand() * 0.05      # Separación pequeña (dipolo débil)
        qE = 0.8e-6 + np.random.randn() * 0.25e-6  # Carga efectiva base
        masa = 1.0
        color = 'green'
        label = 'Sana'
    else:
        # CÉLULA INFECTADA (dx grande)
        clase = 1
        vol_parasito = 0.1 + np.random.rand() * 0.8  # Volumen del parásito (10-90%)
        dx = 0.25 + vol_parasito * 0.3               # Separación mayor (dipolo fuerte)
        qE = 1.4e-6 + np.random.randn() * 0.25e-6    # Carga efectiva mayor (50-100% más)
        masa = 1.0 + vol_parasito * 0.2              # Mayor masa por parásito
        color = 'red'
        label = 'Infectada'
    
    # Simular trayectoria
    x_traj, y_traj, tiempo = simular_celula(clase, dx, qE, masa, dt)
    
    # Guardar resultado
    x_final = x_traj[-1]
    resultados.append({
        'clase': clase,
        'label': label,
        'x_final': x_final,
        'tiempo': tiempo,
        'dx': dx,
        'qE': qE,
        'masa': masa
    })
    
    if clase == 0:
        trayectorias_sanas.append((x_traj, y_traj, color))
    else:
        trayectorias_infectadas.append((x_traj, y_traj, color))
    
    print(f"Célula {i+1:2d}: {label:10s} | dx = {dx:.4f} m | qE = {qE:.2e} C | "
          f"x_final = {x_final:6.3f} m | tiempo = {tiempo:.4f} s")

# ============================================================================
# GRÁFICA 1: Campo eléctrico y potencial (como el original)
# ============================================================================
fig2, ax2 = plt.subplots(figsize=(9, 7))

V_plot = np.clip(V, -500, 500)

# contourf con colormap 'bone'
cf = ax2.contourf(X, Y, V_plot, levels=60, cmap='bone')
plt.colorbar(cf, ax=ax2, label='Potencial eléctrico V')

# Líneas equipotenciales blancas
ax2.contour(X, Y, V_plot, levels=10, colors='white', linewidths=0.5)

# streamplot
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

# ============================================================================
# GRÁFICA 2: Trayectorias de células sanas vs infectadas
# ============================================================================
fig3, ax3 = plt.subplots(figsize=(10, 8))

# Dibujar las placas
ax3.add_patch(Rectangle((x_roja-abs(ancho), -largo_roja/2), 2*abs(ancho), largo_roja,
                         color='red', alpha=0.3, label='Placa +'))
ax3.add_patch(Rectangle((x_azul-abs(ancho), -largo_azul/2), 2*abs(ancho), largo_azul,
                         color='blue', alpha=0.3, label='Placa −'))

# Trayectorias de células sanas (verde)
for x_traj, y_traj, color in trayectorias_sanas:
    ax3.plot(x_traj, y_traj, color='green', alpha=0.5, linewidth=0.8)
    ax3.scatter(x_traj[0], y_traj[0], color='green', s=30, marker='o', alpha=0.7)
    ax3.scatter(x_traj[-1], y_traj[-1], color='darkgreen', s=50, marker='s', alpha=0.8)

# Trayectorias de células infectadas (rojo)
for x_traj, y_traj, color in trayectorias_infectadas:
    ax3.plot(x_traj, y_traj, color='red', alpha=0.5, linewidth=0.8)
    ax3.scatter(x_traj[0], y_traj[0], color='red', s=30, marker='o', alpha=0.7)
    ax3.scatter(x_traj[-1], y_traj[-1], color='darkred', s=50, marker='s', alpha=0.8)

ax3.set_xlim(-2.5, 2.5)
ax3.set_ylim(-2.5, 3)
ax3.set_xlabel(r'$x$ (m)')
ax3.set_ylabel(r'$y$ (m)')
ax3.set_title('Trayectorias de células: SANAS (verde) vs INFECTADAS (rojo)\n'
              'Las infectadas se desvían lateralmente por mayor momento dipolar', 
              fontsize=12)
ax3.legend(loc='upper right')
ax3.grid(True, linestyle='--', alpha=0.3)
ax3.set_aspect('equal')

# Añadir anotación explicativa
ax3.annotate('Campo más fuerte\n(∇E mayor)', xy=(0.5, 0), xytext=(1.2, 0.5),
             arrowprops=dict(arrowstyle='->', color='gray'),
             fontsize=9, color='gray')
ax3.annotate('Las infectadas se desvían\nhacia la placa negativa', 
             xy=(1.0, -1), xytext=(1.5, -1.8),
             arrowprops=dict(arrowstyle='->', color='red'),
             fontsize=9, color='red')

plt.tight_layout()
plt.savefig(os.path.join(carpeta, 'trayectorias_celulas_malaria.png'), dpi=150, bbox_inches='tight')

# ============================================================================
# GRÁFICA 3: Distribución de posiciones finales
# ============================================================================
fig4, ax4 = plt.subplots(figsize=(8, 5))

x_final_sanas = [r['x_final'] for r in resultados if r['clase'] == 0]
x_final_infectadas = [r['x_final'] for r in resultados if r['clase'] == 1]

ax4.hist(x_final_sanas, bins=15, alpha=0.7, color='green', label='Sanas', density=True)
ax4.hist(x_final_infectadas, bins=15, alpha=0.7, color='red', label='Infectadas', density=True)

ax4.axvline(x=x_azul, color='blue', linestyle='--', linewidth=2, label='Placa negativa')
ax4.axvline(x=x_roja, color='red', linestyle='--', linewidth=2, label='Placa positiva')

ax4.set_xlabel('Posición final en x (m)')
ax4.set_ylabel('Densidad de probabilidad')
ax4.set_title('Distribución de posiciones finales\n'
              'Las infectadas tienden a desviarse más hacia la placa negativa',
              fontsize=12)
ax4.legend()
ax4.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig(os.path.join(carpeta, 'distribucion_celulas_malaria.png'), dpi=150, bbox_inches='tight')

# ============================================================================
# IMPRIMIR RESUMEN
# ============================================================================
print("\n" + "="*60)
print("RESUMEN DE LA SIMULACIÓN")
print("="*60)
print(f"Total células simuladas: {n_celulas}")
print(f"Células sanas: {len(x_final_sanas)}")
print(f"Células infectadas: {len(x_final_infectadas)}")
print(f"\nPosición media final - Sanas: {np.mean(x_final_sanas):.4f} m")
print(f"Posición media final - Infectadas: {np.mean(x_final_infectadas):.4f} m")
print(f"Diferencia: {abs(np.mean(x_final_sanas) - np.mean(x_final_infectadas)):.4f} m")

print("\n✅ Gráficos guardados:")
print("   - campo_placas_malaria.png")
print("   - trayectorias_celulas_malaria.png")
print("   - distribucion_celulas_malaria.png")

plt.show()
print("\nListo: Simulación completada")