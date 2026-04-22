# Versión Final Segmentada — Simulación Modular con Exportación CSV

Refactorización modular de `Celulas_Simuladas` que divide la simulación en cuatro módulos independientes y añade la exportación de todas las trayectorias a un archivo CSV, facilitando el análisis posterior de los datos.

---

## Descripción

El código monolítico de `Simulacion_Malaria.py` se reorganiza en cuatro archivos con responsabilidades separadas:

1. **`campoElectrico.py`**: define la geometría de las placas y calcula el campo eléctrico y el potencial sobre la malla de visualización.
2. **`trayectorias.py`**: integra numéricamente las ecuaciones de movimiento de cada célula y exporta los datos a CSV.
3. **`grafico.py`**: genera la figura con el campo de fondo y las trayectorias superpuestas.
4. **`main.py`**: punto de entrada que orquesta los tres módulos anteriores.

La física y los parámetros son idénticos a `Celulas_Simuladas`; únicamente cambia la organización del código y se añade la salida en CSV.

---

## Ejemplo de salida

```
trayectorias_malaria.png
trayectorias_malaria.csv
```

La imagen replica la visualización de `Celulas_Simuladas`: potencial en escala `bone`, líneas equipotenciales blancas, líneas de campo verdes y trayectorias coloreadas (**verde** para células sanas, **rojo** para infectadas).

El CSV contiene una fila por punto de trayectoria con las columnas `particula_id`, `paso`, `x`, `y`, `tipo` y `carga`, lo que permite reproducir o analizar cualquier trayectoria sin re-ejecutar la simulación.

---

## Estructura del código

```
Version_Final_Segmentada/
│
├── main.py                  # Orquestador: campo → trayectorias → gráfico
│
├── campoElectrico.py
│   ├── Constantes de geometría   # DX_PLACA, LARGO_ROJA, LARGO_AZUL, ANCHO, NQ_VIS
│   └── calcular_campo()          # Malla 1000×1000, Ex/Ey/V normalizados (k=1, q=1)
│
├── trayectorias.py
│   ├── Constantes físicas        # K_FIS, Q_FIS, NQ_SIM, Q_SANA, Q_INFECTADA, RUIDO
│   ├── _posiciones_placa()       # Puntos discretos de cada placa (NQ_SIM puntos)
│   ├── simular()                 # Euler explícito, retorna lista de dicts
│   └── exportar_csv()            # Escribe CSV con todas las trayectorias
│
└── grafico.py
    └── graficar()                # contourf + streamplot + trayectorias → PNG
```

---

## Parámetros configurables

### Geometría de las placas (`campoElectrico.py`)

| Parámetro    | Valor  | Descripción |
|---|---|---|
| `DX_PLACA`   | `-1.0` | Posición X de la placa positiva (m) |
| `LARGO_ROJA` | `4`    | Longitud de la placa positiva (m) |
| `LARGO_AZUL` | `1`    | Longitud de la placa negativa (m) |
| `NQ_VIS`     | `500`  | Cargas por placa para el campo de fondo |

### Dinámica de partículas (`trayectorias.py`)

| Parámetro       | Valor    | Descripción |
|---|---|---|
| `K_FIS`         | `9e9`    | Constante de Coulomb (N·m²/C²) |
| `Q_FIS`         | `1e-5`   | Carga de cada punto discreto de la placa (C) |
| `NQ_SIM`        | `150`    | Cargas por placa en la simulación dinámica |
| `N_PARTICULAS`  | `30`     | Número de células simuladas |
| `DT`            | `0.02`   | Paso de tiempo del integrador (s) |
| `MASA`          | `1`      | Masa de cada célula (kg, adimensional) |
| `MAX_STEPS`     | `2000`   | Pasos máximos por trayectoria |
| `V_INI`         | `10`     | Velocidad inicial hacia abajo (m/s) |
| `Q_SANA`        | `0.8e-6` | Carga efectiva de célula sana (C) |
| `Q_INFECTADA`   | `1.4e-6` | Carga efectiva de célula infectada (C) |
| `RUIDO`         | `0.3e-6` | Desviación estándar del ruido en la carga (C) |

---

## Fundamento físico

### Fuerza de Coulomb sobre la célula

La fuerza que ejerce cada carga discreta $q_i$ de la placa sobre una célula de carga $q_{cell}$ en la posición $\vec{r}$ es:

$$\vec{F}_i = k \, q_{cell} \, q_i \frac{\vec{r} - \vec{r}_i}{|\vec{r} - \vec{r}_i|^3}$$

La fuerza total se obtiene por superposición sobre todas las cargas de ambas placas.

### Integración de movimiento (Euler explícito)

$$\vec{v}^{n+1} = \vec{v}^n + \frac{\vec{F}^n}{m} \, \Delta t$$
$$\vec{r}^{n+1} = \vec{r}^n + \vec{v}^{n+1} \, \Delta t$$

### Dielectroforesis (principio subyacente)

La dielectroforesis es la migración de partículas polarizables bajo un campo eléctrico **no uniforme**:

$$\vec{F}_{DEP} \propto \nabla |\vec{E}|^2$$

Los glóbulos rojos infectados con *Plasmodium falciparum* tienen una polarizabilidad eléctrica distinta a los sanos, representada en este modelo como una carga efectiva mayor. Esto provoca una desviación diferencial al atravesar el campo no uniforme, base del principio de separación y detección.

---

## Formato del CSV exportado

El archivo `trayectorias_malaria.csv` tiene una fila por punto de trayectoria:

| Columna        | Tipo    | Descripción |
|---|---|---|
| `particula_id` | entero  | Índice de la célula (0 a N_PARTICULAS−1) |
| `paso`         | entero  | Número de paso de tiempo dentro de la trayectoria |
| `x`            | float   | Posición horizontal de la célula (m) |
| `y`            | float   | Posición vertical de la célula (m) |
| `tipo`         | string  | `"Sana"` o `"Infectada"` |
| `carga`        | float   | Carga efectiva real de esa célula (C, incluye ruido) |

---

## Modelo físico

Las células se introducen desde `y = 2.5` con velocidad inicial uniforme hacia abajo (`vy = -V_INI`). Al entrar en la región de campo no uniforme entre las placas asimétricas, la fuerza de Coulomb actúa de forma diferente sobre células sanas e infectadas:

- **Células sanas** (carga menor, verde): experimentan menor fuerza lateral y atraviesan la región con menor desviación.
- **Células infectadas** (carga mayor, rojo): experimentan mayor fuerza y se desvían más notoriamente hacia una de las placas.

La variación aleatoria con `RUIDO` introduce dispersión realista entre células de la misma categoría. La condición de salida (`y < -2.5` o `|x| > 3`) detiene la integración cuando la célula abandona la región de interés.

---

## Dependencias

```bash
pip install numpy matplotlib
```

| Librería     | Versión mínima | Uso |
|---|---|---|
| `numpy`      | 1.20+          | Álgebra vectorial, mallas, integración numérica |
| `matplotlib` | 3.3+           | Visualización 2D (contourf, streamplot, plot) |

No se requiere ninguna librería adicional para la exportación CSV (módulo `csv` de la biblioteca estándar de Python).

---

## Ejecución

```bash
python Actividades/Version_Final_Segmentada/main.py
```

Los archivos de salida se guardan automáticamente en la misma carpeta del script:

```
Actividades/Version_Final_Segmentada/trayectorias_malaria.png
Actividades/Version_Final_Segmentada/trayectorias_malaria.csv
```

No se abre ventana gráfica (usa el backend `Agg` de matplotlib).
