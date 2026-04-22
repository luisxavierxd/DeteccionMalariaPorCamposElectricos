# Modelo de Identificación — Clasificación de Células Sanas vs Infectadas

Clasificador de aprendizaje automático que predice si una célula es sana o infectada con *Plasmodium falciparum* a partir de las características geométricas de su trayectoria bajo el campo eléctrico no uniforme. Compara cuatro algoritmos de clasificación sobre las componentes principales obtenidas con PCA.

---

## Descripción

El script carga el CSV generado por `Version_Final_Segmentada/main.py`, extrae ocho características observables por partícula —incluyendo la masa, que varía entre células infectadas y sanas—, aplica PCA para encontrar las dos direcciones de máxima varianza y entrena cuatro clasificadores:

1. **Naive Bayes Gaussiano** — modelo probabilístico que asume independencia entre componentes
2. **Árbol de Decisión** — partición recursiva del espacio de componentes por umbrales
3. **KNN (k=5)** — asigna la clase mayoritaria entre los 5 vecinos más cercanos
4. **Regresión Logística** — frontera de decisión lineal en el espacio PC

Para cada modelo imprime en consola la matriz de confusión y el reporte completo de métricas, guarda la matriz como PNG y al final genera un gráfico de barras comparando accuracy y recall macro de los cuatro modelos.

---

## Ejemplo de salida

```
confusion_naive_bayes.png
confusion_arbol.png
confusion_knn.png
confusion_logistica.png
comparacion_modelos.png
```

La consola muestra, para cada modelo, la matriz de confusión en texto, precisión, recall y F1 por clase, y un resumen tabular final con accuracy y recall macro de los cuatro modelos.

---

## Estructura del código

```
ia.py
│
├── Carga del CSV                  # trayectorias_malaria.csv → DataFrame
├── Ingeniería de características  # 7 features por partícula desde xs, ys
├── Escalado MinMax                # Normaliza features a [0, 1]
├── División train / test          # 80 % / 20 % estratificado
├── PCA                            # Reduce a 2 componentes principales (PC1, PC2)
├── Naive Bayes                    # fit → predict → matriz + reporte
├── Árbol de Decisión              # fit → predict → matriz + reporte
├── KNN (k=5)                      # fit → predict → matriz + reporte
├── Regresión Logística            # fit → predict → matriz + reporte
├── Resumen en consola             # Tabla accuracy / recall de los 4 modelos
└── Gráfico de barras              # comparacion_modelos.png
```

---

## Parámetros configurables

### Extracción de características

| Feature        | Descripción |
|---|---|
| `x_inicial`    | Posición X de entrada de la célula (m) |
| `x_final`      | Posición X al salir de la región (m) |
| `y_final`      | Posición Y al salir de la región (m) |
| `n_pasos`      | Número de pasos de integración (duración de la trayectoria) |
| `dx_neto`      | Desplazamiento lateral neto: `x_final − x_inicial` (m) |
| `x_max_abs`    | Máxima desviación lateral absoluta durante la trayectoria (m) |
| `x_std`        | Desviación estándar de las posiciones X (variabilidad lateral) |
| `masa`         | Masa de la célula (adimensional; solo varía en infectadas) |

### PCA y modelos

| Parámetro        | Valor     | Descripción |
|---|---|---|
| `N_COMPONENTES`  | `2`       | Componentes principales retenidos |
| `test_size`      | `0.2`     | Fracción del dataset para prueba |
| `random_state`   | `0`       | Semilla de reproducibilidad |
| `max_depth`      | `3`       | Profundidad máxima del árbol de decisión |
| `n_neighbors`    | `5`       | Vecinos en KNN |
| `max_iter`       | `500`     | Iteraciones máximas en Regresión Logística |

---

## Fundamento teórico

### PCA vs LDA — por qué se eligió PCA

Ambas técnicas reducen la dimensionalidad del espacio de features, pero con objetivos distintos:

| | PCA | LDA |
|---|---|---|
| **Objetivo** | Maximizar varianza total de los datos | Maximizar separación entre clases |
| **Supervisado** | No (ignora etiquetas) | Sí (usa las etiquetas durante el ajuste) |
| **Componentes** | Hasta `n_features` | Hasta `n_clases − 1` |
| **Útil cuando** | Hay múltiples fuentes de variación independientes | Solo importa la frontera entre clases |

LDA fue descartado porque con 2 clases produce un único componente discriminante (LD1), lo que impide capturar los **dos ejes físicos independientes** presentes en este problema: la deflexión lateral (determinada por la carga) y la dinámica temporal (determinada por la masa). Al forzar todo a una dimensión, LDA pierde la información complementaria de la masa.

PCA, en cambio, puede retener **2 componentes** que corresponden a estas dos fuentes de variación ortogonales:
- **PC1** — captura principalmente la variación en deflexión lateral (`dx_neto`, `x_max_abs`, `x_std`, `x_final`), relacionada con la carga
- **PC2** — captura principalmente la variación temporal (`n_pasos`, `y_final`, `masa`), relacionada con la inercia de la célula

### Análisis de Componentes Principales (PCA)

PCA encuentra las $k$ direcciones ortogonales de máxima varianza en el espacio de features:

$$\mathbf{Z} = \mathbf{X} \mathbf{W}$$

donde $\mathbf{W}$ contiene los $k$ autovectores principales de la matriz de covarianza de $\mathbf{X}$, ordenados por varianza explicada descendente. Con $k=2$ se retienen las dos direcciones que concentran más información, eliminando redundancia entre las 8 features correlacionadas.

### Naive Bayes Gaussiano

Aplica el teorema de Bayes asumiendo que LD1 sigue una distribución normal por clase:

$$P(clase \mid z) \propto P(clase) \cdot P(z \mid clase)$$

### Árbol de Decisión

Encuentra un umbral sobre LD1 que minimiza la impureza de Gini:

$$G = 1 - \sum_{c} p_c^2$$

### KNN (k vecinos más cercanos)

Asigna a cada punto de prueba la clase mayoritaria entre sus $k$ vecinos más cercanos sobre el eje LD1.

### Regresión Logística

Modela la probabilidad con una sigmoide sobre LD1:

$$P(Infectada \mid z) = \frac{1}{1 + e^{-(wz + b)}}$$

---

## Modelo de detección

Las células infectadas difieren de las sanas en **dos propiedades físicas simultáneas**:

1. **Carga mayor** (`q_infectada = 1.4e-6 C` vs `q_sana = 0.8e-6 C`): produce mayor deflexión lateral al atravesar el campo no uniforme, reflejada en `dx_neto`, `x_max_abs` y `x_std`.
2. **Masa mayor** (`masa = 1 + U[0.6, 1.0]` vs `masa = 1`): la mayor inercia hace que las infectadas aceleren más despacio, tomando más pasos y saliendo en posiciones `y_final` distintas.

Estas dos fuentes de variación son **ortogonales entre sí** — una célula puede tener mucha carga pero poca masa extra, o viceversa — lo que justifica usar PCA con 2 componentes: PC1 captura la deflexión (carga) y PC2 captura la dinámica (masa). Los cuatro clasificadores reciben ambas componentes y pueden explotar las dos señales físicas simultáneamente.

---

## Dependencias

```bash
pip install numpy pandas matplotlib seaborn scikit-learn
```

| Librería       | Versión mínima | Uso |
|---|---|---|
| `numpy`        | 1.20+          | Operaciones vectoriales sobre trayectorias |
| `pandas`       | 1.3+           | Carga de CSV y construcción de features |
| `matplotlib`   | 3.3+           | Gráfico de barras comparativo |
| `seaborn`      | 0.11+          | Heatmaps de matrices de confusión |
| `scikit-learn` | 1.0+           | LDA, escalado, modelos y métricas |

---

## Ejecución

Primero generar el CSV (si no existe o se quiere actualizar):

```bash
python Actividades/Version_Final_Segmentada/main.py
```

Luego correr el clasificador:

```bash
python Actividades/ModeloIdentificacion/ia.py
```

Los archivos de salida se guardan automáticamente en la misma carpeta:

```
Actividades/ModeloIdentificacion/confusion_naive_bayes.png
Actividades/ModeloIdentificacion/confusion_arbol.png
Actividades/ModeloIdentificacion/confusion_knn.png
Actividades/ModeloIdentificacion/confusion_logistica.png
Actividades/ModeloIdentificacion/comparacion_modelos.png
```

No se abre ventana gráfica (usa el backend `Agg` de matplotlib).
