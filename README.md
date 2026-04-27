# Visualización de Campos Eléctricos — Detección de Malaria

Simulaciones y visualizaciones de campos eléctricos en Python/Matplotlib, desarrolladas como parte de un reto de física aplicada. El objetivo final es modelar el principio de **dielectroforesis** para la detección de malaria en glóbulos rojos infectados.

---

## Propósito

Los glóbulos rojos infectados con malaria (*Plasmodium falciparum*) tienen propiedades dieléctricas distintas a los sanos. Al exponerlos a un **campo eléctrico no uniforme**, experimentan fuerzas dielectroforéticas diferentes, lo que permite separarlos y detectar la infección sin métodos invasivos.

Este repositorio modela ese principio paso a paso, desde un dipolo simple hasta la simulación del movimiento de células en un campo generado por electrodos asimétricos.

---

## Estructura

```
Actividades/
├── Dipolo/
│   ├── Dipolo_Electrico.py               # Campo eléctrico 2D de un dipolo de cargas puntuales
│   ├── Documentacion_Dipolo.md           # Documentación detallada
│   └── campo_electrico_dipolo.png        # Imagen generada
│
├── Placas_Cargadas/
│   ├── Placas_Cargadas.py                # Campo 3D de dos placas con distribución discreta de cargas
│   ├── Documentacion_Placas_Cargadas.md  # Documentación detallada
│   └── campo_placas_paralelas.png        # Imagen generada
│
├── Campos_Malaria/
│   ├── Campos_Malaria.py                 # Campo no uniforme de placas asimétricas + potencial
│   ├── Documentacion_Campo_Malaria.md    # Documentación detallada
│   └── campo_placas_malaria.png          # Imagen generada
│
├── Celulas_Simuladas/
│   ├── Simulacion_Malaria.py             # Trayectorias de células sanas vs infectadas
│   ├── Documentacion_Celulas_Simuladas.md# Documentación detallada
│   └── trayectorias_malaria.png          # Imagen generada
│
├── Version_Final_Segmentada/
│   ├── main.py                           # Punto de entrada
│   ├── campoElectrico.py
│   ├── trayectorias.py                   # Genera trayectorias_malaria.csv
│   └── grafico.py
│
└── ModeloIdentificacion/
    ├── ia.py                             # Clasificador ML (PCA + 4 algoritmos)
    ├── Documentacion_Modelo_Identificacion.md
    └── trayectorias_malaria.csv          # CSV generado por Version_Final_Segmentada
```

---

## Actividades

### 1. Dipolo Eléctrico
Campo vectorial 2D de dos cargas puntuales (+q, −q) usando superposición.
- Malla `meshgrid` con `quiver` normalizado
- Exporta PNG

### 2. Placas Cargadas (campo 3D)
Arreglo de cargas distribuidas en dos placas paralelas del mismo tamaño (positiva y negativa).
- Visualización 3D con `Poly3DCollection` (prismas) + `quiver3`
- Principio de superposición sobre `Nq` cargas por placa
- Máscara para no dibujar flechas dentro de las placas

### 3. Campo No Uniforme — Configuración para Malaria
Dos placas asimétricas (larga positiva, corta negativa) generando el gradiente de campo necesario para dielectroforesis.
- Malla de alta resolución (1000×1000) con `Nq = 500` cargas por placa
- Mapa de potencial (`contourf`), líneas equipotenciales y líneas de campo (`streamplot`)
- Exporta PNG sin abrir ventana gráfica (backend `Agg`)

### 4. Simulación de Células — Trayectorias Sanas vs Infectadas
Integración numérica (Euler explícito) de 30 células bajo la fuerza de Coulomb del campo asimétrico.
- Células sanas (`q = 0.8e-6 C`) y células infectadas (`q = 1.4e-6 C`)
- Ruido aleatorio en la carga para mayor realismo
- Trayectorias superpuestas sobre la visualización del campo de fondo

### 5. Versión Final Segmentada
Refactorización del código en módulos independientes para mayor mantenibilidad. Además de la imagen PNG, genera `trayectorias_malaria.csv` con las trayectorias detalladas por partícula.
- `campoElectrico.py` — cálculo del campo y potencial
- `trayectorias.py` — simulación dinámica de partículas + exportación CSV
- `grafico.py` — visualización y exportación
- `main.py` — punto de entrada principal

### 6. Modelo de Identificación — Clasificador ML
Carga el CSV de trayectorias, extrae 8 características observables por partícula, reduce con PCA a 2 componentes y compara cuatro clasificadores:
- **Naive Bayes Gaussiano**, **Árbol de Decisión**, **KNN (k=5)**, **Regresión Logística**
- Genera matrices de confusión individuales (PNG) y un gráfico de barras comparativo de accuracy y recall macro.

---

## Instalación

```bash
git clone https://github.com/luisxavierxd/CamposElectricosMalaria
cd CamposElectricosMalaria
python -m venv venv
venv\Scripts\activate        # Windows
pip install numpy matplotlib seaborn scikit-learn pandas
```

---

## Uso

Correr siempre desde la **raíz del proyecto**:

```bash
python Actividades/Dipolo/Dipolo_Electrico.py
python Actividades/Placas_Cargadas/Placas_Cargadas.py
python Actividades/Campos_Malaria/Campos_Malaria.py
python Actividades/Celulas_Simuladas/Simulacion_Malaria.py

# Versión segmentada (genera PNG + trayectorias_malaria.csv)
python Actividades/Version_Final_Segmentada/main.py

# Modelo de identificación ML (requiere trayectorias_malaria.csv en ModeloIdentificacion/)
python Actividades/ModeloIdentificacion/ia.py
```

Los archivos PNG se guardan automáticamente en la carpeta de cada script.

---

## Referencias

- Pohl, H.A. (1978). *Dielectrophoresis*. Cambridge University Press.
- Gascoyne, P. et al. — Separación de células infectadas con malaria por dielectroforesis.
- [Separación de proteínas con DEP — ITESM](https://repositorio.tec.mx/items/29a9ce9c-7cae-4165-affc-077f2b5c8e41)
- [Dielectrophoresis — Springer Reference](https://link.springer.com/referenceworkentry/10.1007/978-90-481-9751-4_131)

---

## Autores

- Luis Xavier García Pimentel Ascencio
- Mario Donaciano Castillos Santos
- Angel Raúl Luna Tirado
- Fernando Gómez López
- Camila Ruiz Casas
