# Simulación de Dielectroforesis — Detección de Malaria (modelo físico + ML)

Simulaciones y visualizaciones de campos eléctricos en Python/Matplotlib que **modelan** el principio de **dielectroforesis** aplicado a la separación de glóbulos rojos infectados con malaria. Proyecto de física aplicada (Sistemas Eléctricos).

> **Nota importante:** este proyecto es una **simulación y prueba de concepto computacional**, no un sistema de detección clínica real. Modela cómo *se comportarían* células sanas vs. infectadas bajo un campo eléctrico no uniforme, y clasifica sus trayectorias simuladas con machine learning. No procesa muestras biológicas reales.

---

## Propósito

Los glóbulos rojos infectados con *Plasmodium falciparum* tienen propiedades dieléctricas distintas a los sanos. Bajo un **campo eléctrico no uniforme**, ambos experimentan fuerzas dielectroforéticas diferentes, lo que en principio permitiría separarlos. Este repositorio modela ese principio paso a paso — de un dipolo simple hasta la simulación del movimiento de células bajo un campo generado por electrodos asimétricos — y luego entrena un clasificador que distingue sanas de infectadas a partir de las trayectorias simuladas.

---

## Estructura

```
Actividades/
├── Dipolo/                     # Campo E 2D de un dipolo de cargas puntuales
├── Placas_Cargadas/            # Campo 3D de dos placas con distribución discreta de cargas
├── Campos_Malaria/             # Campo no uniforme de placas asimétricas + potencial
├── Celulas_Simuladas/          # Trayectorias de células sanas vs. infectadas
├── Version_Final_Segmentada/   # Refactor modular; genera trayectorias_malaria.csv
└── ModeloIdentificacion/       # Clasificador ML (PCA + 4 algoritmos) sobre el CSV
```

---

## Actividades

1. **Dipolo eléctrico** — campo vectorial 2D de dos cargas puntuales (+q, −q) por superposición; `meshgrid` + `quiver`.
2. **Placas cargadas (3D)** — arreglo de cargas en dos placas paralelas; visualización 3D con `Poly3DCollection` + `quiver3` y superposición sobre `Nq` cargas por placa.
3. **Campo no uniforme (config. malaria)** — dos placas asimétricas (larga positiva, corta negativa) generan el gradiente necesario para dielectroforesis; malla de alta resolución, mapa de potencial (`contourf`), equipotenciales y líneas de campo (`streamplot`).
4. **Simulación de células** — integración numérica (Euler explícito) de 30 células bajo la fuerza de Coulomb del campo asimétrico; células sanas (`q = 0.8e-6 C`) vs. infectadas (`q = 1.4e-6 C`), con ruido en la carga para realismo.
5. **Versión final segmentada** — código refactorizado en módulos (`campoElectrico.py`, `trayectorias.py`, `grafico.py`, `main.py`); exporta `trayectorias_malaria.csv`.
6. **Modelo de identificación (ML)** — carga el CSV, extrae 8 características por partícula, reduce con **PCA** a 2 componentes y compara cuatro clasificadores: **Naive Bayes Gaussiano**, **Árbol de Decisión**, **KNN (k=5)** y **Regresión Logística**. Genera matrices de confusión y un comparativo de accuracy y recall macro.

---

## Instalación

```bash
git clone https://github.com/luisxavierxd/DeteccionMalariaPorCamposElectricos
cd DeteccionMalariaPorCamposElectricos
python -m venv venv
venv\Scripts\activate        # Windows  (macOS/Linux: source venv/bin/activate)
pip install numpy matplotlib seaborn scikit-learn pandas
```

## Uso

Correr siempre desde la **raíz del proyecto**:

```bash
python Actividades/Dipolo/Dipolo_Electrico.py
python Actividades/Placas_Cargadas/Placas_Cargadas.py
python Actividades/Campos_Malaria/Campos_Malaria.py
python Actividades/Celulas_Simuladas/Simulacion_Malaria.py

# Versión segmentada (genera PNG + trayectorias_malaria.csv)
python Actividades/Version_Final_Segmentada/main.py

# Modelo ML (requiere trayectorias_malaria.csv en ModeloIdentificacion/)
python Actividades/ModeloIdentificacion/ia.py
```

Los PNG se guardan en la carpeta de cada script.

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

## Licencia

Distribuido bajo la licencia incluida en el archivo `LICENSE`.
