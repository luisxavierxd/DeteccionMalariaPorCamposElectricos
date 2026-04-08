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
│   └── Dipolo_Electrico.py       # Campo eléctrico de un dipolo de cargas puntuales (2D)
└── Placas_Cargadas/
    └── Placas_Cargadas.py        # Campo 3D de dos placas con distribución de cargas
```

---

## Actividades

### 1. Dipolo Eléctrico
Campo vectorial 2D de dos cargas puntuales (+q, −q) usando superposición.
- Malla `meshgrid` con `quiver` normalizado
- Exporta PDF vectorial

### 2. Placas Cargadas (campo 3D)
Arreglo de cargas distribuidas en dos placas paralelas (positiva y negativa).
- Visualización 3D con `Poly3DCollection` (prismas) + `quiver3`
- Principio de superposición sobre `Nq` cargas por placa
- Máscara para no dibujar flechas dentro de las placas

### 3. Campo No Uniforme *(próximamente)*
Electrodos asimétricos generando gradiente de campo → condición necesaria para dielectroforesis.

### 4. Simulación de Células *(próximamente)*
Trayectoria de glóbulos rojos sanos vs infectados bajo la fuerza dielectroforética:

$$\vec{F}_{DEP} \propto \nabla|\vec{E}|^2$$

---

## Instalación

```bash
git clone https://github.com/luisxavierxd/CamposElectricosMalaria
cd CamposElectricosMalaria
python -m venv venv
venv\Scripts\activate        # Windows
pip install numpy matplotlib
```

---

## Uso

Correr siempre desde la **raíz del proyecto**:

```bash
python Actividades/Dipolo/Dipolo_Electrico.py
python Actividades/Placas_Cargadas/Placas_Cargadas.py
```

Los archivos PDF/PNG se guardan automáticamente en la carpeta de cada script.

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
