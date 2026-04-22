# -*- coding: utf-8 -*-
"""Clasificación de células sanas vs infectadas a partir de trayectorias.

Carga el CSV generado por Version_Final_Segmentada/main.py, extrae
características por partícula y compara cuatro clasificadores siguiendo
la misma estructura del clasificador original de géneros de películas.
"""

import os
import sys
sys.stdout.reconfigure(encoding='utf-8')
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import MinMaxScaler
from sklearn.naive_bayes import GaussianNB
from sklearn.tree import DecisionTreeClassifier
from sklearn.neighbors import KNeighborsClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.decomposition import PCA
from sklearn.metrics import accuracy_score, confusion_matrix, recall_score, classification_report

carpeta  = os.path.dirname(os.path.abspath(__file__))
ruta_csv = os.path.join(carpeta, 'trayectorias_malaria.csv')

# ── Cargar CSV ────────────────────────────────────────────────────────────────
df = pd.read_csv(ruta_csv)

# ── Ingeniería de características (una fila por partícula) ────────────────────
# Las variables de entrada son solo datos observables de la trayectoria;
# la carga real no se incluye porque en un experimento real es desconocida.
filas = []
for pid, grupo in df.groupby('particula_id'):
    xs = grupo['x'].values
    ys = grupo['y'].values
    filas.append({
        'x_inicial':  xs[0],
        'x_final':    xs[-1],
        'y_final':    ys[-1],
        'n_pasos':    len(xs),
        'dx_neto':    xs[-1] - xs[0],
        'x_max_abs':  np.max(np.abs(xs)),
        'x_std':      np.std(xs),
        'masa':       grupo['masa'].iloc[0],
        'tipo':       grupo['tipo'].iloc[0],
    })

DF = pd.DataFrame(filas)
DF.head()

# ── Separar X e y ─────────────────────────────────────────────────────────────
X = DF.drop(columns=['tipo'])
y = DF['tipo']

y.value_counts()

# ── Escalado MinMax ───────────────────────────────────────────────────────────
scaler    = MinMaxScaler()
X_scaled  = pd.DataFrame(scaler.fit_transform(X), columns=X.columns)
DF_scaled = pd.concat([X_scaled, y.reset_index(drop=True)], axis=1)

# ── División entrenamiento / prueba (70 / 30) ─────────────────────────────────
X_train, X_test, y_train, y_test = train_test_split(
    X_scaled, y,
    test_size=0.2,
    stratify=y,
    random_state=0,
    shuffle=True,
)

X_train = X_train.reset_index(drop=True)
X_test  = X_test.reset_index(drop=True)
y_train = pd.Series(y_train).reset_index(drop=True)
y_test  = pd.Series(y_test).reset_index(drop=True)

# ── Análisis de Componentes Principales (PCA) ─────────────────────────────────
# Con carga Y masa como ejes físicos independientes, PCA captura ambas fuentes
# de variación en componentes ortogonales.
N_COMPONENTES = 2
pca     = PCA(n_components=N_COMPONENTES)
cols_pc = [f'PC{i+1}' for i in range(N_COMPONENTES)]
X_train = pd.DataFrame(pca.fit_transform(X_train), columns=cols_pc)
X_test  = pd.DataFrame(pca.transform(X_test),      columns=cols_pc)

print("\nVarianza explicada por PCA:")
for i, v in enumerate(pca.explained_variance_ratio_):
    print(f"  PC{i+1}: {v:.1%}")
print(f"  Total acumulado: {pca.explained_variance_ratio_.sum():.1%}")

CLASES = sorted(y.unique())


def _imprimir_confusion(m, nombre):
    ancho = 12
    print(f"\n  Matriz de confusion - {nombre}")
    print(f"  {'':>{ancho}} " + "  ".join(f"{c:>{ancho}}" for c in CLASES) + "  <- Prediccion")
    for i, clase_real in enumerate(CLASES):
        fila = "  ".join(f"{m[i, j]:>{ancho}}" for j in range(len(CLASES)))
        print(f"  {clase_real:>{ancho}} {fila}")


def _graficar_confusion(m, titulo, nombre_archivo):
    fig, ax = plt.subplots(figsize=(5, 4))
    sns.heatmap(m, annot=True, cmap='Greens', fmt='d',
                xticklabels=CLASES, yticklabels=CLASES, ax=ax)
    ax.set_xlabel('Predicción')
    ax.set_ylabel('Real')
    ax.set_title(titulo)
    plt.tight_layout()
    ruta = os.path.join(carpeta, nombre_archivo)
    plt.savefig(ruta, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  Guardada: {ruta}")


resultados        = {}
recall_resultados = {}

# ── Algoritmo de Clasificación de Bayes ───────────────────────────────────────
modelo_nb = GaussianNB()
modelo_nb.fit(X_train, y_train)
y_pred_nb = modelo_nb.predict(X_test)

acc_nb    = accuracy_score(y_test, y_pred_nb)
rec_nb    = recall_score(y_test, y_pred_nb, average='macro')
resultados['Naive Bayes']        = acc_nb
recall_resultados['Naive Bayes'] = rec_nb

m_nb = confusion_matrix(y_test, y_pred_nb, labels=CLASES)
print("\n" + "="*60)
print("NAIVE BAYES")
print("="*60)
_imprimir_confusion(m_nb, 'Naive Bayes')
print(classification_report(y_test, y_pred_nb, target_names=CLASES))
_graficar_confusion(m_nb, f'Naive Bayes (exactitud {acc_nb:.0%})', 'confusion_naive_bayes.png')

# ── Árbol de Decisiones ───────────────────────────────────────────────────────
arbol = DecisionTreeClassifier(max_depth=3, random_state=0)
arbol.fit(X_train, y_train)
var_predecidas = arbol.predict(X_test)

acc_arbol = accuracy_score(y_test, var_predecidas)
rec_arbol = recall_score(y_test, var_predecidas, average='macro')
resultados['Árbol de decisión']        = acc_arbol
recall_resultados['Árbol de decisión'] = rec_arbol

m_arbol = confusion_matrix(y_test, var_predecidas, labels=CLASES)
print("\n" + "="*60)
print("ARBOL DE DECISION")
print("="*60)
_imprimir_confusion(m_arbol, 'Árbol de decisión')
print(classification_report(y_test, var_predecidas, target_names=CLASES))
_graficar_confusion(m_arbol, f'Árbol de decisión (exactitud {acc_arbol:.0%})', 'confusion_arbol.png')

# ── K Vecinos más cercanos (KNN) ──────────────────────────────────────────────
knn_clf = KNeighborsClassifier(n_neighbors=5)
knn_clf.fit(X_train, y_train)
y_pred_knn = knn_clf.predict(X_test)

acc_knn = accuracy_score(y_test, y_pred_knn)
rec_knn = recall_score(y_test, y_pred_knn, average='macro')
resultados['KNN (k=5)']        = acc_knn
recall_resultados['KNN (k=5)'] = rec_knn

m_knn = confusion_matrix(y_test, y_pred_knn, labels=CLASES)
print("\n" + "="*60)
print("KNN (k=5)")
print("="*60)
_imprimir_confusion(m_knn, 'KNN k=5')
print(classification_report(y_test, y_pred_knn, target_names=CLASES))
_graficar_confusion(m_knn, f'KNN k=5 (exactitud {acc_knn:.0%})', 'confusion_knn.png')

# ── Regresión Logística ───────────────────────────────────────────────────────
modelo_lr = LogisticRegression(random_state=0, max_iter=500)
modelo_lr.fit(X_train, y_train)
y_pred_lr = modelo_lr.predict(X_test)

acc_lr = accuracy_score(y_test, y_pred_lr)
rec_lr = recall_score(y_test, y_pred_lr, average='macro')
resultados['Regresión logística']        = acc_lr
recall_resultados['Regresión logística'] = rec_lr

m_lr = confusion_matrix(y_test, y_pred_lr, labels=CLASES)
print("\n" + "="*60)
print("REGRESION LOGISTICA")
print("="*60)
_imprimir_confusion(m_lr, 'Regresión logística')
print(classification_report(y_test, y_pred_lr, target_names=CLASES))
_graficar_confusion(m_lr, f'Regresión logística (exactitud {acc_lr:.0%})', 'confusion_logistica.png')

# ── Resumen comparativo en consola ────────────────────────────────────────────
print("\n" + "="*60)
print(f"{'Modelo':<25} {'Accuracy':>10} {'Recall':>10}")
print("-"*60)
for nombre in resultados:
    print(f"  {nombre:<23} {resultados[nombre]:>10.2%} {recall_resultados[nombre]:>10.2%}")
print("="*60)
mejor = max(resultados, key=resultados.get)
print(f"\nMejor modelo por accuracy: {mejor} ({resultados[mejor]:.2%})")

# ── Gráfico de barras comparativo ─────────────────────────────────────────────
nombres    = list(resultados.keys())
accuracies = list(resultados.values())
recalls    = list(recall_resultados.values())

x     = np.arange(len(nombres))
ancho = 0.35

fig, ax = plt.subplots(figsize=(9, 5))
barras_acc = ax.bar(x - ancho/2, accuracies, ancho, label='Accuracy',      color='steelblue')
barras_rec = ax.bar(x + ancho/2, recalls,    ancho, label='Recall (macro)', color='seagreen')

for barra in barras_acc:
    ax.text(barra.get_x() + barra.get_width() / 2,
            barra.get_height() + 0.01,
            f'{barra.get_height():.0%}',
            ha='center', va='bottom', fontsize=9)

for barra in barras_rec:
    ax.text(barra.get_x() + barra.get_width() / 2,
            barra.get_height() + 0.01,
            f'{barra.get_height():.0%}',
            ha='center', va='bottom', fontsize=9)

ax.set_ylabel('Puntuación')
ax.set_title('Comparación de modelos: Accuracy vs Recall')
ax.set_xticks(x)
ax.set_xticklabels(nombres, rotation=12, ha='right')
ax.set_ylim(0, 1.15)
ax.legend()
ax.grid(axis='y', linestyle='--', alpha=0.4)

plt.tight_layout()
ruta_barras = os.path.join(carpeta, 'comparacion_modelos.png')
plt.savefig(ruta_barras, dpi=150, bbox_inches='tight')
plt.close()
print(f"\nGráfico comparativo guardado: {ruta_barras}")
