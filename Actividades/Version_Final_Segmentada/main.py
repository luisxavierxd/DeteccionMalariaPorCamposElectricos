import os

from campoElectrico import calcular_campo
from trayectorias   import simular, exportar_csv
from grafico        import graficar

carpeta = os.path.dirname(os.path.abspath(__file__))


def main():
    print("Calculando campo eléctrico de fondo...")
    xx, yy, X, Y, Ex, Ey, V = calcular_campo()

    print("Simulando trayectorias celulares...")
    trayectorias = simular()

    ruta_csv = os.path.join(carpeta, 'trayectorias_malaria.csv')
    exportar_csv(trayectorias, ruta_csv)

    ruta_png = os.path.join(carpeta, 'trayectorias_malaria.png')
    graficar(xx, yy, X, Y, Ex, Ey, V, trayectorias, ruta_png)


if __name__ == '__main__':
    main()
