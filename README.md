# Optimizacion del Mecanismo de Exoesqueleto de Mano

Optimizacion de los parametros dimensionales de un mecanismo de exoesqueleto de mano para rehabilitacion de dedos, buscando replicar el movimiento natural de agarre del dedo indice a partir de datos de captura de movimiento.

## Descripcion del Proyecto

Este proyecto implementa un proceso de optimizacion para un mecanismo de exoesqueleto de mano (mecanismos acoplados de 5 barras y 4 barras) que actua sobre las articulaciones IFP, IFD y la punta del dedo. El objetivo es encontrar las dimensiones optimas de los eslabones (18 parametros de diseno) que permitan al mecanismo reproducir la trayectoria natural del dedo durante el agarre, medida con un sistema de captura de movimiento.

## Archivos del Repositorio

| Archivo | Descripcion |
|---------|-------------|
| `exo_17.py` | Script principal de optimizacion en Python. Utiliza Optuna (optimizacion bayesiana), Evolucion Diferencial (scipy), distancia de Chamfer, alineacion rigida basada en SVD, filtrado Savitzky-Golay y penalizacion anti-gancho. Optimiza 18 parametros de diseno. |
| `CinematicaExoFinal.m` | Script de MATLAB con el analisis de posicion cinematica original del mecanismo (mecanismos acoplados de 5 barras y 4 barras para las articulaciones IFP, IFD y punta). |
| `mocap_indice_120pts.csv` | Datos de captura de movimiento del dedo indice (120 puntos con angulos MCP, PIP, DIP). |
| `ReporteModeloCinematicoExoPavel.pdf` | Reporte tecnico del modelo cinematico del exoesqueleto. |
| `result_exo_17.png` | Grafica de resultados de la optimizacion. |

## Como Ejecutar la Optimizacion

### Requisitos

- Python 3.8+
- Dependencias:

```bash
pip install numpy scipy optuna matplotlib pandas
```

### Ejecucion

```bash
python exo_17.py
```

El script ejecutara la optimizacion utilizando Optuna para la busqueda global de hiperparametros y Evolucion Diferencial como optimizador interno, minimizando la distancia de Chamfer entre la trayectoria generada por el mecanismo y los datos de captura de movimiento.

## Enfoque de Optimizacion

El proceso de optimizacion combina las siguientes tecnicas:

1. **Optuna (Optimizacion Bayesiana):** Busqueda inteligente del espacio de 18 parametros de diseno mediante muestreo bayesiano.
2. **Evolucion Diferencial (scipy):** Optimizador local utilizado dentro de cada trial de Optuna para refinar los parametros.
3. **Distancia de Chamfer:** Metrica de comparacion entre la trayectoria generada por el mecanismo y los datos de referencia de captura de movimiento.
4. **Alineacion rigida (SVD):** Alineacion optima de las trayectorias antes de calcular la distancia.
5. **Filtrado Savitzky-Golay:** Suavizado de las trayectorias para reducir ruido.
6. **Penalizacion anti-gancho:** Restriccion que evita configuraciones no deseadas del mecanismo.

## Base de Datos de Captura de Movimiento

Los datos de captura de movimiento utilizados en este proyecto provienen de la base de datos Pizzoleto:

- **Repositorio:** [Base-de-Datos-Pizzoleto](https://github.com/Utemita/Base-de-Datos-Pizzoleto)
- **Dataset en Zenodo:** [DOI: 10.5281/zenodo.20561611](https://doi.org/10.5281/zenodo.20561611)

## Autora

**Adriana Elorza Ramos**
Universidad Tecnologica de la Mixteca
