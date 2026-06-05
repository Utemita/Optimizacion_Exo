# Analisis de Modificacion de la Articulacion DIP del Exoesqueleto de Mano

## 1. Planteamiento del Problema

La articulacion interfalangica distal (DIP) del exoesqueleto de mano no posee movimiento independiente. En el modelo cinematico actual, implementado en `CinematicaExoFinal.m`, el angulo de la falange distal se calcula como un offset constante respecto a la falange medial:

```matlab
THETAfd(j) = THETAfm(j) + THETAauxfd;
```

Donde `THETAauxfd = 38.78` grados es una constante fija.

**Consecuencia directa:** La falange distal esta rigidamente unida a la falange medial con un angulo relativo constante. Esto significa que:

- El angulo relativo DIP-PIP es siempre 38.78 grados sin importar la posicion del mecanismo.
- La articulacion DIP no presenta flexion/extension durante el ciclo de agarre.
- En el modelo CAD, la falange distal se mueve solidariamente con la medial, sin articulacion propia.
- Los datos de captura de movimiento (MOCAP) muestran un rango de movimiento DIP de aproximadamente 20 grados (de -11 a 9 grados), lo cual no es reproducido por el modelo actual.

## 2. Analisis de Causa Raiz

### Cadena Cinematica Completa: Motor a DIP

La cadena cinematica del exoesqueleto transmite el movimiento desde un unico motor hasta las tres falanges del dedo a traves de la siguiente secuencia:

| Etapa | Mecanismo | Entrada | Salida | DOF Independiente |
|-------|-----------|---------|--------|-------------------|
| 1 | Par de engranes | theta_motor | theta_1, theta_2 | Si (relacion fija reng=2) |
| 2 | Mecanismo de 5 barras #1 | theta_1, theta_2 | Punto P (pxP, pyP) | Si |
| 3 | Mecanismo de 4 barras #1 | theta_1 | theta_fp (falange proximal) | Si |
| 4 | Mecanismo de 5 barras #2 | Punto P, S1, S2, M4 | Punto P2 (pxP2, pyP2) | Si |
| 5 | Mecanismo de 4 barras #2 | P2, S2, IFP | theta_fm (falange medial) | Si |
| 6 | **Offset constante** | theta_fm | **theta_fd = theta_fm + 38.78** | **NO** |

### Diagrama de la Cadena de Transmision

```
MOTOR (1 DOF)
    |
    v
[Par de Engranes] ---> reng = 2
    |          |
    v          v
theta_1      theta_2
    |          |
    v          v
[Mec. 4B #1] [Mec. 5B #1]
    |              |
    v              v
theta_fp        Punto P
    |              |
    |   [Mec. 5B #2]
    |        |
    v        v
   S2      Punto P2
    |        |
    v        v
    [Mec. 4B #2]
         |
         v
    theta_fm (falange medial)
         |
         v
    theta_fd = theta_fm + CONSTANTE  <--- PROBLEMA: Sin mecanismo
```

### Identificacion del Problema

Los mecanismos 1-5 proporcionan movimiento independiente y variable a las falanges proximal y medial. Sin embargo, la falange distal simplemente hereda el angulo de la medial mas un offset fijo. No existe un mecanismo fisico que transmita movimiento diferencial a la articulacion DIP.

### Ecuacion Problematica (linea en CinematicaExoFinal.m):

```matlab
% Calculando la posicion de la falange distal
THETAfd(j) = THETAfm(j) + THETAauxfd;  % THETAauxfd = 38.78 (constante)
```

El angulo relativo entre la falange distal y la medial es:

```
DIP_relativo = THETAfd - THETAfm = THETAauxfd = 38.78 = constante
```

Esto viola la cinematica natural del dedo, donde la articulacion DIP se flexiona de manera acoplada (pero no rigida) con la articulacion PIP.

## 3. Solucion Propuesta - Tercer Mecanismo de 4 Barras

### Concepto de Diseno

Se propone agregar un tercer mecanismo de 4 barras que proporcione movimiento independiente a la falange distal. Este mecanismo se monta sobre la falange medial y su entrada se deriva del movimiento RELATIVO entre la falange proximal y la falange medial.

### Principio de Funcionamiento

Un "eslabon remoto" (Link11) conecta un punto en la falange proximal (cercano al soporte S2) con un punto pivote (soporte S3) en la falange medial. Conforme la articulacion PIP se flexiona (theta_fm cambia respecto a theta_fp), este eslabon remoto cambia su angulo relativo a la falange medial, proporcionando la entrada tipo manivela para el tercer mecanismo de 4 barras.

### Componentes del Tercer Mecanismo

- **Bancada (d3):** Distancia a lo largo de la falange medial entre el soporte S3 y la articulacion IFD (punto de salida)
- **Manivela/Entrada (a3 = Link11):** Eslabon remoto desde el soporte S2 (en falange proximal) hasta el pivote de entrada en S3 (en falange medial)
- **Acoplador (b3 = Link12):** Eslabon de conexion intermedio
- **Balancin/Salida (c3 = Link13):** Determina el angulo de la falange distal respecto a la medial

### Diagrama de Topologia del Mecanismo

```
    FALANGE PROXIMAL                    FALANGE MEDIAL                    FALANGE DISTAL
    ==================                  ====================              ================
         |                                   |                                 |
         |  S2                               |  S3              IFD            |
         |  (soporte en                      |  (soporte en     (articulacion  |
         |   proximal)                       |   medial)         DIP)          |
         |   *                               |   *               *             |
         |   |                               |   |               |             |
         |   |                               |   |    d3         |             |
         |   |         Link11 (a3)           |   |<------------>|             |
         |   |==============================>|   |               |             |
         |   |   (eslabon remoto)            |   |               |             |
         |                                   |   |--Link12(b3)-->|             |
         |                                   |   |               |             |
         |                                   |   |  Link13(c3)   |             |
         |                                   |   |<--------------|             |
         |                                   |                   |             |
    ==================                  ====================     |  ============
                                                                 |
                                                            theta_fd (variable!)
```

### Vista Esquematica del Mecanismo de 4 Barras #3

```
                S2 (en proximal)         S3 (en medial)          IFD
                |                        |                       |
  PROXIMAL      |------- Link11 -------->|                       |
  FALANGE       |       (a3=25mm)        |    FALANGE MEDIAL     |
                |                        |---d3---[4-bar #3]--->| DISTAL
                |                        |  a3 = Link11          |
                                         |  b3 = Link12          |
                                         |  c3 = Link13          |
                                         |  d3 = dist IFP-S3     |
                                         |       (ground)        |
```

### Ventaja del Diseno

Al utilizar el movimiento relativo entre las falanges proximal y medial como entrada, el tercer mecanismo logra:
1. Movimiento independiente de la falange distal sin actuadores adicionales.
2. Acoplamiento cinematico natural (DIP se flexiona cuando PIP se flexiona).
3. Relacion de transmision variable (no lineal) controlada por la geometria del mecanismo.
4. Compatible con la topologia existente (no requiere modificar mecanismos 1-5).

## 4. Ecuaciones Cinematicas del Tercer Mecanismo

### 4.1 Sistema de Referencia Flotante en la Falange Medial

El tercer mecanismo de 4 barras se analiza en un sistema de referencia local fijado a la falange medial:

- **Origen:** Articulacion IFP (punto de union entre falanges proximal y medial)
- **Eje X local:** Direccion de la falange medial (angulo `theta_fm` en el sistema global)
- **Eje Y local:** Perpendicular a la falange medial (sentido antihorario)

La transformacion entre el sistema global y el sistema local de la falange medial es:

```
x_local = (x_global - pxIFP)*cos(theta_fm) + (y_global - pyIFP)*sin(theta_fm)
y_local = -(x_global - pxIFP)*sin(theta_fm) + (y_global - pyIFP)*cos(theta_fm)
```

### 4.2 Angulo de Entrada (Crank) del Tercer Mecanismo

La entrada del mecanismo se calcula a partir de la posicion del soporte S2 (en falange proximal) relativa al soporte S3 (en falange medial).

**Posicion del soporte S3 sobre la falange medial:**

S3 se ubica a una distancia `dsm3` del IFP a lo largo de la falange medial, con una altura `hsm3` perpendicular:

```
THETAps3 = THETAfm - atand(hsm3/dsm3)
rs3 = sqrt(hsm3^2 + dsm3^2)
pxS3 = pxIFP + rs3*cos(deg2rad(THETAps3))
pyS3 = pyIFP + rs3*sin(deg2rad(THETAps3))
```

**Angulo de la manivela en el sistema global:**

```
THETA_crank3 = atan2d(pyS2 - pyS3, pxS2 - pxS3)
```

**Angulo de la manivela en el sistema local de la falange medial:**

```
theta_crank3_local = deg2rad(THETA_crank3) - thetafm
```

Este angulo cambia cuando la articulacion PIP se flexiona, ya que S2 esta fijo a la falange proximal mientras S3 esta fijo a la medial.

### 4.3 Resolucion del Mecanismo de 4 Barras (Metodo de Media Tangente)

Se utiliza el mismo metodo de resolucion que los mecanismos de 4 barras existentes (#1 y #2), basado en la sustitucion de media tangente `t = tan(theta/2)`.

**Interpretacion geometrica del eslabon de bancada (ground link):** En la formulacion local del tercer mecanismo, `d3 = sqrt(hsm3^2 + dsm3^2)` representa la distancia en linea recta desde el origen del sistema de referencia local (ubicado en IFP) hasta el pivote de la manivela en S3. Es decir, `d3` es el eslabon de bancada (ground link) que conecta el origen de referencia con el pivote del crank. El pivote de salida del balancin (rocker) se conecta a la articulacion DIP, la cual se encuentra a una distancia `fm` del IFP a lo largo de la falange medial. La relacion angular entre el balancin y la orientacion real de la falange distal se absorbe en el offset geometrico `THETAaux_fd3`.

**Parametros del mecanismo:**
- `a3 = Link11` (longitud de la manivela/eslabon remoto)
- `b3 = Link12` (longitud del acoplador)
- `c3 = Link13` (longitud del balancin/salida)
- `d3 = sqrt(hsm3^2 + dsm3^2)` (longitud de bancada)
- `theta_base3 = 0` (bancada alineada con eje local X de falange medial)
- `theta_input3 = theta_crank3_local` (angulo de entrada calculado en 4.2)

**Ecuaciones de resolucion:**

```
k1_3 = d3*cos(theta_base3) + a3*cos(theta_crank3_local)
k2_3 = d3*sin(theta_base3) + a3*sin(theta_crank3_local)
k3_3 = k1_3^2 + k2_3^2 + c3^2 - b3^2
```

**Coeficientes de la ecuacion cuadratica:**

```
A1_3 = -k3_3 - 2*k1_3*c3
B1_3 = 4*k2_3*c3
C1_3 = 2*k1_3*c3 - k3_3
```

**Discriminante y solucion:**

```
disc3 = B1_3^2 - 4*A1_3*C1_3
```

Si `disc3 >= 0`:
```
tan_theta4_3 = (-B1_3 - sqrt(disc3)) / (2*A1_3)    % Configuracion abierta
theta4_3 = 2*atan(tan_theta4_3)
```

Si `disc3 < 0`: El mecanismo no puede ensamblarse en esa posicion. Se utiliza el valor de fallback (offset constante original).

### 4.4 Angulo de Salida de la Falange Distal

El angulo absoluto de la falange distal en el sistema global es:

```
theta_fd = theta_fm + theta4_3 + THETAaux_fd3
```

Donde:
- `theta_fm`: angulo de la falange medial (en grados o radianes segun contexto)
- `theta4_3`: angulo de salida del tercer mecanismo de 4 barras (en el sistema local de la falange medial)
- `THETAaux_fd3`: offset geometrico entre el angulo del balancin y la orientacion real de la falange distal (parametro de diseno, valor inicial: 30 grados)

**Nota fundamental:** `theta_fd` ya NO es constante respecto a `theta_fm`. Conforme la articulacion PIP se flexiona (theta_fm cambia respecto a theta_fp), el angulo de la manivela `theta_crank3_local` cambia, lo cual acciona el mecanismo de 4 barras y produce un `theta4_3` variable.

### 4.5 Consideraciones sobre la Aproximacion

La formulacion presentada emplea varias simplificaciones que conviene explicitar:

**Longitud de manivela fija (`a3 = Link11`):** En la implementacion, `a3 = Link11` se utiliza como una longitud de manivela constante en las ecuaciones del mecanismo de 4 barras. Sin embargo, la distancia real `|S2 - S3|` varia con la flexion del PIP, ya que S2 esta fijo a la falange proximal y S3 esta fijo a la falange medial. Esta es una aproximacion valida cuando: (a) la variacion en `|S2 - S3|` es pequena en relacion con Link11, o (b) el optimizador ajusta los parametros para minimizar esta discrepancia. Para una formulacion mas rigurosa, se podria calcular `a3_actual(j) = |S2(j) - S3(j)|` en cada paso, pero esto convierte al mecanismo de 4 barras en un mecanismo de geometria variable que requiere solucion iterativa. La aproximacion de longitud fija es aceptable como punto de partida para el diseno, dado que el optimizador encontrara conjuntos de parametros donde el mecanismo se ensambla correctamente en todo el rango de movimiento.

**Eslabon de bancada `d3` y relaciones geometricas:** El eslabon de bancada `d3` representa la distancia IFP-S3 en el marco de referencia local. La salida del balancin se conecta a la articulacion DIP (que esta a una distancia `fm` del IFP a lo largo de la falange medial). Estas relaciones geometricas se manejan mediante el angulo de offset `THETAaux_fd3`, que absorbe la diferencia angular entre la orientacion del balancin y la orientacion real de la falange distal.

## 5. Tabla de Parametros Nuevos

| Parametro | Simbolo | Significado Fisico | Valor Inicial | Unidad |
|-----------|---------|-------------------|---------------|--------|
| Link11 | a3 | Longitud del eslabon remoto (S2 en proximal hasta pivote en medial S3) | 25 | mm |
| Link12 | b3 | Longitud del acoplador del tercer mecanismo de 4 barras | 20 | mm |
| Link13 | c3 | Longitud del balancin (salida) del tercer mecanismo de 4 barras | 15 | mm |
| dsm3 | - | Distancia horizontal desde IFP hasta el soporte S3 a lo largo de la falange medial | 12 | mm |
| hsm3 | - | Altura del soporte S3 perpendicular a la falange medial | 10 | mm |
| THETAaux_fd3 | - | Offset angular entre el balancin del tercer mecanismo y la falange distal | 30 | grados |

**Nota:** La bancada del tercer mecanismo se calcula como `d3 = sqrt(hsm3^2 + dsm3^2)`, por lo que no es un parametro independiente.

### Relacion con Parametros Existentes

Los parametros existentes que mantienen su funcionalidad sin modificacion son:

| Parametro | Valor | Funcion |
|-----------|-------|---------|
| THETAauxfd | 38.78 | Se conserva como valor de FALLBACK cuando el tercer mecanismo no converge |
| fp | 49 | Largo de falange proximal (define posicion de S2) |
| fm | 26 | Largo de falange medial (define geometria del tercer mecanismo) |
| fd | 24 | Largo de falange distal (no cambia) |
| hsp | 17 | Altura del soporte proximal (define S2) |
| dsp | 18 | Distancia del soporte proximal (define S2) |

## 6. Criterios de Verificacion

Para validar que el tercer mecanismo funciona correctamente, se deben cumplir los siguientes criterios:

### 6.1 Rango de Movimiento

- La articulacion DIP debe tener un rango de movimiento independiente de **minimo 20 grados** y **maximo 40 grados** durante el ciclo completo de cierre.
- Los datos MOCAP muestran un rango DIP de aproximadamente 20 grados (de -11 a 9 grados despues de invertir la convencion de signos).

### 6.2 Monotonicidad

- El movimiento de la articulacion DIP debe ser **monotono** (sin inversiones de direccion).
- La falange distal debe flexionarse de manera progresiva durante el cierre del dedo.
- No se admiten oscilaciones o movimientos erraticos.

### 6.3 Relacion de Acoplamiento

- La relacion entre el rango DIP y el rango PIP debe estar entre **0.3 y 0.7**.
- Fisiologicamente, la articulacion DIP se mueve aproximadamente 2/3 del rango del PIP.
- Rango PIP tipico: 0-74 grados (de datos MOCAP).
- Rango DIP esperado: 20-50 grados.

### 6.4 Convergencia del Mecanismo

- El discriminante `disc3 = B1_3^2 - 4*A1_3*C1_3` debe ser **no negativo** en todo el rango de operacion.
- Si el discriminante es negativo en alguna posicion, los parametros deben ajustarse.
- El fallback (offset constante) solo debe activarse en condiciones excepcionales, no durante operacion normal.

### 6.5 Continuidad

- El angulo `theta4_3` debe ser una funcion continua de la posicion de entrada.
- No se admiten saltos discontinuos en el angulo de la falange distal.

## 7. Codigo MATLAB Modificado

La implementacion completa se encuentra en el archivo `CinematicaExoModificada.m`. A continuacion se muestra el fragmento clave que reemplaza el calculo de offset constante:

### Parametros nuevos (al inicio del script):

```matlab
% Tercer mecanismo de 4 barras (articulacion DIP)
Link11 = 25;         % Eslabon remoto (S2 -> pivot en medial)
Link12 = 20;         % Acoplador del tercer mecanismo
Link13 = 15;         % Balancin del tercer mecanismo
dsm3 = 12;           % Distancia soporte S3 a IFP sobre falange medial
hsm3 = 10;           % Altura del soporte S3 perpendicular a falange medial
THETAaux_fd3 = 30;   % Offset angular para falange distal (grados)
```

### Asignaciones de variables:

```matlab
% Tercer mecanismo de 4 barras
a3 = Link11;
b3 = Link12;
c3 = Link13;
d3 = sqrt(hsm3^2 + dsm3^2);  % Hipotenusa del soporte S3 en falange medial
```

### Reemplazo dentro del bucle (despues del calculo de pyIFD):

```matlab
% TERCER MECANISMO DE 4 BARRAS (Articulacion DIP)
% Calculo del soporte S3 sobre la falange medial
THETAps3(j) = THETAfm(j) - atand(hsm3/dsm3);
thetaps3(j) = deg2rad(THETAps3(j));
rs3 = sqrt(hsm3^2 + dsm3^2);
pxS3(j) = pxIFP(j) + rs3*cos(thetaps3(j));
pyS3(j) = pyIFP(j) + rs3*sin(thetaps3(j));

% Angulo de entrada del tercer mecanismo (enlace remoto desde S2)
THETA_crank3(j) = atan2d(pyS2(j) - pyS3(j), pxS2(j) - pxS3(j));
theta_crank3_local(j) = deg2rad(THETA_crank3(j)) - thetafm(j);

% Bancada del tercer mecanismo en sistema local de la falange medial
theta_base3 = 0;

% Resolucion del mecanismo de 4 barras #3
k1_3 = d3*cos(theta_base3) + a3*cos(theta_crank3_local(j));
k2_3 = d3*sin(theta_base3) + a3*sin(theta_crank3_local(j));
k3_3 = k1_3^2 + k2_3^2 + c3^2 - b3^2;
A1_3 = -k3_3 - 2*k1_3*c3;
B1_3 = 4*k2_3*c3;
C1_3 = 2*k1_3*c3 - k3_3;

disc3 = B1_3^2 - 4*A1_3*C1_3;
if disc3 < 0
    warning('Discriminante negativo en mecanismo 3, paso %d', j);
    THETAfd(j) = THETAfm(j) + THETAauxfd;  % Fallback al offset constante
else
    tantheta4_3 = (-B1_3 - sqrt(disc3)) / (2*A1_3);  % Config abierta
    theta4_3(j) = 2*atan(tantheta4_3);
    THETA4_3(j) = rad2deg(theta4_3(j));
    % Angulo de la falange distal en sistema global
    THETAfd(j) = THETAfm(j) + THETA4_3(j) + THETAaux_fd3;
end
```

### Verificacion despues del bucle:

```matlab
% VERIFICACION: Rango de movimiento de la articulacion DIP
DIP_relative = THETAfd - THETAfm;
fprintf('Rango de movimiento DIP (relativo a medial): %.2f a %.2f grados\n', ...
        min(DIP_relative), max(DIP_relative));
fprintf('Excursion DIP: %.2f grados\n', max(DIP_relative) - min(DIP_relative));
if max(DIP_relative) - min(DIP_relative) < 5
    warning('La articulacion DIP tiene menos de 5 grados de excursion. Revisar parametros.');
end
```

## 8. Impacto en Optimizacion Python (exo_18.py)

### 8.1 Nuevos Parametros en el Vector de Optimizacion

Se agregan 5 nuevos parametros al vector de diseno:

| Indice | Parametro | Valor Inicial | Limites Inferiores | Limites Superiores |
|--------|-----------|---------------|--------------------|--------------------|
| 17 | Link11 | 0.025 m | 0.015 | 0.040 |
| 18 | Link12 | 0.020 m | 0.010 | 0.035 |
| 19 | Link13 | 0.015 m | 0.008 | 0.030 |
| 20 | dsm3 | 0.012 m | 0.005 | 0.020 |
| 21 | THETAaux_fd3 | 0.5236 rad (30 deg) | 0.0 | 1.2217 (70 deg) |

**Nota:** `hsm3` se puede fijar como constante o derivarse de otros parametros para limitar el espacio de busqueda. Se recomienda fijarlo inicialmente en 10 mm.

### 8.2 Total de Parametros

- Parametros originales: 17
- Parametros nuevos: 5
- **Total: 22 parametros**

### 8.3 Modificacion en `run_kinematics`

En la funcion `run_kinematics` de `exo_18.py`, se debe reemplazar:

```python
# ANTES (offset constante):
theta_fd = theta_fm + theta_aux_fd
```

Con el calculo del tercer mecanismo de 4 barras:

```python
# DESPUES (tercer mecanismo de 4 barras):
# Posicion de S3 en falange medial
rs3 = np.sqrt(hsm3**2 + dsm3**2)
theta_ps3 = theta_fm_deg - np.degrees(np.arctan2(hsm3, dsm3))
px_s3 = px_ifp + rs3 * np.cos(np.radians(theta_ps3))
py_s3 = py_ifp + rs3 * np.sin(np.radians(theta_ps3))

# Angulo de entrada (crank) en sistema local de la medial
theta_crank3_global = np.arctan2(py_s2 - py_s3, px_s2 - px_s3)
theta_crank3_local = theta_crank3_global - theta_fm_rad

# Resolucion 4-bar #3
k1_3 = d3 + a3 * np.cos(theta_crank3_local)
k2_3 = a3 * np.sin(theta_crank3_local)
k3_3 = k1_3**2 + k2_3**2 + c3**2 - b3**2
A1_3 = -k3_3 - 2*k1_3*c3
B1_3 = 4*k2_3*c3
C1_3 = 2*k1_3*c3 - k3_3

disc3 = B1_3**2 - 4*A1_3*C1_3
if disc3 >= 0:
    tan_t4_3 = (-B1_3 - np.sqrt(disc3)) / (2*A1_3)
    theta4_3 = 2*np.arctan(tan_t4_3)
    theta_fd = theta_fm_deg + np.degrees(theta4_3) + theta_aux_fd3_deg
else:
    theta_fd = theta_fm + theta_aux_fd  # Fallback
```

### 8.4 Funcion de Aptitud (Fitness)

La funcion de aptitud no requiere modificaciones estructurales. Ya evalua la trayectoria de la punta del dedo (pxPF, pyPF), que se calcula a partir de `theta_fd`. Al modificar el calculo de `theta_fd`, la funcion de aptitud automaticamente optimizara los nuevos parametros para que la trayectoria se ajuste a los datos MOCAP.

### 8.5 Consideraciones para la Optimizacion

1. **Restriccion de Grashof:** Los parametros Link11, Link12, Link13, dsm3 deben satisfacer la condicion de ensamble del mecanismo de 4 barras en todo el rango de operacion.
2. **Penalizacion por discriminante negativo:** Si `disc3 < 0` en alguna posicion durante la simulacion, la funcion de aptitud debe penalizar fuertemente ese conjunto de parametros.
3. **Inicializacion:** Se recomienda inicializar los nuevos parametros con los valores de la tabla de la Seccion 5 para garantizar convergencia inicial.
4. **Espacio de busqueda ampliado:** Con 22 parametros, se recomienda aumentar el numero de iteraciones de Optuna para asegurar buena exploracion del espacio.
