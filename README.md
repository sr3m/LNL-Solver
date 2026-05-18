---
title: "Guía gráfica del LNL-Solver"
author: "Willmers Isaías Hernández Cornielle"
course: "Análisis Numérico"
lang: es
---

# Guía gráfica del LNL-Solver

**Linear & Non-Linear Solver**

**Cómo importar, usar y entender `LNL_Solver.py`**

**LNL-Solver** significa **Linear & Non-Linear Solver**. En Python se usa mediante el archivo importable `LNL_Solver.py`, que funciona como una pequeña librería personal para resolver sistemas lineales y sistemas no lineales desde notebooks de Python.

```mermaid
flowchart LR
    A["LNL_Solver.py<br/>archivo módulo"] --> B["LinearSystem<br/>sistemas Ax=b"]
    B --> C["NonLinearSystem<br/>sistemas F(x,y)=0"]
    B --> D["show_solution_set<br/>impresión ordenada"]
    C --> D
```

## Índice

- [Mapa general](#mapa-general)
- [Importar el solver en un notebook](#importar-el-solver-en-un-notebook)
- [Flujo para sistemas lineales](#flujo-para-sistemas-lineales)
- [Métodos lineales disponibles](#métodos-lineales-disponibles)
- [Criterios de parada lineales](#criterios-de-parada-lineales)
- [Qué devuelve un método lineal](#qué-devuelve-un-método-lineal)
- [Historial de iteraciones](#historial-de-iteraciones)
- [Flujo para sistemas no lineales](#flujo-para-sistemas-no-lineales)
- [Newton-Raphson](#newton-raphson)
- [Punto fijo](#punto-fijo-en-sistemas-no-lineales)
- [Cómo leer `show_solution_set`](#cómo-leer-show_solution_set)
- [Diagnóstico rápido de errores](#diagnóstico-rápido-de-errores)
- [Recetas de uso](#recetas-de-uso)
- [Ejercicio guiado: Solver + `numpy` + `matplotlib`](#ejercicio-guiado-solver--numpy--matplotlib)
- [Checklist antes de entregar](#checklist-antes-de-entregar)

---

## Mapa general

```mermaid
flowchart TB
    M["Módulo<br/><code>LNL_Solver.py</code>"]
    L["Clase<br/><code>LinearSystem</code><br/>Ax=b"]
    N["Clase<br/><code>NonLinearSystem</code><br/>F(x,y)=0"]
    S["Función<br/><code>show_solution_set</code><br/>salida legible"]
    J["Jacobi"]
    G["Gauss-Seidel"]
    D["Matriz inversa"]
    NR["Newton-Raphson"]
    PF["Punto fijo"]

    M --> L
    M --> N
    M --> S
    L --> J
    L --> G
    L --> D
    N --> NR
    N --> PF
```

| Objeto | Para qué sirve |
|---|---|
| `LinearSystem` | Resolver sistemas lineales \(Ax=b\), revisar dominancia diagonal y usar Jacobi, Gauss-Seidel o método directo. |
| `NonLinearSystem` | Resolver sistemas no lineales de dos variables mediante Newton-Raphson o punto fijo. |
| `show_solution_set` | Imprimir una solución con nombres de variables, número de iteraciones y error. |

> **Idea central.** En vez de reescribir las clases en cada notebook, se guarda el solver en un archivo `.py` y se importa cuando se necesita resolver un sistema.

---

## Importar el Solver en un Notebook

### Caso recomendado: el notebook está en la misma carpeta

```mermaid
flowchart TB
    subgraph F["Carpeta de trabajo"]
        A["mi_notebook.ipynb"]
        B["LNL_Solver.py"]
    end
    A --> C["El import funciona directamente"]
    B --> C
```

La celda de importación debe verse así:

```python
import numpy as np
from LNL_Solver import *
```

> **Detalle importante.** En Python se importa el nombre del módulo, no el nombre del archivo con extensión. Además, aunque el nombre visual es **LNL-Solver**, el módulo usa guion bajo porque Python no permite guiones en nombres importados con `from`.

Correcto:

```python
from LNL_Solver import *
```

Incorrecto:

```python
from LNL_Solver.py import *
```

### Si el notebook está en otra carpeta

```mermaid
flowchart LR
    A["Notebook en otra carpeta"] --> B["Agregar ruta con sys.path"]
    B --> C["Importar LNL_Solver"]
```

```python
import sys
import numpy as np

sys.path.append(
    r"/ruta/a/la/carpeta/donde/esta/LNL_Solver"
)

from LNL_Solver import *
```

Si editas `LNL_Solver.py` y el notebook no reconoce los cambios, reinicia el kernel o recarga el módulo:

```python
import LNL_Solver
from importlib import reload

reload(LNL_Solver)
from LNL_Solver import *
```

---

## Flujo Para Sistemas Lineales

```mermaid
flowchart TB
    A["1. Definir A"] --> B["2. Definir b"]
    B --> C["3. Elegir x0"]
    A --> D["4. Crear LinearSystem"]
    B --> D
    C --> D
    D --> E["5. Revisar convergencia"]
    D --> F["6. Resolver"]
    F --> G["7. Mostrar resultado"]
```

### Plantilla mínima

```python
import numpy as np
from LNL_Solver import *

A = np.array([
    [-2, 1, 0, 0],
    [1, -2, 1, 0],
    [0, 1, -2, 1],
    [0, 0, 1, -2]
])

b = np.array([[1], [2], [-7], [-1]])
x0 = [0, 0, 0, 0]

system = LinearSystem(A, b, x0)
system.convergence()

solution, iterations, error = system.solved_by_gauss_seidel(
    'max_iteration',
    max_iteration=5
)

show_solution_set(solution, iterations, 'x', error)
```

| Línea | Significado |
|---|---|
| `A = np.array(...)` | Matriz de coeficientes del sistema. |
| `b = np.array(...)` | Vector de términos independientes. |
| `x0 = [...]` | Punto inicial para métodos iterativos. |
| `LinearSystem(A,b,x0)` | Construye el objeto que sabe resolver el sistema. |
| `convergence()` | Indica si la matriz es diagonalmente dominante. |
| `solved_by_gauss_seidel(...)` | Ejecuta Gauss-Seidel hasta cumplir el criterio elegido. |

---

## Métodos Lineales Disponibles

```mermaid
flowchart TB
    A["LinearSystem"]
    A --> B["solved_by_jacobi<br/>usa valores viejos"]
    A --> C["solved_by_gauss_seidel<br/>usa valores nuevos"]
    A --> D["solved_by_matrix_direct_method<br/>usa A^{-1}b"]
```

| Método | Cuándo usarlo | Llamada típica |
|---|---|---|
| Jacobi | Comparar métodos iterativos o cuando se quiere una actualización por bloque. | `solved_by_jacobi(...)` |
| Gauss-Seidel | Método iterativo usual; aprovecha los valores recién calculados. | `solved_by_gauss_seidel(...)` |
| Directo | Para obtener una solución exacta numérica si \(A\) es invertible. | `solved_by_matrix_direct_method()` |

| Método | Idea matemática | Detalle computacional |
|---|---|---|
| Jacobi | Calcula cada \(x_i^{(k+1)}\) usando solo valores de \(x^{(k)}\). | Es más fácil de comparar y paralelizar, pero suele necesitar más iteraciones. |
| Gauss-Seidel | Calcula \(x_i^{(k+1)}\) usando los valores nuevos apenas están disponibles. | Suele converger más rápido que Jacobi cuando la matriz se presta para iteración. |
| Directo | Calcula \(x=A^{-1}b\). | No genera historial iterativo; sirve como referencia si \(A\) es cuadrada e invertible. |

Jacobi:

$$
x_i^{(k+1)}=
\frac{1}{a_{ii}}
\left(
b_i-\sum_{j\ne i}a_{ij}x_j^{(k)}
\right)
$$

Gauss-Seidel:

$$
x_i^{(k+1)}=
\frac{1}{a_{ii}}
\left(
b_i-\sum_{j<i}a_{ij}x_j^{(k+1)}
-\sum_{j>i}a_{ij}x_j^{(k)}
\right)
$$

> **Una llamada sin argumentos hace una sola iteración.**

```python
system.solved_by_gauss_seidel()
```

Si quieres que el método se detenga automáticamente, debes pasar un criterio de parada:

```python
solution, iterations, error = system.solved_by_gauss_seidel(
    'residue',
    tol=1e-5
)
```

---

## Criterios de Parada Lineales

```mermaid
flowchart LR
    A["Iteración nueva<br/>x_{k+1}"] --> B["Calcular error<br/>y residuo"]
    B --> C{"¿Cumple criterio?"}
    C -->|Sí| D["Devolver solución"]
    C -->|No| E["Seguir iterando"]
    E --> A
```

| Criterio | Qué revisa | Ejemplo |
|---|---|---|
| `max_iteration` | Se detiene al llegar al número máximo de iteraciones. | `max_iteration=20` |
| `absolute_error` | Norma euclidiana de \(x_{k+1}-x_k\). | `tol=1e-4` |
| `max_absolute_error` | Mayor cambio componente a componente. | `tol=1e-4` |
| `relative_error` | Error absoluto dividido entre la norma del punto actual. | `tol=1e-5` |
| `residue` | Norma del residuo \(Ax-b\). | `tol=1e-5` |
| `combined` | Exige simultáneamente error absoluto máximo y residuo. | `absolute_tol=1e-4`, `residue_tol=1e-5` |

```python
# Ejemplos de criterios
system.solved_by_gauss_seidel('max_iteration', max_iteration=10)
system.solved_by_gauss_seidel('absolute_error', tol=1e-5)
system.solved_by_gauss_seidel('relative_error', tol=1e-5)
system.solved_by_gauss_seidel('residue', tol=1e-5)
system.solved_by_gauss_seidel('max_absolute_error', tol=1e-5)
system.solved_by_gauss_seidel(
    'combined',
    absolute_tol=1e-5,
    residue_tol=1e-5
)
```

### Cómo escoger el criterio lineal

| Necesidad | Criterio recomendado | Por qué |
|---|---|---|
| Solo repetir un número fijo de veces | `max_iteration` | Sirve para tareas donde piden exactamente \(n\) iteraciones. |
| Comparar convergencia entre métodos | `absolute_error` | Mide el tamaño del salto \(x_{k+1}-x_k\). |
| Controlar el peor cambio individual | `max_absolute_error` | Útil si ninguna variable debe moverse más que cierta tolerancia. |
| Variables con escalas muy distintas | `relative_error` | Normaliza el salto por el tamaño de la solución actual. |
| Comprobar que \(Ax=b\) se cumple | `residue` | Mide directamente el defecto \(Ax-b\). |
| Entrega más rigurosa | `combined` | Exige simultáneamente poco cambio y bajo residuo. |

> **Recomendación de uso.** Para reportes, usa `residue` o `combined`. Para estudiar la velocidad del método, usa `absolute_error` y grafica `system.history`. Para un ejercicio donde se piden cinco iteraciones, usa `max_iteration`.

---

## Qué Devuelve un Método Lineal

```mermaid
flowchart TB
    A["Llamada al método"] --> B["Devuelve una tupla"]
    B --> C["solution<br/>vector solución"]
    B --> D["iterations<br/>número de vueltas"]
    B --> E["error<br/>error final"]
```

```python
solution, iterations, error = system.solved_by_gauss_seidel(
    'combined',
    absolute_tol=1e-5,
    residue_tol=1e-5
)

show_solution_set(solution, iterations, 'x', error)
```

> **Historial de iteraciones.** Cuando usas un criterio de parada, el solver guarda información en `system.history`. Esto sirve para hacer tablas o gráficas de convergencia.

```python
for item in system.history:
    print(item['iteration'], item['solution'], item['residue'])
```

---

## Historial de Iteraciones

```mermaid
flowchart TB
    A["Método con criterio de parada"] --> B["system.history<br/>lista de diccionarios"]
    B --> C["Tabla de iteraciones"]
    B --> D["Gráfica de convergencia"]
    B --> E["Diagnóstico del método"]
```

`history` es la memoria de la ejecución. Cada elemento de la lista corresponde a una iteración y guarda la aproximación calculada junto con las medidas de error disponibles.

Se reinicia cada vez que ejecutas de nuevo un método con criterio de parada.

### Qué guarda en sistemas lineales

| Clave | Contenido |
|---|---|
| `iteration` | Número de iteración: \(1,2,3,\ldots\). |
| `solution` | Vector \(x_k\) obtenido en esa iteración. |
| `absolute_error` | \(\|x_k-x_{k-1}\|_2\). |
| `max_absolute_error` | \(\max_i |x_{k,i}-x_{k-1,i}|\). |
| `relative_error` | \(\|x_k-x_{k-1}\|_2/\|x_k\|_2\). |
| `residue` | \(\|Ax_k-b\|\), por defecto con norma euclidiana. |

```python
solution, iterations, error = system.solved_by_gauss_seidel(
    'combined',
    absolute_tol=1e-5,
    residue_tol=1e-5
)

# Primer registro del historial
system.history[0]

# Último registro del historial
system.history[-1]
```

Forma típica de un registro lineal:

```python
{
    'iteration': 5,
    'solution': array([0.750002, 0.833331, 2.083330]),
    'absolute_error': 0.000041,
    'max_absolute_error': 0.000036,
    'relative_error': 0.000017,
    'residue': 0.000098
}
```

### Convertir el historial en tabla

```python
import pandas as pd

history_table = pd.DataFrame(system.history)

# Separar el vector solución en columnas x_1, x_2, ...
solution_columns = pd.DataFrame(
    history_table['solution'].tolist(),
    columns=['x_1', 'x_2', 'x_3']
)

history_table = pd.concat(
    [history_table.drop(columns=['solution']), solution_columns],
    axis=1
)

history_table.round(6)
```

| iter | \(x_1\) | \(x_2\) | \(x_3\) | err. abs. | residuo |
|---:|---:|---:|---:|---:|---:|
| 1 | 2.000000 | 1.000000 | -1.500000 | 2.692582 | 14.500000 |
| 2 | -1.000000 | 1.900000 | 4.550000 | 6.808818 | 28.750000 |
| 3 | -7.050000 | 1.900000 | 16.200000 | 13.125738 | 70.250000 |
| ... | ... | ... | ... | ... | ... |

> **Lectura de la tabla.** Si el error y el residuo bajan, el método está convergiendo. Si crecen o empiezan a oscilar, cambia el orden de las ecuaciones, revisa la dominancia diagonal o usa el método directo como comparación.

### Graficar el historial

```python
import matplotlib.pyplot as plt

errors = [item['absolute_error'] for item in system.history]
residues = [item['residue'] for item in system.history]
iterations = range(1, len(system.history) + 1)

plt.semilogy(iterations, errors, marker='o', label='Error absoluto')
plt.semilogy(iterations, residues, marker='s', label='Residuo')
plt.xlabel('Iteración')
plt.ylabel('Magnitud')
plt.grid(True, alpha=0.3)
plt.legend()
plt.show()
```

---

## Flujo Para Sistemas No Lineales

```mermaid
flowchart TB
    A["1. Definir F(x,y)"] --> B["2. Definir Jacobiano"]
    B --> C["3. Elegir punto inicial"]
    A --> D["4. Crear NonLinearSystem"]
    B --> D
    C --> D
    D --> E["5A. Newton-Raphson"]
    D --> F["5B. Punto fijo"]
    E --> G["6. Solución, iteraciones, error"]
    F --> G
```

### Ejemplo Newton-Raphson

Sistema:

$$
\begin{cases}
-x^2+x+0.75-y=0,\\
x^2-5xy-y=0.
\end{cases}
$$

```python
import numpy as np
from LNL_Solver import *

functions = [
    lambda x, y: -x**2 + x + 0.75 - y,
    lambda x, y: x**2 - 5*x*y - y
]

jacobian = [
    [lambda x, y: -2*x + 1, lambda x, y: -1],
    [lambda x, y: 2*x - 5*y, lambda x, y: -5*x - 1]
]

system = NonLinearSystem(functions, [1.2, 1.2])

solution, iterations, error = system.solved_by_newton_raphson(
    jacobian,
    tol=0.3/100
)

show_solution_set(solution, iterations, ['x', 'y'], error)
```

---

## Newton-Raphson

```mermaid
flowchart TB
    A["Punto actual<br/>(x_k,y_k)"] --> B["Evaluar F(x_k,y_k)"]
    B --> C["Evaluar J(x_k,y_k)"]
    C --> D["Actualizar<br/>x_{k+1}=x_k-J(x_k)^{-1}F(x_k)"]
    B --> D
    D --> E{"¿Error menor que tolerancia?"}
    E -->|Sí| F["Devolver solución"]
    E -->|No| A
```

Newton-Raphson suele converger en pocas iteraciones, pero necesita Jacobiano invertible. Si aparece un error de matriz singular, cambia el punto inicial o revisa el Jacobiano.

### Criterios no lineales

| Criterio | Interpretación | Parámetro |
|---|---|---|
| `relative_error` | Error relativo aproximado entre iteraciones. | `tol` |
| `max_absolute_error` | Máximo cambio por variable. | `tol` |
| `absolute_error` | Norma euclidiana del cambio. | `tol` |
| `euclidean_error` | Alias del error absoluto. | `tol` |
| `residue` | Norma de \(F(x,y)\). | `tol` |
| `combined` | Combina máximo cambio y residuo. | dos tolerancias |

En Newton-Raphson, `relative_error`, `absolute_error` y `max_absolute_error` miran el cambio entre dos aproximaciones. En cambio, `residue` mira qué tan cerca está el punto de cumplir \(F(x,y)=0\). Si necesitas una verificación fuerte, usa `combined`.

```python
# El criterio predeterminado es relative_error
system.solved_by_newton_raphson(jacobian, tol=1e-5)

# Criterios alternativos
system.solved_by_newton_raphson(jacobian, stop_criteria='max_absolute_error', tol=1e-6)
system.solved_by_newton_raphson(jacobian, stop_criteria='absolute_error', tol=1e-6)
system.solved_by_newton_raphson(jacobian, stop_criteria='residue', tol=1e-8)
system.solved_by_newton_raphson(
    jacobian,
    stop_criteria='combined',
    absolute_tol=1e-6,
    residue_tol=1e-8
)
```

### Historial no lineal

| Método | Qué queda guardado en `history` |
|---|---|
| Newton-Raphson | `iteration`, `solution`, `error`, `relative_error`, `max_absolute_error`, `absolute_error`, `residue`. |
| Punto fijo | `iteration`, `solution`, `error`. |

```python
solution, iterations, error = system.solved_by_newton_raphson(
    jacobian,
    stop_criteria='residue',
    tol=1e-8
)

for row in system.history:
    print(row['iteration'], row['solution'], row['residue'])
```

Registro típico de Newton-Raphson:

```python
{
    'iteration': 5,
    'solution': array([1.372065, 0.239502]),
    'error': 0.000000,
    'relative_error': 0.000009,
    'max_absolute_error': 0.000000,
    'absolute_error': 0.000000,
    'residue': 0.000000
}
```

---

## Punto Fijo en Sistemas No Lineales

```mermaid
flowchart LR
    A["Sistema original<br/>F(x,y)=0"] --> B["Despejar<br/>x=g1(x,y), y=g2(x,y)"]
    B --> C["Iterar<br/>con valores actualizados"]
```

```python
fixed_point_functions = [
    lambda x, y: np.sqrt(-y + x + 0.75),
    lambda x, y: (x**2) / (1 + 5*x)
]

system = NonLinearSystem(functions, [1.2, 1.2])

solution, iterations, error = system.solved_by_fixed_point(
    fixed_point_functions,
    tol=1e-4,
    max_iteration=100
)

show_solution_set(solution, iterations, ['x', 'y'], error)
```

> **Punto fijo depende muchísimo del despeje.** El mismo sistema puede converger o divergir dependiendo de cómo se despejen las ecuaciones. Si no converge, prueba otra formulación \(g_1,g_2\) o usa Newton-Raphson.

---

## Cómo Leer `show_solution_set`

```mermaid
flowchart LR
    A["Entrada<br/>solution, iterations, error"] --> B["Da formato<br/>a cada variable"]
    B --> C["Imprime<br/>resultado final"]
```

```python
# Para variables x_1, x_2, x_3...
show_solution_set(solution, iterations, 'x', error)

# Para variables con nombre propio
show_solution_set(solution, iterations, ['x', 'y'], error)

# Para poner un título al resultado
show_solution_set(
    solution,
    iterations,
    ['x', 'y'],
    error,
    system_name='Usando Newton-Raphson'
)
```

Forma de la salida:

```text
Conjunto solución | Usando Newton-Raphson
x = 1.372065
y = 0.239502
Resultados obtenidos a 5 iteraciones
Error = 8.935694e-06 %
```

---

## Diagnóstico Rápido de Errores

| Síntoma | Causa probable | Solución |
|---|---|---|
| `ModuleNotFoundError` | El notebook no ve `LNL_Solver.py`. | Coloca el archivo en la misma carpeta o usa `sys.path.append`. |
| `from LNL_Solver.py import *` falla | Se incluyó la extensión `.py`. | Usa `from LNL_Solver import *`. |
| Matriz singular en Newton | Jacobiano no invertible en el punto actual. | Cambia el punto inicial o revisa derivadas. |
| Gauss-Seidel no converge | La matriz no garantiza convergencia. | Reordena ecuaciones, pivotea o usa otro método. |
| Cero en diagonal | División por cero en método iterativo. | Reordena filas o reformula el sistema. |
| Resultados viejos tras editar el solver | Jupyter ya había cargado el módulo. | Reinicia kernel o usa `reload`. |

```mermaid
flowchart LR
    A["No<br/>from LNL_Solver.py import *"] --> B["Sí<br/>from LNL_Solver import *"]
```

---

## Recetas de Uso

### Receta A: resolver \(Ax=b\) con Gauss-Seidel

```python
import numpy as np
from LNL_Solver import *

A = np.array([
    [4, -1, 0],
    [-1, 4, -1],
    [0, -1, 4]
])

b = np.array([[15], [10], [10]])
x0 = [0, 0, 0]

system = LinearSystem(A, b, x0)
system.convergence()

solution, iterations, error = system.solved_by_gauss_seidel(
    'combined',
    absolute_tol=1e-6,
    residue_tol=1e-8
)

show_solution_set(solution, iterations, 'x', error)
```

### Receta B: comparar Gauss-Seidel con método directo

```python
iter_solution, iter_n, iter_error = system.solved_by_gauss_seidel(
    'residue',
    tol=1e-8
)

direct_solution, direct_n, direct_error = system.solved_by_matrix_direct_method()

show_solution_set(iter_solution, iter_n, 'x', iter_error)
show_solution_set(direct_solution, direct_n, 'x', direct_error)
```

### Receta C: graficar el historial de convergencia

```python
import matplotlib.pyplot as plt

errors = [item['residue'] for item in system.history]

plt.semilogy(range(1, len(errors) + 1), errors, marker='o')
plt.xlabel('Iteración')
plt.ylabel('Residuo')
plt.grid(True)
plt.show()
```

---

## Ejercicio Guiado: Solver + `numpy` + `matplotlib`

El objetivo es resolver y visualizar un sistema no lineal completo usando `LNL_Solver.py`. La idea es ver el flujo real de trabajo: importar librerías, definir funciones con `numpy`, resolver con Newton-Raphson, usar el historial y producir gráficas con `matplotlib`.

Usaremos el sistema:

$$
x^2+y^2=4,
\qquad
x^2-y=1.
$$

En forma \(F(x,y)=0\):

$$
f_1(x,y)=x^2+y^2-4,
\qquad
f_2(x,y)=x^2-y-1.
$$

```mermaid
flowchart LR
    A["numpy<br/>mallas y arreglos"] --> B["LNL-Solver<br/>LNL_Solver.py"]
    B --> C["matplotlib<br/>curvas e historial"]
```

### Paso 1: importar el módulo y las librerías

```python
import numpy as np
import matplotlib.pyplot as plt
from LNL_Solver import *
```

> **Importación correcta.** El import correcto es `from LNL_Solver import *`. No se escribe la extensión `.py` dentro del import. El nombre visual del proyecto puede llevar guion, pero el módulo de Python usa guion bajo.

### Paso 2: definir \(F(x,y)\) y el Jacobiano

```python
functions = [
    lambda x, y: x**2 + y**2 - 4,
    lambda x, y: x**2 - y - 1
]

jacobian = [
    [lambda x, y: 2*x, lambda x, y: 2*y],
    [lambda x, y: 2*x, lambda x, y: -1]
]
```

Estas funciones sirven tanto para resolver como para graficar porque usan operaciones compatibles con `numpy`. Pueden recibir números individuales o arreglos creados con `np.meshgrid`.

### Paso 3: resolver desde un punto inicial

```python
initial_point = [1, 1]
system = NonLinearSystem(functions, initial_point)

solution, iterations, residue = system.solved_by_newton_raphson(
    jacobian,
    stop_criteria='residue',
    tol=1e-8,
    max_iteration=30
)

print(solution)
print(iterations)
print(f"Residuo final = {residue:.3e}")
```

Salida esperada:

```text
[1.517490 1.302776]
4
Residuo final = 3.166e-10
```

El punto \((1.517490,\ 1.302776)\) está en la intersección de la circunferencia y la parábola. El residuo pequeño confirma que \(F(x,y)\) está muy cerca de cero.

### Paso 4: probar varios puntos iniciales

```python
initial_points = [
    [1, 1],
    [-1, 1],
    [2, 0],
    [-0.2, 0.2],
    [0, 2]
]

results = []

for point in initial_points:
    try:
        system = NonLinearSystem(functions, point)
        solution, iterations, residue = system.solved_by_newton_raphson(
            jacobian,
            stop_criteria='residue',
            tol=1e-8,
            max_iteration=30
        )
        results.append([point, solution, iterations, residue])
    except np.linalg.LinAlgError:
        results.append([point, 'Jacobiano singular', None, None])
```

| Punto inicial | Solución | Iter. | Residuo |
|---|---|---:|---:|
| `[1,1]` | `(1.517490, 1.302776)` | 4 | \(3.166\times10^{-10}\) |
| `[-1,1]` | `(-1.517490, 1.302776)` | 4 | \(3.166\times10^{-10}\) |
| `[2,0]` | `(1.517490, 1.302776)` | 6 | \(2.809\times10^{-15}\) |
| `[-0.2,0.2]` | `(-1.517490, 1.302776)` | 7 | \(3.225\times10^{-10}\) |
| `[0,2]` | Jacobiano singular | -- | -- |

> **Jacobiano singular.** Newton-Raphson necesita invertir \(J(x,y)\). Si el punto inicial produce un Jacobiano no invertible, el paso \(J^{-1}F\) no existe. En ese caso, se prueba otro punto inicial o se revisa la formulación del sistema.

### Paso 5: graficar las curvas y las trayectorias

```python
x = np.linspace(-2.5, 2.5, 500)
y = np.linspace(-2.5, 2.5, 500)
X, Y = np.meshgrid(x, y)

plt.contour(X, Y, functions[0](X, Y), levels=[0])
plt.contour(X, Y, functions[1](X, Y), levels=[0])
plt.xlabel('x')
plt.ylabel('y')
plt.axis('equal')
plt.grid(True)
plt.show()
```

Para dibujar las trayectorias de Newton, se usa el historial:

```python
path = np.vstack([
    np.array(initial_point),
    *[row['solution'] for row in system.history]
])

plt.plot(path[:, 0], path[:, 1], '--o')
```

![Curvas del sistema, puntos iniciales y trayectorias](figures_whsolver_guide/newton_paths.png)

### Paso 6: usar `history` para estudiar convergencia

Después de ejecutar el método, `system.history` contiene un diccionario por iteración:

```python
history = system.history

for row in history:
    print(
        row['iteration'],
        row['solution'],
        row['absolute_error'],
        row['residue']
    )
```

También se puede convertir en tabla:

```python
import pandas as pd

history_table = pd.DataFrame(system.history)
history_table[['iteration', 'absolute_error', 'relative_error', 'residue']]
```

Y graficar el comportamiento de los errores:

```python
iterations = [row['iteration'] for row in system.history]
absolute_errors = [row['absolute_error'] for row in system.history]
relative_errors = [row['relative_error'] for row in system.history]
residues = [row['residue'] for row in system.history]

plt.semilogy(iterations, absolute_errors, marker='o', label='Error absoluto')
plt.semilogy(iterations, relative_errors, marker='s', label='Error relativo')
plt.semilogy(iterations, residues, marker='^', label='Residuo')
plt.xlabel('Iteración')
plt.ylabel('Magnitud')
plt.grid(True, which='both')
plt.legend()
plt.show()
```

![Historial de errores y residuo](figures_whsolver_guide/newton_history.png)

La escala logarítmica permite ver caídas muy grandes en pocas iteraciones. Si el residuo baja rápidamente, el punto calculado satisface cada vez mejor el sistema original.

### Paso 7: mapa de convergencia por punto inicial

Un uso más avanzado combina `numpy`, `matplotlib` y el solver: generar muchos puntos iniciales y colorearlos según la raíz a la que llegan.

```python
x0 = np.linspace(-2, 2, 85)
y0 = np.linspace(-2, 2, 85)
X0, Y0 = np.meshgrid(x0, y0)
points = np.stack((X0, Y0), axis=-1).reshape(-1, 2)

roots = []
root_number = []

for point in points:
    try:
        system = NonLinearSystem(functions, point)
        solution, iterations, residue = system.solved_by_newton_raphson(
            jacobian,
            stop_criteria='residue',
            tol=1e-7,
            max_iteration=18
        )

        # Clasificar la solución encontrada.
        # Si el residuo no es pequeño, marcar como no convergente.

    except np.linalg.LinAlgError:
        root_number.append(-1)
```

![Mapa de convergencia por punto inicial](figures_whsolver_guide/newton_basins.png)

> **Qué enseña el mapa.** Newton-Raphson es local: el punto inicial importa. Una zona verde y una zona violeta indican regiones que convergen a raíces distintas. Las marcas negras señalan puntos donde el método no converge con la tolerancia o el número máximo de iteraciones elegidos.

---

## Checklist Antes de Entregar

```mermaid
flowchart TB
    A["1. LNL_Solver.py está en la carpeta correcta"] --> B["2. Importaste con<br/>from LNL_Solver import *"]
    B --> C["3. Definiste A,b,x0<br/>o F,J,p0"]
    B --> D["4. Elegiste método<br/>y criterio de parada"]
    E["6. Revisaste convergencia<br/>o residuo"] --> D
    D --> F["5. Mostraste solución,<br/>iteraciones y error"]
```

> **Regla de oro.** Primero formula bien el problema matemático. Después codifica matrices, funciones y derivadas. El solver automatiza las iteraciones, pero la calidad del resultado depende de que el sistema esté bien planteado.

---

## Resumen Visual Final

```mermaid
flowchart LR
    A["Importar<br/>from LNL_Solver import *"] --> B["Elegir<br/>lineal o no lineal"]
    B --> C["Construir<br/>objeto del sistema"]
    C --> D["Resolver<br/>con método y criterio"]
    D --> E["Presentar<br/>solución, iteraciones, error"]
```
