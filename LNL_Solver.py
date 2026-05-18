import numpy as np
class LinearSystem:
    """Representa un sistema lineal de la forma Ax = b.

    La clase permite resolver el sistema con métodos iterativos como
    Gauss-Seidel y Jacobi, o con el método directo de la matriz inversa.
    Los métodos iterativos pueden usar criterios de parada como número
    máximo de iteraciones, error absoluto, error relativo, residuo,
    criterio combinado y error absoluto máximo.

    Parámetros:
        coeff_matrix: matriz de coeficientes A.
        indep_terms_matrix: vector o matriz columna de términos independientes b.
        initial_point: aproximación inicial para los métodos iterativos.

    Estructura para usar criterios de parada:
        sistema.solved_by_gauss_seidel(stop_criteria='max_iteration', max_iteration=10)
        sistema.solved_by_jacobi(stop_criteria='absolute_error', tol=1e-5)
        sistema.solved_by_jacobi(stop_criteria='relative_error', tol=1e-5)
        sistema.solved_by_gauss_seidel(stop_criteria='residue', tol=1e-5)
        sistema.solved_by_jacobi(stop_criteria='max_absolute_error', tol=1e-5)
        sistema.solved_by_gauss_seidel(stop_criteria='combined', absolute_tol=1e-5, residue_tol=1e-5)

    Criterios disponibles:
        'max_iteration': detiene el método al llegar al máximo de iteraciones.
        'absolute_error': usa la norma euclidiana del error entre iteraciones consecutivas.
        'relative_error': usa el error relativo entre iteraciones consecutivas.
        'residue': usa la norma del residuo Ax - b.
        'max_absolute_error': usa la norma supremo del error.
        'combined': combina error absoluto máximo y residuo.

    Nota:
        Gauss-Seidel y Jacobi requieren una aproximación inicial. El método
        directo no es iterativo y no usa tolerancia ni número de iteraciones.
    """
    
    def __init__(self, coeff_matrix, indep_terms_matrix, initial_point) -> None:
        self.coeff_matrix = coeff_matrix
        self.indep_terms_matrix = indep_terms_matrix
        self.initial_point = np.array(initial_point, dtype=float).reshape(-1)
        self.n = 0
        self.error = None
    
    # Verificamos si la matriz es diagonalmente dominante
    def is_diagonally_dominant(self, matrix = None):
        coeff_matrix = self.coeff_matrix
        matrix = coeff_matrix
        true_false_list = list()
        for i,j in enumerate(range(len(matrix))):
            abs_row_elements = list()
            
            d_element = matrix[i][j]
            for k in range(len(matrix)):
                if i == k:
                    continue
                else:
                    abs_row_elements.append(abs(matrix[i][k]))
            
            true_false_list.append(abs(d_element) > sum(abs_row_elements))

        if all(true_false_list):
            return True
        else:
            return False

    # Verifico si es convergente el método
    def convergence(self):
        # matrix = coeff_matrix
        if self.is_diagonally_dominant():
            print("La matriz es diagonalmente dominante, se garantiza la convergencia")
        else:
            print("** NOTA IMPORTANTE ** \n La matriz no es diagonalmente dominante, no se garantiza la convergencia")

    # Método iterativo
    def solved_by_gauss_seidel(self):
        coeff_matrix = np.array(self.coeff_matrix, dtype=float)
        indep_terms_matrix = np.array(self.indep_terms_matrix, dtype=float).reshape(-1)
        initial_point = np.array(self.initial_point, dtype=float).reshape(-1)
        self.initial_point = initial_point
        
        if np.any(np.isclose(np.diag(coeff_matrix), 0)):
            raise ValueError('No se puede aplicar Gauss-Seidel porque hay ceros en la diagonal')
        for j in range(0, len(coeff_matrix)):
            summ_val = indep_terms_matrix[j]
            for i in range(0, len(coeff_matrix)):
                if (j!=i):
                    summ_val = summ_val - coeff_matrix[j][i] * initial_point[i]
                    # print('-----')
                    # print(f"Row {j+1}, Column {i+1}")
                    # print(summ_val, coeff_matrix[j][i], initial_point[i])
                    # print('-----')
                    # print()
                self.initial_point[j] = summ_val/coeff_matrix[j][j]
        self.n += 1
        # return initial_point
        return self.initial_point
    
    def solved_by_jacobi(self):
        coeff_matrix = np.array(self.coeff_matrix, dtype=float)
        indep_terms_matrix = np.array(self.indep_terms_matrix, dtype=float).reshape(-1)
        initial_point = np.array(self.initial_point, dtype=float).reshape(-1)
        new_point = initial_point.copy()
        for j in range(0, len(coeff_matrix)):
            summ_val = indep_terms_matrix[j]
            for i in range(0, len(coeff_matrix)):
                if j != i:
                    summ_val = summ_val - coeff_matrix[j][i] * initial_point[i]
            new_point[j] = summ_val / coeff_matrix[j][j]
        self.initial_point = new_point
        self.n += 1
        return self.initial_point
    
    def solved_by_matrix_direct_method(self, *args, **kwargs):
        self.n = 0
        self.error = 0

        if args or kwargs:
            raise ValueError('El método directo no es iterativo; no se debe usar stop_criteria, tol ni max_iteration')

        try:
            coeff_matrix = np.array(self.coeff_matrix, dtype=float)
            indep_terms_matrix = np.array(self.indep_terms_matrix, dtype=float)
            det_coeff_matrix = np.linalg.det(coeff_matrix)
            if np.isclose(det_coeff_matrix, 0):
                return 'El determinante es 0, el sistema tiene infinitas soluciones o no tiene', self.n, self.error

            inv_coeff_matrix = np.linalg.inv(coeff_matrix)
            self.initial_point = inv_coeff_matrix @ indep_terms_matrix

        except np.linalg.LinAlgError:
            return f"La matriz no es invertible", self.n, self.error
        
        except ValueError:
            return f"La matriz no es invertible, revisa que sea cuadrada", self.n, self.error
        
        return self.initial_point, self.n, self.error


    def __getattribute__(self, name):
        method = object.__getattribute__(self, name)

        if not (name.startswith('solved_by_') and callable(method)) or name == 'solved_by_matrix_direct_method':
            return method

        def method_with_stop_criteria(*args, **kwargs):
            if not args and not kwargs:
                return method()

            return self._solve_with_stop_criteria(method, *args, **kwargs)

        return method_with_stop_criteria

    def _solve_with_stop_criteria(self, method, stop_criteria='max_iteration', max_iteration=100, tol=1e-5, absolute_tol=1e-5, residue_tol=1e-5, residue_norm=2):
        criterias = ['max_iteration', 'absolute_error', 'relative_error', 'residue', 'combined', 'max_absolute_error']
        if stop_criteria not in criterias:
            raise ValueError(f"El criterio de parada '{stop_criteria}' no está soportado")

        A = np.asarray(self.coeff_matrix, dtype=float)
        b = np.asarray(self.indep_terms_matrix, dtype=float).reshape(-1)
        self.history = []
        self.n = 0
        self.error = None

        for iteration in range(1, max_iteration + 1):
            previous = np.asarray(self.initial_point, dtype=float).reshape(-1)
            current = np.asarray(method(), dtype=float).reshape(-1)

            absolute_error = np.linalg.norm(current - previous)
            max_absolute_error = np.max(np.abs(current - previous))
            current_norm = np.linalg.norm(current)
            relative_error = absolute_error / current_norm if current_norm != 0 else absolute_error
            residue = np.linalg.norm(A @ current - b, ord=residue_norm)

            self.history.append({
                'iteration': iteration,
                'solution': current.copy(),
                'absolute_error': absolute_error,
                'max_absolute_error': max_absolute_error,
                'relative_error': relative_error,
                'residue': residue,
            })

            stop = False
            stop = stop or (stop_criteria == 'max_iteration' and iteration >= max_iteration)
            stop = stop or (stop_criteria == 'absolute_error' and absolute_error < tol)
            stop = stop or (stop_criteria == 'max_absolute_error' and max_absolute_error < tol)
            stop = stop or (stop_criteria == 'relative_error' and relative_error < tol)
            stop = stop or (stop_criteria == 'residue' and residue < tol)
            stop = stop or (stop_criteria == 'combined' and max_absolute_error < absolute_tol and residue < residue_tol)

            if stop or iteration >= max_iteration:
                if stop_criteria == 'max_iteration':
                    self.error = None
                elif stop_criteria == 'absolute_error':
                    self.error = absolute_error
                elif stop_criteria == 'max_absolute_error':
                    self.error = max_absolute_error
                elif stop_criteria == 'relative_error':
                    self.error = relative_error
                elif stop_criteria == 'residue':
                    self.error = residue
                elif stop_criteria == 'combined':
                    self.error = {'max_absolute_error': max_absolute_error, 'residue': residue}
                return current, self.n, self.error

        return current, self.n, self.error

class NonLinearSystem:
    """Sistema no lineal de dos variables.

    functions: lista con las funciones del sistema, escritas como funciones de x, y.
        Ejemplo: [lambda x, y: f1, lambda x, y: f2]

    initial_point: punto inicial de la forma [x_0, y_0].
        Ejemplo: [1.2, 1.2]

    jacobian: lista de listas con las derivadas parciales.
        Cada derivada también debe escribirse como lambda x, y: ...
        Ejemplo: [[lambda x, y: df1_dx, lambda x, y: df1_dy], [lambda x, y: df2_dx, lambda x, y: df2_dy]]

    fixed_point_functions: lista con las funciones despejadas para punto fijo.
        Se usan los valores más actualizados. Por ejemplo, para y se usa x nuevo.
        Ejemplo: [lambda x, y: g1, lambda x, y: g2]

    Clase:
        NonLinearSystem(functions, initial_point)
            Crea el sistema no lineal.
            Ejemplo: system = NonLinearSystem(functions, [1.2, 1.2])    
    
    Métodos:

        evaluate_functions(point)
            Evalúa las funciones del sistema en un punto [x, y].
            Ejemplo: system.evaluate_functions([1.2, 1.2])

        evaluate_jacobian(jacobian, point)
            Evalúa la matriz Jacobiana en un punto [x, y].
            Ejemplo: system.evaluate_jacobian(jacobian, [1.2, 1.2])

        calculate_error(current_point, previous_point)
            Calcula el error relativo aproximado entre dos iteraciones.
            Ejemplo: system.calculate_error([1.3, 0.2], [1.2, 1.2])

        solved_by_fixed_point(fixed_point_functions, tol=1e-4, max_iteration=100)
            Resuelve el sistema usando punto fijo con los valores más actualizados.
            Ejemplo: solution, iterations, error = system.solved_by_fixed_point(fixed_point_functions, tol=1e-4)

        solved_by_newton_raphson(jacobian, tol=1e-4, max_iteration=100, stop_criteria='relative_error')
            Resuelve el sistema usando el método de Newton-Raphson.
            Criterios: 'relative_error', 'max_absolute_error', 'absolute_error', 'euclidean_error', 'residue', 'combined'.
            Ejemplo: solution, iterations, error = system.solved_by_newton_raphson(jacobian, stop_criteria='residue', tol=1e-8)
    """
    def __init__(self, functions, initial_point):
        self.functions = functions
        self.initial_point = np.array(initial_point, dtype=float)
        self.n = 0
        self.error = None
        self.history = []

    def evaluate_functions(self, point):
        x = point[0]
        y = point[1]
        values = []

        for function in self.functions:
            values.append(function(x, y))

        return np.array(values, dtype=float)

    def evaluate_jacobian(self, jacobian, point):
        x = point[0]
        y = point[1]
        values = []

        for row in jacobian:
            jacobian_row = []

            for function in row:
                jacobian_row.append(function(x, y))

            values.append(jacobian_row)

        return np.array(values, dtype=float)

    def calculate_error(self, current_point, previous_point):
        absolute_error = np.linalg.norm(current_point - previous_point)
        current_norm = np.linalg.norm(current_point)

        if current_norm != 0:
            return (absolute_error / current_norm) * 100
        else:
            return absolute_error

    def solved_by_fixed_point(self, fixed_point_functions, tol=1e-4, max_iteration=100):
        current_point = np.array(self.initial_point, dtype=float)
        self.history = []
        self.n = 0
        self.error = None

        for iteration in range(1, max_iteration + 1):
            previous_point = current_point.copy()
            new_point = previous_point.copy()
            for position, function in enumerate(fixed_point_functions):
                x = new_point[0]
                y = new_point[1]
                new_point[position] = function(x, y)

            current_point = np.array(new_point, dtype=float)
            self.error = self.calculate_error(current_point, previous_point)
            self.n = iteration

            self.history.append({
                'iteration': iteration,
                'solution': current_point.copy(),
                'error': self.error
            })
            
            if self.error < tol:
                self.initial_point = current_point
                return current_point, self.n, self.error

        self.initial_point = current_point
        return current_point, self.n, self.error

    def solved_by_newton_raphson(self, jacobian : object, tol=1e-4, max_iteration=100, stop_criteria='relative_error', absolute_tol=1e-6, residue_tol=1e-8)  -> tuple[list, int, float]:
        criterias = ['relative_error', 'max_absolute_error', 'absolute_error', 'euclidean_error', 'residue', 'combined']
        if stop_criteria not in criterias:
            raise ValueError(f"El criterio de parada '{stop_criteria}' no está soportado")

        current_point = np.array(self.initial_point, dtype=float)
        self.history = []
        self.n = 0
        self.error = None

        for iteration in range(1, max_iteration + 1):
            previous_point = current_point.copy()
            function_values = self.evaluate_functions(previous_point)
            jacobian_values = self.evaluate_jacobian(jacobian, previous_point)
            inverse_jacobian = np.linalg.inv(jacobian_values)

            current_point = previous_point - inverse_jacobian @ function_values
            
            difference = current_point - previous_point
            
            relative_error = self.calculate_error(current_point, previous_point)
            max_absolute_error = np.max(np.abs(difference))
            absolute_error = np.linalg.norm(difference)
            residue = np.linalg.norm(self.evaluate_functions(current_point))
            
            self.n = iteration

            if stop_criteria == 'relative_error':
                self.error = relative_error
                stop = relative_error < tol
            elif stop_criteria == 'max_absolute_error':
                self.error = max_absolute_error
                stop = max_absolute_error < tol
            elif stop_criteria == 'absolute_error' or stop_criteria == 'euclidean_error':
                self.error = absolute_error
                stop = absolute_error < tol
            elif stop_criteria == 'residue':
                self.error = residue
                stop = residue < tol
            elif stop_criteria == 'combined':
                self.error = {'max_absolute_error': max_absolute_error, 'residue': residue}
                stop = max_absolute_error < absolute_tol and residue < residue_tol

            self.history.append({
                'iteration': iteration,
                'solution': current_point.copy(),
                'error': self.error,
                'relative_error': relative_error,
                'max_absolute_error': max_absolute_error,
                'absolute_error': absolute_error,
                'residue': residue
            })

            if stop:
                self.initial_point = current_point
                return current_point, self.n, self.error

        self.initial_point = current_point
        return current_point, self.n, self.error

# Imprimir el conjunto solutión
def show_solution_set(solution, iterations, variable_name='x', error=None, system_name=''):
    solution = np.array(solution).reshape(-1)
    if system_name == '':
        print(f'\nConjunto solución')
    else:
        print(f'\nConjunto solución | {system_name}')

    for n,sol in enumerate(solution):
        if isinstance(variable_name, list) and n < len(variable_name):
            print(f"{variable_name[n]} = {sol}")
        else:
            print(f"{variable_name}_{n+1} = {sol}")
    print(f"Resultados obtenidos a {iterations} iteraciones")
    if error is not None:
        print(f"Error = {error} %")
