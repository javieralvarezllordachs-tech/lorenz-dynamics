# lorenz-dynamics

Código de mi Trabajo de Fin de Grado (Grado en Física, Universidad Complutense de Madrid) sobre la predictibilidad del sistema de Lorenz.

**Estado:** en desarrollo (versión del 2026-10-28).

## Qué hay

`codigo/rk4.py`:

- `lorenz(x, r)`: el sistema de Lorenz, con sigma = 10 y b = 8/3 fijos y r = 28 por defecto.
- `paso_euler(f, x, h)` y `paso_rk4(f, x, h)`: un paso de Euler y uno de Runge–Kutta 4.
- `integrador(g, f, x0, T, n)`: n pasos de tamaño fijo h = T/n.
- `errores`, `pendiente` y `pendienteslocales`: el error en el instante final para varios pasos y su pendiente en escala log-log, que da el orden del método.

## Qué falta

- La comprobación del orden en Lorenz: tabla de pendientes y figura `orden.png`.

## Cómo ejecutarlo

Con Python 3.11, en la carpeta del repositorio:

    python -m venv .venv
    .venv\Scripts\Activate.ps1
    python -m pip install -r requirements.txt
    python codigo/rk4.py

En Linux o macOS, el entorno se activa con `source .venv/bin/activate`.
