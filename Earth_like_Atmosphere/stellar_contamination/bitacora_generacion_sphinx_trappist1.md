# Bitácora de generación de contaminación estelar con SPHINX para TRAPPIST-1

## Objetivo

Generar una familia de archivos de contaminación estelar para TRAPPIST-1 usando espectros estelares de la grilla SPHINX, manteniendo la misma convención de coberturas superficiales ya presente en esta carpeta para los archivos basados en PHOENIX.

Los productos finales son los archivos:

- `sphinx_TRAPPIST-1_contam_fspot0.01_ffac0.08.txt`
- `sphinx_TRAPPIST-1_contam_fspot0.01_ffac0.54.txt`
- `sphinx_TRAPPIST-1_contam_fspot0.01_ffac0.70.txt`
- `sphinx_TRAPPIST-1_contam_fspot0.08_ffac0.08.txt`
- `sphinx_TRAPPIST-1_contam_fspot0.08_ffac0.54.txt`
- `sphinx_TRAPPIST-1_contam_fspot0.08_ffac0.70.txt`
- `sphinx_TRAPPIST-1_contam_fspot0.26_ffac0.08.txt`
- `sphinx_TRAPPIST-1_contam_fspot0.26_ffac0.54.txt`
- `sphinx_TRAPPIST-1_contam_fspot0.26_ffac0.70.txt`

## Organización final

La carpeta quedó separada en insumos y productos:

- [sphinx_data](/C:/Proyectos/Astro/gdaespec/Earth_like_Atmosphere/stellar_contamination/sphinx_data)
  - contiene las grillas SPHINX crudas `Teff_*`
- [generate_sphinx_trappist1_contamination.ipynb](/C:/Proyectos/Astro/gdaespec/Earth_like_Atmosphere/stellar_contamination/generate_sphinx_trappist1_contamination.ipynb)
  - notebook reproducible que genera y plotea las curvas
- [Earth_like_Atmosphere/stellar_contamination](/C:/Proyectos/Astro/gdaespec/Earth_like_Atmosphere/stellar_contamination)
  - contiene los productos `sphinx_TRAPPIST-1_contam_*` junto a los `TRAPPIST-1_contam_*` ya existentes

## Parámetros físicos usados

Se adoptó el caso de TRAPPIST-1 definido en el flujo de trabajo de esta carpeta:

- `T_s = 2566 K`
- `log g = 5.25`
- `[M/H] = 0.0`
- `C/O = 0.5`

Para las heterogeneidades estelares se usó:

- `T_spot = 0.86 T_s`
- `T_fac = T_s + 100 K`

Las fracciones de cobertura superficial fueron las mismas ya utilizadas en esta carpeta:

- `f_spot = {0.01, 0.08, 0.26}`
- `f_fac = {0.08, 0.54, 0.70}`

## Grilla espectral

Las curvas se evaluaron sobre la grilla compartida en [waves.txt](/C:/Proyectos/Astro/gdaespec/Earth_like_Atmosphere/waves.txt), de modo que los productos SPHINX sean directamente comparables con los archivos de contaminación ya existentes en esta carpeta.

## Procedimiento

### 1. Indexado de la grilla SPHINX

El notebook indexa los archivos `Teff_*` de `sphinx_data/` usando como claves:

- `Teff`
- `logg`
- `logZ`
- `C/O`

Esto permite localizar rápidamente el espectro estelar apropiado para la fotosfera quieta, las manchas y las fáculas.

### 2. Lectura y limpieza de los espectros

Cada archivo de SPHINX se trata como un espectro de dos columnas:

- longitud de onda en micrones
- flujo monocromático estelar `F_lambda`

Durante la lectura se conservan únicamente valores:

- finitos
- positivos
- ordenados crecientemente en longitud de onda

### 3. Regrillado espectral

Los espectros SPHINX se interpolan linealmente sobre la grilla común de `waves.txt`.

Este paso asegura que:

- todos los espectros estelares estén definidos en la misma malla
- la construcción de `epsilon(lambda)` sea consistente entre configuraciones
- los archivos finales puedan compararse punto a punto con los productos PHOENIX

### 4. Interpolación en temperatura efectiva

Como la temperatura de fotosfera, manchas o fáculas no siempre coincide exactamente con un nodo de la grilla SPHINX, el notebook interpola linealmente en `Teff` a `logg`, metalicidad y `C/O` fijos.

En la práctica, esto se usa para construir:

- `F_phot(lambda)`
- `F_spot(lambda)`
- `F_fac(lambda)`

### 5. Construcción del espectro integrado del disco estelar

Para cada par `(f_spot, f_fac)`, el flujo del disco durante el tránsito se modela como:

```math
F_{\mathrm{disk}}(\lambda) =
f_{\mathrm{quiet}} F_{\mathrm{phot}}(\lambda) +
f_{\mathrm{spot}} F_{\mathrm{spot}}(\lambda) +
f_{\mathrm{fac}} F_{\mathrm{fac}}(\lambda),
```

con

```math
f_{\mathrm{quiet}} = 1 - f_{\mathrm{spot}} - f_{\mathrm{fac}}.
```

### 6. Cálculo del factor de contaminación estelar

El factor de contaminación se definió como:

```math
\epsilon(\lambda) =
\frac{F_{\mathrm{phot}}(\lambda)}{F_{\mathrm{disk}}(\lambda)}.
```

Esta es la cantidad que multiplica el espectro de transmisión observado cuando se quiere modelar el efecto de una superficie estelar heterogénea.

## Formato de salida

Cada archivo `sphinx_TRAPPIST-1_contam_*` se guarda con dos columnas:

- columna 1: longitud de onda en micrones
- columna 2: `epsilon(lambda)`

Esto sigue la convención práctica usada por los archivos `TRAPPIST-1_contam_*` preexistentes.

## Notebook reproducible

El flujo quedó encapsulado en el notebook:

- [generate_sphinx_trappist1_contamination.ipynb](/C:/Proyectos/Astro/gdaespec/Earth_like_Atmosphere/stellar_contamination/generate_sphinx_trappist1_contamination.ipynb)

Ese notebook:

1. carga las grillas crudas desde `sphinx_data/`
2. genera los nueve archivos `sphinx_TRAPPIST-1_contam_*`
3. al final plotea todas las curvas `epsilon(lambda)` generadas

## Procedencia de las grillas SPHINX

Las grillas locales `Teff_*` fueron descargadas por el usuario desde Zenodo y corresponden a modelos SPHINX para enanas M.

Referencia de Zenodo usada como fuente del grid SPHINX:

- Iyer, R. Aishwarya, Line, R. Michael, Muirhead, S. Philip, Fortney, J. Jonathan, & Gharib-Nezhad, Ehsan. *The SPHINX M-dwarf Spectral Grid. I. Benchmarking New Model Atmospheres to Derive Fundamental M-Dwarf Properties*. Zenodo, versión v2, DOI: `10.5281/zenodo.7416042`.
- Enlace: [https://zenodo.org/records/7416042](https://zenodo.org/records/7416042)

## Nota de trazabilidad

Esta bitácora documenta el flujo realmente implementado en el repositorio tras la reorganización de la carpeta:

- separación de insumos en `sphinx_data/`
- generación reproducible mediante notebook
- productos finales guardados junto al resto de archivos de contaminación estelar

Si en el futuro se reemplaza la submalla local `Teff_*` por una versión distinta del grid SPHINX, conviene actualizar esta bitácora indicando:

- la nueva referencia
- el rango de parámetros descargado
- la fecha de actualización
