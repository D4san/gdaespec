# Synthetic Transmission Spectra Corpus (General, Instrument-Agnostic)

## Propósito
Entrenar y evaluar un **descontaminador** (IA) de espectros de **transmisión** bajo **ruido alto** y posible **contaminación estelar**, cubriendo dos familias atmosféricas:
- **Primarias**: gas de relleno (**fill gas**) **H₂–He**  
- **Secundarias**: fill gas **N₂**

El inventario molecular activo del set base es reducido pero informativo: **H₂O, CO₂, CH₄, CO**.  
Además, **cada molécula** admite el estado **“ausente”** (**NaN**) y se incluye un **escenario sin atmósfera** (espectro plano).

---

## Rango espectral e instrumento
El **forward model** se genera en continuo (p. ej., 0.4–12 μm) y luego se aplica **convolución/remuestreo** a los perfiles instrumentales que se necesiten (R bajo/medio/alto).  
Este README es **agnóstico del instrumento**: define la física base y la grilla de parámetros.

---

## Física y grillas

### 1) Perfil térmico (isotermo)
- **Temperaturas**: $T \in \{250, 300, 350, 400\}\,\text{K}$ (paso de 50 K)

> Justificación: en transmisión a presiones mbar, la dependencia fina en $T(p)$ suele ser débil; el perfil isotermo es una primera aproximación estándar.

### 2) Presión de referencia
- **$P_{\rm ref}$** (bar): $\{10,\;5,\;1,\;0.5,\;0.1\}$

> $P_{\rm ref}$ fija el anclaje geométrico y, junto con $T$, la masa molecular media $\mu$ y la gravedad superficial $g$, modula la amplitud efectiva. La escala de altura obedece $H \propto T/(\mu\, g)$.

### 3) Familias composicionales
- **primary\_H2**: fill gas **H₂–He**; especies activas: **H₂O, CO₂, CH₄, CO**  
  - Continuos siempre activos: **Rayleigh (H₂/He)** y **CIA H₂–H₂ / H₂–He**.
- **secondary\_N2**: fill gas **N₂**; especies activas: **H₂O, CO₂, CH₄, CO**  
  - Continuos siempre activos: **Rayleigh (N₂)** y **CIA N₂–N₂**.  
- **flat\_noatm**: caso **sin atmósfera** (espectro plano; sin líneas ni continuos).

> **Importante**: No se incluye O₂ ni ninguna CIA asociada a O₂.

### 4) Abundancias (VMR) por especie activa
Para **cada especie** $X \in \{\text{H}_2\text{O},\ \text{CO}_2,\ \text{CH}_4,\ \text{CO}\}$:
- Niveles **discretos**
