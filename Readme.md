# Entregable — Señales discretas x(t)

Programa: `señales_discretas.py` (generación, muestreo, Nyquist-Shannon, aliasing, FFT y lectura con el cursor).
Comprobación automática: `test_senales_discretas.py` (`python test_senales_discretas.py`).

Al ejecutar el programa se abre **una ventana con dos figuras de 4×2**, una fila por
señal: a la izquierda la curva continua con las muestras superpuestas, a la derecha el
espectro. La primera figura muestrea a 100 Hz (cumple Nyquist-Shannon) y la segunda a
25 Hz (aliasing intencional). **Al mover el cursor por cualquier panel aparece una caja
negra con los datos del punto más cercano:** en los paneles de tiempo, `t`, `x(t)`, la
amplitud de la señal y las frecuencias observadas; en los de frecuencia, `f`, `|X(f)|` y
las frecuencias observadas. Antes de abrirse las ventanas, el programa imprime por consola
la relación `fs` / `f_max` de cada señal y a qué frecuencia se pliega cada componente.

Las señales se definen con el modelo

```text
x(t) = A · sin(2πFt + φ)
```

- `A` = amplitud (magnitud máxima de la señal).
- `F` = frecuencia en Hz (ciclos por segundo).
- `t` = tiempo en segundos.
- `φ` = ángulo de fase (en qué punto del ciclo está la señal al empezar).

La observación dura `T = 2 s`. Los paneles de tiempo muestran el primer segundo
(`T_PANEL = 1 s`) para que las muestras se vean separadas; el cálculo y la FFT usan los
2 s completos, es decir 200 muestras a 100 Hz y 50 a 25 Hz.

---

## 2.1 Identificación de las señales

Las cuatro señales se generaron con los valores iniciales siguientes. La duración de la
observación es común a las cuatro (`T = 2 s`); la amplitud, la frecuencia y la fase son
los parámetros de cada una y están identificados con un nombre en el programa, de modo
que se modifican en la propia definición.

### Señal 1

- **Nombre:** S1: Senoidal simple
- **Amplitud:** A = 1.5
- **Frecuencia:** F = 10 Hz (φ = 0°)
- **Tiempo:** 0 a 2 s → 20 ciclos completos
- **Ángulo de fase:** φ = 0° (arranca en el cero, subiendo)
- **Color:** azul
- **Modelo matemático:** x₁(t) = 1.5 · sin(2π·10·t) = 1.5 · sin(20πt)

### Señal 2

- **Nombre:** S2: Multitono
- **Amplitud:** A₁ = 1.0 (a 5 Hz) y A₂ = 0.8 (a 20 Hz)
- **Frecuencia:** F₁ = 5 Hz y F₂ = 20 Hz
- **Tiempo:** 0 a 2 s → 10 ciclos de 5 Hz y 40 de 20 Hz
- **Ángulo de fase:** φ₁ = 0° y φ₂ = 90° (0.8 · sin(2π·20t + 90°) = 0.8 · cos(2π·20t))
- **Color:** verde
- **Modelo matemático:** x₂(t) = 1.0 · sin(2π·5·t) + 0.8 · cos(2π·20·t)

### Señal 3

- **Nombre:** S3: Modulada en amplitud (AM)
- **Amplitud:** portadora A = 1.0 con modulante de amplitud 0.5
- **Frecuencia:** portadora F_c = 15 Hz, modulante F_m = 2 Hz
- **Tiempo:** 0 a 2 s → 30 ciclos de portadora
- **Ángulo de fase:** φ_c = 0° y φ_m = 90° (la modulante es un coseno)
- **Color:** naranja
- **Modelo matemático:** x₃(t) = (1 + 0.5 · cos(2π·2·t)) · sin(2π·15·t)

El factor entre paréntesis es la envolvente: vale 1.5 en los máximos del coseno y 0.5 en
los mínimos, de modo que la portadora de 15 Hz aparece y desaparece cuatro veces en los
2 s. Al multiplicar, el producto genera tres componentes: 15 Hz (la portadora) y dos
bandas laterales en 13 Hz y 17 Hz, cada una de amplitud 0.25. **Por eso la frecuencia
máxima de esta señal no es 15 Hz sino 17 Hz**, aunque 15 Hz sea la portadora.

### Señal 4

- **Nombre:** S4: Decaimiento exponencial
- **Amplitud:** A = 2.0 al inicio (el máximo absoluto de la señal)
- **Frecuencia:** F = 30 Hz
- **Tiempo:** 0 a 2 s → 60 ciclos, con envolvente e^(−1.2t)
- **Ángulo de fase:** φ = 0°
- **Color:** púrpura
- **Modelo matemático:** x₄(t) = 2 · e^(−1.2t) · sin(2π·30·t)

La envolvente decae con una constante de tiempo de 1/1.2 ≈ 0.83 s: la señal vale 2 al
comenzar y solo 0.18 en t = 2 s. Representa un golpe o una vibración que se amortigua.

---

## 2.2 Aplicación del teorema de Nyquist-Shannon

El teorema dice que una señal limitada en banda de frecuencia máxima `f_max` se
reconstruye sin pérdida a partir de sus muestras si

```text
fs >= 2 · f_max
```

`fs/2` es el **límite de Nyquist** de esa medición: por encima de él, dos frecuencias
distintas dan exactamente las mismas muestras y el problema no tiene solución única.

La `f_max` de cada señal es su frecuencia oscilatoria más alta, contando las bandas
laterales de la modulación:

| Señal | f_max | fs mínima (2·f_max) | fs usada | Límite de Nyquist | ¿Cumple? |
|---|---:|---:|---:|---:|---|
| S1 Senoidal | 10 Hz | 20 Hz | 100 Hz | 50 Hz | Sí |
| S2 Multitono | 20 Hz | 40 Hz | 100 Hz | 50 Hz | Sí |
| S3 AM | 17 Hz | 34 Hz | 100 Hz | 50 Hz | Sí |
| S4 Decaimiento | 30 Hz | 60 Hz | 100 Hz | 50 Hz | Sí |

Para el conjunto de las cuatro señales el criterio se aplica a la más rápida:
`f_max = 30 Hz → fs >= 60 Hz`. Los 100 Hz elegidos superan ese mínimo con un margen de
1.67 veces, de modo que ninguna de las cuatro señales se pliega.

Un matiz sobre S4: la envolvente exponencial hace que la señal **no sea estrictamente
limitada en banda** —su espectro es una lorentziana que nunca llega a cero—, así que en
un caso real no basta con mirar 30 Hz. Habría que medir el ancho de banda efectivo y
colocar un filtro antialias antes del conversor analógico-digital. Aquí 30 Hz se usa como
frecuencia oscilatoria característica para ilustrar el teorema.

**El criterio se aplica a la frecuencia más alta de la señal completa, no a su
frecuencia fundamental ni a la portadora.** S3 es el ejemplo: con 15 Hz de portadora y
100 Hz de muestreo todo parecería bien, pero la componente real de 17 Hz es la que
justifica el mínimo de 34 Hz.

---

## 2.3 Análisis del aliasing

El caso de aliasing se genera bajando `fs` a 25 Hz, cuyo límite de Nyquist es 12.5 Hz.
No cumple para tres de las cuatro señales. La frecuencia a la que se pliega cada
componente se calcula con

```text
f_alias = |f − fs · round(f / fs)|
```

que es la función `f_alias` del programa, y coincide con la posición de los picos del
espectro medido:

| Señal | Componente original | f_alias a 25 Hz | Cómo se ve |
|---|---:|---:|---|
| S1 | 10 Hz | 10 Hz | **No hay aliasing**: está bajo el límite de 12.5 Hz. |
| S2 | 5 Hz | 5 Hz | Correcta, no se pliega. |
| S2 | 20 Hz | **5 Hz** | Cae justo encima de la componente real de 5 Hz: las dos se confunden. |
| S3 | 13 Hz | 12 Hz | Se pliega contra el límite. |
| S3 | 15 Hz | 10 Hz | La portadora aparece como una vibración más lenta. |
| S3 | 17 Hz | 8 Hz | Se pliega al otro lado y crea un pico donde no hay nada. |
| S4 | 30 Hz | **5 Hz** | La vibración más rápida se lee como la más lenta del conjunto. |

**Qué se pierde y qué no.** El fallo de Nyquist-Shannon no se ve en los datos: el número
de puntos es el correcto y la frecuencia de muestreo también. Lo que se pierde es la
**capacidad de distinguir**. Dos señales de 5 Hz y 20 Hz, muestreadas a
25 Hz, producen exactamente la misma secuencia de números, y el programa no puede saber
cuál de las dos era la real. En S2 el defecto se ve en el dominio del tiempo: la
componente de 20 Hz se pliega sobre la de 5 Hz y la señal observada es una oscilación de
5 Hz con la amplitud y la fase alteradas. En el tiempo de S4 lo grave es peor: una
oscilación de 30 Hz se ve como una de 5 Hz, seis veces más lenta. Un profesional que
mira el panel de tiempo de la figura de 25 Hz cree que el motor o el sistema es
mucho más lento de lo que es, y podría diagnosticar una pieza que gira a 30 Hz como si
estuviera casi parada.

**Por qué no se puede arreglar después.** La FFT solo ve los mismos números, así que no
puede recuperar la información: si dos frecuencias originales producen las mismas
muestras, el espectro discreto no permite decidir cuál era la original. La única
prevención es muestrear **más rápido que el doble de la frecuencia más alta** o un filtro
antialias antes de digitalizar.

---

## 2.4 Análisis mediante FFT

`muestrear` calcula la FFT de las muestras con `numpy.fft`, toma la mitad positiva del
espectro y escala la magnitud con `2/N` para obtener un espectro unilateral. La resolución
es `fs/N = 1/T = 0.5 Hz`, de modo que los picos caen exactamente sobre bins de 0.5 Hz.

**Con 100 Hz (cumple Nyquist).** Cada pico aparece en la frecuencia real, con la amplitud
esperada:

| Señal | Picos medidos a 100 Hz | Amplitudes |
|---|---|---|
| S1 | 10 Hz | 1.5 |
| S2 | 5 Hz y 20 Hz | 1.0 y 0.8 |
| S3 | 13 Hz, 15 Hz y 17 Hz | 0.25, 1.0 y 0.25 |
| S4 | 30 Hz (con bins vecinos por la envolvente) | 0.76 |

S3 confirma con números lo que dice la envolvente: la portadora de 15 Hz domina, y a sus
lados aparecen dos picos iguales y pequeños en 13 y 17 Hz, exactamente la mitad de la
amplitud de modulación.

Dos detalles honestos de esta tabla:

- El pico de S4 mide 0.76 y no 2.0. La escala `2/N` reproduce la amplitud de una señal
  **periódica**, y S4 no lo es: su energía está concentrada en el primer segundo y se
  reparte en los bins vecinos, como se ve en la cola del espectro. Es una limitación de
  aplicar una DFT a una señal transitoria, no un error del muestreo.
- La resolución de 0.5 Hz es gruesa. Con 2 s de observación, 30 Hz y 5 Hz caen en bins
  diferentes, pero dos componentes separadas por menos de 1 Hz no se distinguirían. Para
  mejorar la resolución hace falta más tiempo de observación, no más frecuencia de
  muestreo.

**Con 25 Hz (aliasing).** El espectro es el mismo, con las frecuencias plegadas: el pico
de S4 está en 5 Hz con 0.76 de magnitud, el de la portadora de S3 en 10 Hz y sus bandas
laterales en 8 y 12 Hz, y en S2 solo queda un pico en 5 Hz de 1.28 — la suma en
cuadratura de 1.0 y 0.8, porque las dos componentes se pliegan una sobre la otra con 90°
de diferencia de fase. La forma del espectro es engañosa: **parece un sistema lento y
limpio, con una sola vibración dominante**, cuando en realidad hay tres componentes rápidas
superpuestas. El diagnóstico sería erróneo en las dos direcciones: se atribuiría a una
frecuencia que no existe y se pasarían por alto las que sí.

---

## 2.5 Aplicación profesional

### Campo seleccionado

**Campo:** Ingeniería de mantenimiento — monitoreo de vibraciones de un motor
(automotriz o industrial) mediante un acelerómetro.

### Aplicación

Un acelerómetro —un sensor piezoeléctrico o uno MEMS que devuelve una aceleración— se
pega al bloque del motor o a una carcasa y su salida es una tensión proporcional a la
vibración. Esa tensión es analógica y continua: un tren de oscilaciones de amplitud
variable alrededor de cero. Para guardarla en un ordenador, mostrarla en una pantalla o
mandarla por una red, hay que digitalizarla, y digitalizar significa muestrear.

Las cuatro señales de este programa son los cuatro patrones que aparecen en la práctica:

| Señal del programa | Patrón real en un motor | Qué lo produce |
|---|---|---|
| S1: senoidal simple | La vibración de primer orden de una pieza que gira a ritmo constante. | El desequilibrio o el desalineado de un eje. |
| S2: multitono | Dos o tres mecanismos a la vez. | Dos mecanismos con rpm distintas, o un motor y su accesorio. |
| S3: modulada en AM | La misma vibración cuya amplitud respira. | Una leva, un engranaje con dientes o un pistón que carga y descarga cada vuelta. |
| S4: decaimiento exponencial | Un golpe puntual que se amortigua. | El impacto de un balancín o de un rodamiento al partir. |

La decisión de diseño es la tasa de muestreo, y se toma a partir de la frecuencia más
alta que hay que vigilar, no por costumbre. Si el régimen del motor llega a 3000 rpm, el
primer orden está a 50 Hz y sus armónicos a 100, 150 o 200 Hz: un equipo que muestrea a
1 kHz los ve sin problema, y uno que muestrea a 200 Hz solo puede representar hasta
100 Hz. Los equipos de adquisición industriales vienen con un filtro
antialias conmutado precisamente porque el error de muestreo no avisa: la señal guardada
parece válida, solo que movida de sitio.

### Explicación para un profesional

Si a usted le entregaran únicamente el dibujo de las curvas de la segunda figura —donde
las muestras se toman 25 veces por segundo— y se le pidiera un diagnóstico, casi todas
las conclusiones serían incorrectas. **La forma de la curva no dice nada por sí sola.** Lo
que dice es hasta dónde puede mirar el instrumento.

Con 100 mediciones por segundo, el equipo registra una muestra cada 0.01 s y su límite
de información está en 50 Hz. Todo lo que ocurra más rápido que eso queda registrado sin
ambigüedad: la FFT los devuelve en su sitio y usted puede leer "esta máquina vibra a
15 Hz con bandas laterales en 13 y 17", que es información de diagnóstico directa.

Con 25 mediciones por segundo, el mismo equipo **no ha fallado**: ha recogido el mismo
número de puntos, con la misma limpieza, y ha mentido igual. Su límite baja a 12.5 Hz.
Una vibración de 20 Hz entra y sale como una de 5 Hz; una de 30 Hz, también. Para usted,
que no ve el muestreo, el motor parece más lento de lo que es. La consecuencia práctica es
concreta: si el motor ha desarrollado una resonancia a 30 Hz, el registro con 25 Hz la
esconde detrás de un pico de 5 Hz que no corresponde con nada en la máquina, el técnico
busca una avería donde no la hay y la resonancia sigue creciendo hasta dañar el
equipo. Y el caso inverso también es real: un pico de 5 Hz espurio puede hacer creer que
hay un problema de alineación que en realidad no existe.

La lectura correcta de estas gráficas lleva siempre tres preguntas en el mismo orden.
**¿Qué estoy viendo?** Las componentes de la señal, con su frecuencia y su amplitud. **¿Hasta
dónde puede mirar el instrumento?** La mitad de la frecuencia de muestreo, el límite de
Nyquist. **¿He mirado más allá?** Si la frecuencia sospechada está por encima de ese
límite, el número que tengo delante no describe la máquina: describe el instrumento. En
la primera figura la respuesta es "todo visible, nada fuera de alcance". En la segunda, la
respuesta es "hay cuatro señales y solo puedo ver la mitad del espectro", y por eso los
picos están desplazados.

En un sistema real esto se traduce en una regla práctica: se mide a una frecuencia de
muestreo **claramente superior** al doble de la frecuencia más alta que interesa vigilar —
no al doble de la que resulta más cómoda— y se coloca un filtro antialias antes del
conversor. Sobre esa base, la FFT es lo que convierte una curva en un diagnóstico, porque
da la firma de cada mecanismo, y el muestreo es lo que decide si esa firma es real o un
reflejo del instrumento.

---

## Conclusiones

El programa construye cuatro señales, las muestrea a dos frecuencias de muestreo y
compara sus espectros. A 100 Hz, por encima del mínimo nominal de 60 Hz que impone la
componente más rápida del conjunto, cada pico de la FFT aparece en su frecuencia real y
con la amplitud esperada. A 25 Hz, con un límite de Nyquist de 12.5 Hz, las componentes
de 20 y 30 Hz se pliegan a 5 Hz y las bandas laterales de la señal modulada caen en 8, 10
y 12 Hz: el espectro resultante describe un sistema más lento del que es.

Las tres conclusiones del trabajo:

1. **El muestreo decide qué es visible.** La misma señal medida con más frecuencia de
   muestreo da información distinta, y la diferencia no está en la señal sino en la
   decisión de adquisición.
2. **El aliasing no se detecta en los datos, se descarta en el diseño.** Las muestras
   son válidas y limpias; lo que falta es la capacidad de distinguir frecuencias, y eso
   solo se sabe mirando la relación entre `fs` y `f_max` antes de medir.
3. **Nyquist-Shannon se aplica a la componente más alta, no a la principal.** En S3 la
   portadora está en 15 Hz, pero la banda lateral de 17 Hz es la que fija el mínimo de
   muestreo, y en S4 la envolvente obliga a ir más allá de los 30 Hz nominales.

Queda pendiente para una implementación industrial: medir el ancho de banda efectivo del
sensor, conmutar filtros antialias, considerar el ruido y la resolución del conversor, y
repetir el análisis con ventanas más largas cuando la resolución de 0.5 Hz no alcance
para separar componentes cercanos.
