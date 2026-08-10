---
title: "Cinco señales publicadas para actualizar metanálisis, aplicadas y evaluadas: un estudio exploratorio sobre cuándo pueden usarse siquiera"
author: |
  Javier Núñez\
  \small Investigador independiente · ORCID [0009-0003-9770-4986](https://orcid.org/0009-0003-9770-4986) · jnverbel@gmail.com
date: 2026-08-10
lang: es
keywords: [actualización de revisiones sistemáticas, metanálisis, síntesis de evidencia, metodología de la investigación, reproducibilidad]
---

# Resumen

**Antecedentes.** Entre 1999 y 2007 se publicaron cinco métodos estadísticos
para decidir si una revisión sistemática ha quedado desactualizada. Se han
ejecutado uno junto a otro exactamente una vez, sobre 80 revisiones Cochrane,
donde dos no marcaron nada y los tres que sí discriminaron concordaron con un
Kappa de 0,14. No parece existir ninguna implementación reutilizable de
ninguno de ellos, de modo que la comparación no se ha repetido.

**Métodos.** Implementamos los cinco en un paquete de R de código abierto y los
aplicamos a evidencia histórica en dos brazos. El brazo A barre los 17
metanálisis de la colección `metadat` que traen datos por estudio y años de
publicación, y ejecuta los cinco detectores contra tres definiciones
operativas de que la estimación combinada se ha movido. El brazo B reúne 6.686
pares de versiones consecutivas de las 4.132 revisiones Cochrane con más de
una versión, cosechadas de Europe PMC, y les adjunta el desenlace editorial
registrado en las conclusiones de los autores de cada versión, lo que da
soporte a los tres detectores que solo necesitan estimaciones combinadas.

**Resultados.** En el brazo A, 168 de 185 cortes (91 %) tenían un metanálisis
previo ya significativo, de modo que `barrowman` y `simulation` —que exigen un
previo no significativo— solo pudieron interrogarse en 4 y 5 de las 17
revisiones. Esto era invisible para la comparación publicada porque su cohorte
se seleccionó por no ser significativa. El criterio de efecto del método
Ottawa es un cociente de reducciones relativas del riesgo cuyo denominador
tiende a cero a medida que el efecto previo se acerca al nulo, así que dispara
en el 64 % de las muestras que no contienen cambio alguno; su especificidad
cae a 0,14 en una revisión nula donde los otros cuatro se mantienen en 1,00. La
mitad de estabilidad del método de suficiencia se define como la pendiente de
un ajuste por mínimos cuadrados sobre la serie acumulada de efectos; esa
pendiente carece de distribución nula válida y disparó en 209 de 300 muestras
de evidencia inmóvil. En el brazo B, 560 pares llevan a la vez una estimación
combinada comparable y las conclusiones de los autores en ambos extremos; un
cribado automático separa los cambios probables de conclusión con una razón de
11,3 entre estratos.

**Conclusiones.** Los cinco se reparten en dos clases de problema, y la
distinción importa. Dos sencillamente no están definidos una vez que el
metanálisis previo es significativo, que es el caso en el 91 % de los cortes
aquí: eso es una restricción de dominio y no un defecto, pero significa que
normalmente ni siquiera pueden interrogarse, y la única comparación publicada
no pudo verlo porque su cohorte estaba seleccionada por la única condición bajo
la cual sí pueden. Los otros dos fallos sí son defectos: un criterio es
inestable por construcción justo en las revisiones a las que su método apunta,
y un estadístico carece de distribución nula válida. Son hallazgos
exploratorios sobre 17 revisiones sin conjunto retenido, contrastados contra
objetivos operativos que observan el movimiento de la estimación combinada y no
lo que ningún equipo de revisión hizo. El programa, los barridos y el corpus son
reproducibles a partir de fuentes públicas.

---

# 1. Introducción

Un metanálisis es una fotografía, no un hecho permanente. Alrededor de una
cuarta parte de las revisiones sistemáticas quedan desactualizadas a los dos
años de publicadas, y la mitad a los cinco y medio [@shojania2007]. Decidir
*cuándo* una revisión concreta necesita actualizarse sigue siendo en buena
medida una cuestión de juicio, y la lista de comprobación de consenso para ese
juicio [@garner2016] trata la estimación combinada como una entrada entre
varias.

Se han propuesto cinco métodos estadísticos para informar la parte
cuantitativa de esa decisión, publicados entre 1999 y 2007: el metanálisis
acumulativo recursivo [@ioannidis1999], los indicadores de suficiencia y
estabilidad [@mullen2001], el criterio de razón de participantes de Barrowman
[@barrowman2003], el método Ottawa [@shojania2007] y la simulación de potencia
prospectiva [@sutton2007].

Se han ejecutado uno junto a otro exactamente una vez. Pattanittum y
colaboradores [@pattanittum2012] aplicaron los cinco a 80 revisiones Cochrane
y reportaron que dos no marcaron absolutamente nada, mientras que los tres que
sí discriminaron —Ottawa marcando 34 revisiones, el metanálisis acumulativo
recursivo 7 y Barrowman 7— concordaron con un Kappa de 0,14, indistinguible
del azar. Solo una revisión fue marcada por los tres a la vez.

Ese resultado lleva más de una década sin examinarse, y la razón es prosaica.
Cada método se publicó como una descripción en un artículo, y no pudimos
encontrar ninguna implementación reutilizable de ninguno. Una búsqueda en los
metadatos de los 24.734 paquetes de CRAN, hecha el 2026-08-09, no devuelve
ninguna implementación del método Ottawa, de Barrowman, ni del metanálisis
acumulativo recursivo como diagnóstico de actualización; cada coincidencia
sobre esos nombres es un falso positivo, y todas quedan adjudicadas una por
una en `inst/cran-search/` del paquete. La afirmación es *no encontramos
ninguna*, no *no existe ninguna*: la búsqueda lee metadatos de paquetes y no
código fuente, y el propio `metafor` no coincide ni con «cumulative
meta-analysis» ni con «fail-safe» en sus propios metadatos, aun exportando
`cumul()` y `fsn()`.

La consecuencia es que un equipo que quiera saber en qué señal de
actualización confiar tiene primero que reimplementar cinco métodos a partir
de prosa, y cada reimplementación es una oportunidad nueva de diferir de todas
las demás.

Este artículo reporta qué ocurre cuando los cinco se implementan una vez, en
abierto, y se ejecutan contra historia real. Los hallazgos son en su mayoría
negativos, y el más útil trata de **cuándo pueden usarse siquiera los métodos**
antes que de si aciertan.

# 2. Qué establece este estudio y qué no

Enunciamos las cotas antes que los resultados, porque califican cada cifra que
sigue.

**Los objetivos de evaluación no son desenlaces.** El brazo A contrasta cada
detector con tres definiciones operativas de que la estimación combinada se
movió entre un punto de corte y un objetivo posterior. Ninguna de las tres
observa lo que un equipo de revisión hizo, si cambió alguna recomendación, ni
si alguien resultó perjudicado por actuar sobre la estimación antigua. Una
sensibilidad de 0,32 significa por tanto «coincidió con un criterio declarado
sobre la estimación el 32 % de las veces», no «acertó el 32 % de las veces».
Dos de los tres comparten un numerador idéntico y difieren solo en qué error
estándar los divide, así que son una distancia en dos escalas y no
comprobaciones independientes; el acuerdo entre ellos es lo esperable y el
desacuerdo es el caso informativo.

**El brazo A son 17 revisiones y no hay nada retenido.** Las 17 son lo que
sobrevive de las 110 tablas de `metadat`, y 54 de las 93 exclusiones caen por
una sola causa: el conjunto de datos no registra el año de publicación por
estudio. Los que sí lo registran tienden hacia los clásicos bien curados, de
modo que las 17 son una muestra de conveniencia cuya selección puede
correlacionar con lo que se mide, en una dirección que desconocemos. Toda
cifra del brazo A se genera y se enuncia sobre esas mismas 17.

**Dos de los cinco detectores no son reproducciones literales.** La mitad de
estabilidad del método de suficiencia se sustituye por un estadístico de punto
de cambio, por las razones del §5.3, y la sustitución va en el nombre de la
función. El detector de simulación simula efectos y no participantes, porque
no hay datos a nivel de participante disponibles. Ambas desviaciones están
medidas; ninguna está oculta.

**El brazo B no alcanza los datos por estudio.** El texto completo de Cochrane
está tras un muro de pago, y las tablas de análisis están ausentes del texto
completo disponible en abierto en Europe PMC (comprobado en cinco revisiones:
ninguna tabla de datos en ninguna). El grano más fino disponible a escala es
por tanto la estimación combinada, lo que da soporte a `rcma`, `ottawa` y
`barrowman`, y no a `sufficiency_changepoint` ni a `simulation`. Eso es una
propiedad del acceso a los datos, no del diseño.

# 3. Métodos

## 3.1 Implementación

Los cinco detectores están implementados en `staleness`, un paquete de R con
licencia MIT construido enteramente sobre `metafor` [@viechtbauer2010] para el
ajuste de modelos. Cada detector es una función pura con una firma común y
devuelve un veredicto, una señal y un registro de detalle. Un detector que no
puede responder devuelve `not_applicable` y no `current`, de modo que
abstenerse nunca se contabiliza como acierto.

La retroevaluación (*backtesting*) reajusta el metanálisis en cada punto de
corte a partir de los estudios disponibles en esa fecha, así que ningún
detector ve evidencia que todavía no existía. Los cortes demasiado próximos al
final de la serie como para tener un futuro observable se marcan como
censurados y se excluyen de las métricas, en lugar de contrastarse contra un
objetivo que no puede conocerse.

Cuando un procedimiento publicado resulta ambiguo, la lectura adoptada se
declara y se justifica a partir de ejemplos trabajados en la literatura que lo
aplica, en vez de elegirse por conveniencia. El §5.2 recoge el caso donde esto
más importó.

## 3.2 Brazo A: barrido histórico con datos por estudio

Se cribó cada tabla de `metadat` exigiendo año de publicación por estudio, una
medida de efecto de dos grupos construible con `escalc()`, al menos ocho
estudios, un solo efecto por par estudio-año, y al menos tres cortes no
censurados. Sobrevivieron diecisiete. Cada una se retroevaluó con un único
juego de parámetros (`horizon = 3`, `window = 5`, `min_k = 3`) para que la
comparación fuera equivalente, y los cinco detectores se ejecutaron en cada
corte anual.

## 3.3 Brazo B: cadenas de versiones con un desenlace editorial

Cochrane asigna a cada revisión un identificador estable y a cada versión un
sufijo dentro del DOI: `10.1002/14651858.CD005343.pub7` es la revisión
`CD005343`, versión 7. Agrupar por identificador y ordenar por sufijo
reconstruye por tanto el historial completo de versiones de una revisión solo
a partir de metadatos, sin ninguna consulta adicional; Crossref expone además,
de forma independiente, las relaciones `update-to` y `updated-by` con fechas.

Cosechamos 17.247 versiones de revisiones de la API REST abierta de Europe
PMC, que cubren 9.862 revisiones distintas. De ellas, 4.132 tienen más de una
versión y aportan por tanto los 6.686 pares consecutivos que representan una
actualización; el resto tiene una sola versión y no puede aportar ninguno. La
sección de conclusiones de los autores está presente en el resumen gratuito de
cada versión, lo que da la posición editorialmente registrada en ambos
extremos de un par.

Comparar estimaciones combinadas a lo largo de un par exige que ambas citen el
mismo desenlace. Exigimos la misma medida de efecto y una similitud textual de
al menos 0,80 entre las palabras que preceden a cada estimación; 746 pares lo
cumplen, de los cuales 560 llevan además conclusiones en ambos extremos.

Nada del brazo B toca cochranelibrary.com. Toda la entrada la sirven Europe
PMC o Crossref bajo sus términos normales, que es lo que permite depositar el
corpus junto a este artículo en vez de describirlo y no poder compartirlo.

## 3.4 El brazo B está retenido para los detectores, y no para el cribado

La cronología importa y es verificable, no una afirmación retrospectiva. El
comportamiento de los detectores quedó congelado en el commit `f59bb41`
(2026-08-09 19:30 −0500); el corpus del brazo B se introdujo por primera vez en
el commit `c7e9505` (2026-08-10 07:56 −0500). Desde ese segundo commit,
`git diff` reporta **cero líneas de código cambiadas** en cualquier archivo
fuente de detector: las únicas ediciones son de documentación. Ningún
parámetro, ningún umbral y ninguna interpretación de un criterio publicado
ambiguo se decidió mirando el brazo B, porque ninguno pudo haberlo sido.

El brazo B es por tanto un conjunto genuinamente retenido **para los
detectores**. No lo es para el cribado de conclusiones, que se desarrolló y se
calibró dentro del corpus y está validado solo en la medida en que su muestra
estratificada ha sido codificada. Son dos afirmaciones distintas y no las
mezclamos:

| | estatus |
|---|---|
| Brazo A | desarrollo y exploratorio |
| Brazo B, evaluación de detectores | retenido |
| Brazo B, cribado de conclusiones | desarrollo, parcialmente validado |

La consecuencia es que una evaluación confirmatoria de los detectores no exige
reunir evidencia nueva. Exige especificar el análisis antes de ejecutarlo: qué
detectores son evaluables, el estimando de cada uno, las exclusiones, el
tratamiento de `not_applicable`, y qué cuenta como éxito — escrito contra un
commit congelado y ejecutado una sola vez. Esa evaluación no se reporta aquí, y
todo lo que corramos después irá etiquetado como post hoc.

# 4. Reproducibilidad

Tres barridos regeneran toda afirmación cuantitativa de más abajo a partir de
fuentes públicas: `inst/applicability/` produce el brazo A,
`inst/calibration/` regenera cada figura de calibración desde sus semillas
originales **y reconstruye los estadísticos que fueron reemplazados**, de modo
que el antes y el después puede volver a derivarse en lugar de aceptarse de
palabra, y `corpus/` reconstruye el brazo B en unos diez minutos.

El paquete lleva 1.176 pruebas y el corpus 35, con la prueba de mutación como
estándar: para cada guarda, se modifica el código y se exige que la prueba se
ponga en rojo. Se reporta porque detectó errores que la lectura no detectó,
incluidos cuatro en la tubería del corpus descritos en el §5.5.

# 5. Resultados

## 5.1 Dos de los cinco detectores son estructuralmente incapaces de responder

A lo largo de las 17 revisiones del brazo A y sus 185 cortes anuales, **168
cortes (91 %) tenían un metanálisis previo ya significativo**. En 11 de las 17
revisiones, todos los cortes lo tenían.

Tanto `barrowman` como `simulation` exigen un previo no significativo: el
primero pregunta cuántos participantes harían falta para alcanzar
significación, y el segundo simula la potencia del siguiente lote para
lograrla. Ninguna de las dos preguntas está definida una vez que el previo ya
es significativo. Por eso solo pueden interrogarse en **4 y 5 de las 17
revisiones**, respectivamente.

| Detector | ¿Publicado? | Puede responder | Sensibilidad media | Revisiones con cero |
|---|---|---|---|---|
| `ottawa` | tal cual | 17 de 17 | 0,320 | 8 |
| `rcma` | tal cual | 17 de 17 | 0,309 | 9 |
| `barrowman` | tal cual | **4 de 17** | 0,167 | 2 de 3 |
| `simulation` | nivel de efecto | **5 de 17** | 0,167 | 2 de 3 |
| `sufficiency_changepoint` | **sustituto** | 17 de 17 | 0,041 | 11 |

La segunda columna no es decoración. `sufficiency_changepoint` ejecuta un
estadístico que su fuente nunca describió, así que **0,041 es una propiedad de
nuestro sustituto y no aporta ninguna información sobre el método de
suficiencia publicado** —cuyo propio estadístico de estabilidad se muestra en
el §5.3 que carece de distribución nula válida y por tanto no tiene
sensibilidad significativa que reportar—. `simulation` simula efectos y no
participantes. Ninguna de las dos filas puede leerse como evidencia sobre el
procedimiento al que su nombre alude, y ninguna conclusión de más abajo lo
hace.

Esto no es un hallazgo sobre que esos métodos estén equivocados. Es un
hallazgo sobre cuándo pueden usarse, y era **invisible para la comparación
publicada porque sus 80 revisiones se seleccionaron por tener un resultado
combinado no significativo** — la única cohorte en la que esos dos detectores
siempre pueden interrogarse.

La cifra no depende de la configuración de la retroevaluación: ejecutado con
`horizon = 6`, el barrido devuelve las mismas 17 revisiones, los mismos 168 de
185 cortes y la misma cobertura. Eso es esperable antes que tranquilizador
—que un previo ya fuera significativo es un hecho sobre la evidencia en un
corte— pero significa que el hallazgo es una propiedad de las revisiones y no
de esta parametrización.

## 5.2 El criterio de efecto de Ottawa es inestable justo donde el método debe usarse

El método Ottawa enuncia sus señales cuantitativas como «cambios en la
significación estadística o cambios relativos en la magnitud del efecto de al
menos el 50 %» [@shojania2007]. La formulación no dice *de qué*, y la elección
importa: las dos lecturas dan detectores materialmente distintos.

Dos aplicaciones independientes la resuelven del mismo modo, sobre datos que
pueden comprobarse. Pattanittum et al. [@pattanittum2012], tabla 1, calcula el
cambio sobre las *reducciones* relativas del riesgo, `(1 − RR_new) /
(1 − RR_prev)`; de las diez revisiones con el mayor indicador Ottawa del
apéndice de ese estudio, las diez disparan bajo esa lectura y ninguna dispara
bajo el cociente de las razones de riesgo. Mickenautsch y Yengopal
[@mickenautsch2013] trabajan dos ejemplos más —RR de 2,10 a 1,51 y RR de 2,61
a 1,66— describiendo ambos, con sus propias palabras, como «un cambio en el
tamaño relativo del efecto de más del 50 %». Ninguno lo cumple sobre el
cociente de los efectos (0,719 y 0,636) ni como cambio porcentual de la
estimación (28 % y 36 %); una cuarta lectura candidata, el cambio expresado
sobre la estimación *nueva*, da 39 % y 57 %, y por tanto supera el umbral en un
ejemplo y no en el otro. Solo el cociente de reducciones del riesgo encaja en
los cuatro ejemplos trabajados de los dos artículos.

La corrección arrastra un hallazgo. El denominador `1 − RR_prev` tiende a cero
a medida que el efecto previo se acerca a la ausencia de efecto, de modo que el
criterio es inestable precisamente sobre los metanálisis nulos a los que el
método apunta. Sobre evidencia simulada que **no contiene cambio alguno**, la
señal de efecto dispara en el **64 % de las muestras bajo un efecto nulo** y en
el **0 %** una vez que el efecto es real y preciso. Confirmado fuera de
simulación en `metadat::dat.laopaiboon2015`, una revisión nula donde la
especificidad de este detector cae a **0,14** mientras los otros cuatro se
mantienen en 1,00.

Esto no es un artefacto de implementación; se sigue del criterio tal como está
escrito. Explica además la comparación publicada: Ottawa marcó 34 de 80
revisiones donde el metanálisis acumulativo recursivo y Barrowman marcaron 7
cada uno, sobre una cohorte de metanálisis nulos por criterio de inclusión.

No lo corregimos. Corregirlo supondría implementar un método distinto, y el
propósito del ejercicio es averiguar cómo se comportan los publicados. Lo que
sí corregimos es el silencio: el veredicto ahora reporta el denominador por el
que dividió, y señala cuándo ese denominador está cerca de cero.

## 5.3 El estadístico de estabilidad de la suficiencia carece de distribución nula válida

La mitad de estabilidad del método de suficiencia se describe como la
«pendiente absoluta de la regresión lineal ajustada sobre los efectos de
tratamiento acumulados frente al incremento de información»
[@pattanittum2012]. Tomada al pie de la letra, la regla es degenerada —sobre
datos continuos esa pendiente nunca es exactamente cero— así que alguna regla
de significación tiene que ocupar su lugar.

Se probaron y se midieron tres implementaciones:

1. **La prueba *t* de mínimos cuadrados.** Una media acumulada está
   autocorrelacionada casi perfectamente por construcción y converge sobre el
   efecto combinado por la ley de los grandes números, de modo que la prueba
   detecta *convergencia* y reporta *inestabilidad*. Sobre 300 muestras de
   evidencia genuinamente inmóvil devolvió `out_of_date` **209 veces**, donde
   `rcma` y `ottawa` no lo devolvieron ninguna.

2. **Permutar el orden de los estudios, conservando la pendiente.** Esto lleva
   la tasa de falsas alarmas a 16 de 300, nominal en promedio, pero la
   pendiente de una serie acumulada está dominada por sus primeros puntos, así
   que un cambio *tardío* no puede superar la nula de permutación: la potencia
   frente a diez estudios nuevos con RR de 0,30 fue de **1 entre 200**, y caía
   a cero a medida que el desplazamiento crecía. Peor aún, la nula de
   permutación es en sí misma inválida cuando las varianzas de los estudios
   cambian a lo largo del tiempo calendario: con ensayos pequeños al principio,
   grandes después y ninguna deriva en absoluto, disparó en el **28 %** de las
   muestras, y sobre un esquema de 20 ensayos pequeños seguidos de 10 grandes,
   en el **42 %**.

3. **Un estadístico de punto de cambio** —la mayor diferencia estandarizada
   entre los estudios anteriores y posteriores a cualquier división—
   contrastado contra la misma nula de permutación por orden. Falsas alarmas
   15 de 300 (5,0 %); potencia de 200 de 200 frente a diez estudios nuevos en
   cada uno de cuatro tamaños de efecto; y la tasa de falsas alarmas sin
   deriva bajo heterocedasticidad baja a 5,3 %, 6,7 % y 6,7 % en los esquemas
   que arriba produjeron 28 % y 42 %.

El tercero es lo que ejecuta el paquete, y el detector se llama
`sufficiency_changepoint()` y no `sufficiency()` para que la sustitución sea
visible en el punto de llamada y en toda tabla de resultados. La pendiente
publicada se sigue calculando y devolviendo como diagnóstico; no decide nada.

El sustituto tampoco está uniformemente calibrado, y reportamos dónde falla con
el mismo cuidado con que reportamos dónde falla el original. A lo largo de once
regímenes de varianza la tasa de falsas alarmas va del 2,4 % (conservadora) al
11,1 %, aproximadamente el doble de lo nominal, con el peor caso en una rampa
de precisión monótona y suave de 650 a 1. El patrón no es la dirección de la
tendencia —crecimiento y decrecimiento son igual de malos— sino su suavidad,
monotonía y rango: las mismas varianzas de `dat.bcg` reordenadas de forma
creciente llevan la tasa de 5,0 % a 9,4 %, que es la demostración más limpia
de que lo que importa es la disposición y no los números.

Un segundo límite importa más para la interpretación, y es el único punto en
que este artículo corrige una cifra que el propio paquete había publicado. El
estadístico es un máximo sobre puntos de división, y la división que aísla un
único estudio está entre ellos, así que un solo estudio discordante que llegue
último puede cargarlo. Medido por `inst/persistence/persistence.R` sobre 400
réplicas por celda, contra una línea base sin atípico alguno:

| k | sin atípico | un estudio a 5 EE, último | un estudio a 5 EE, en medio |
|---|---|---|---|
| 20 | 0,050 | 0,035 | 0,005 |
| 30 | 0,052 | 0,207 | 0,003 |
| 40 | 0,043 | 0,395 | 0,000 |

El efecto es real y crece con `k` —con 40 estudios, un único estudio final
discordante fuerza `unstable` unas ocho veces más a menudo que el azar— pero
está confinado a los estudios que llegan **últimos**. Uno situado en mitad de
la serie no hace prácticamente nada, porque la división que lo aísla deja
bloques grandes a ambos lados.

Una versión anterior de esta afirmación, en la documentación del paquete,
reportaba tasas mucho más altas (0,995 con `k` = 40) y sostenía que un estudio
así fuerza `unstable` «casi siempre a partir de `k` = 25 aproximadamente». No
pudimos reproducir esas cifras. Una reimplementación independiente que coincide
exactamente con el estadístico que envía el paquete en 200 entradas aleatorias,
y que reproduce la fila de control sin atípico, devuelve entre 0,44 y 0,48 con
`k` = 40 bajo todas las lecturas de «cinco errores estándar» que probamos, y se
satura cerca de 0,55 incluso con veinte. La tabla documentada no tenía ningún
guion generador en el repositorio, que es por lo que sobrevivió; las cifras de
arriba sí lo tienen, y la afirmación cualitativa —que en la cola esto es en
buena medida un detector de un solo atípico— sobrevive con ellas, mientras que
la magnitud no.

### Una regla de persistencia arregla la mayor parte

Si un desplazamiento tiene que persistir para contar, la división que aísla un
solo estudio deja de estar disponible. Exigir que ambos lados de una división
tengan al menos `r` estudios da, promediando sobre `k` = 20, 30 y 40:

| r | deriva detectada | un atípico, último | un atípico, en medio | sin cambio |
|---|---|---|---|---|
| 1 (el que se envía) | 0,802 | 0,212 | 0,003 | 0,048 |
| 2 | 0,815 | 0,226 | 0,003 | 0,051 |
| 3 | 0,825 | 0,147 | 0,006 | 0,054 |
| 5 | 0,845 | **0,083** | 0,013 | 0,052 |

Exigir cinco es estrictamente mejor sobre esta evidencia: la potencia frente a
un desplazamiento tardío real sube un poco (de 0,802 a 0,845), el falso
positivo por un solo atípico cae a menos de la mitad (de 0,212 a 0,083), y la
calibración sin cambio queda intacta (de 0,048 a 0,052). Con `k` = 40 la tasa
del atípico baja de 0,395 a 0,090.

No lo hemos adoptado. La medición son cuatro regímenes sobre evidencia simulada
con un solo esquema de varianza, y el paquete envía el estadístico que estos
resultados cuestionan en vez del que favorecen, porque cambiarlo con la fuerza
de un único experimento repetiría el error que esta misma sección documenta.
Queda enunciado como la siguiente prueba a correr, junto con el guion que la
corre.

## 5.4 Un desenlace editorial es alcanzable, para tres de los cinco detectores

El brazo B adjunta a cada par de versiones las conclusiones de los autores tal
como quedaron registradas en ambos extremos. De 6.686 pares consecutivos,
4.530 llevan conclusiones en ambos extremos una vez eliminados los registros
duplicados, y 560 llevan además estimaciones combinadas que pueden compararse.

Se aplicó a todos los pares un cribado que combina similitud textual con el
vocabulario controlado de la propia Cochrane —los niveles de certeza de GRADE,
los enunciados estandarizados de matización, y la aparición o desaparición de
«evidencia insuficiente»—, y se extrajo una muestra de 120 estratificada a lo
largo del rango de puntuación, codificada a ciegas contra la definición de
French et al. [@french2005], para quienes una conclusión cambiada es la que
«altera la sustancia o el significado de una sección o altera la
interpretación», con las reescrituras de estilo explícitamente excluidas.

| Estrato | n | Mayor | Menor | Ninguno | % mayor (IC 95 %) |
|---|---|---|---|---|---|
| Alto | 40 | 34 | 6 | 0 | 85 % (71–93) |
| Medio | 40 | 26 | 14 | 0 | 65 % (50–78) |
| Bajo | 40 | 3 | 20 | 17 | 8 % (3–20) |

El cribado separa con una razón de 11,3, que es lo que lo hace utilizable para
ordenar los pares restantes. La muestra se estratificó en vez de tomarse de la
cabeza del ordenamiento precisamente para que esto pudiera medirse: una
muestra de los N primeros habría mostrado una tasa alta y no habría dicho nada
sobre lo que el cribado deja pasar.

**La tasa resultante del 52 % es seis veces el 9 % que French et al. midieron
sobre 254 revisiones actualizadas, y la reportamos como cifra de trabajo antes
que como hallazgo.** Hay cuatro explicaciones disponibles y ninguna ha sido
comprobada. La mayor es probablemente la selección: el denominador de French
son revisiones actualizadas, mientras que el nuestro son pares que llevan una
estimación cuantitativa comparable en ambos extremos —una revisión que reporta
un efecto combinable dos veces tiene algo que puede moverse, y una que dice
«evidencia insuficiente» de principio a fin no—. La reestructuración de cómo
Cochrane redacta sus conclusiones hacia 2011 es una segunda candidata, el
mayor volumen del reporte moderno por desenlace una tercera, y un codificador
único en lugar de dos investigadores humanos independientes una cuarta.
Distinguir «los corpus difieren» de «el codificador es liberal» requiere o bien
codificación humana de una submuestra, o bien acceso a la codificación
original de French como patrón externo; ninguna de las dos está disponible
todavía.

## 5.5 Lo que las pruebas detectaron y la lectura no

Cuatro defectos de la tubería del corpus fueron silenciosos, y tres de ellos
produjeron salidas bien formadas que habrían entrado en el análisis sin que
nadie lo notara.

Un DOI sin sufijo se leyó como versión 1; dado que Europe PMC guarda varios
registros de índice bajo el mismo DOI pelado, 551 revisiones tenían varios
registros de «versión 1» y 681 pares (9 %) tenían la misma versión en ambos
extremos. Lo que lo destapó fue una **comprobación de rango**: el intervalo
entre versiones salió de −13 años, y un intervalo negativo no puede existir.

142 pares tenían un resumen idéntico byte a byte en ambos extremos. Eso no es
lo mismo que una actualización cuyas conclusiones se sostuvieron —ese caso es
real y frecuente: 662 pares tienen conclusiones idénticas con la sección de
resultados principales visiblemente reescrita— sino el mismo registro llegando
al índice dos veces.

El recuento de pares utilizables para el análisis de detectores se reportó
inicialmente como 1.907 y es 746, porque nada comprobaba que las dos versiones
citaran el mismo desenlace; el 42 % nombra desenlaces claramente distintos y el
12 % ni siquiera comparte medida de efecto. Salió a la luz a partir de un caso
leído a simple vista.

Y 14 pares llevaban un intervalo que no contiene su propia estimación
puntual. Aquí el analizador tenía razón y **los resúmenes publicados están
equivocados**: `MD −2,76; IC 95 % 3,57 a 1,96` con los signos negativos
perdidos, `RD 0,03; IC 95 % −0,01 a −0,07` corriendo al revés. Se rechazan en
vez de repararse, porque adivinar qué dígito está mal inventa datos, y tampoco
se reordenan, porque un signo perdido en ambos extremos no se arregla
intercambiándolos.

Reportamos todo esto porque un corpus ensamblado a partir de metadatos
bibliográficos invita exactamente a esta clase de error, y porque las
comprobaciones que los detectaron fueron sobre todo comprobaciones de rango y
no pruebas sofisticadas.

# 6. Discusión

Los resultados son de dos clases y no deben leerse como una sola. `barrowman` y
`simulation` **no están definidos** una vez que el metanálisis previo es
significativo, que es el caso en el 91 % de los cortes de esta muestra. Eso es
una propiedad de las preguntas que formulan —cuántos participantes harían falta
para alcanzar significación, y qué potencia tiene el siguiente lote para
lograrla— y ninguna de las dos es un defecto. Sí significa que, sobre evidencia
tomada como viene, normalmente no pueden interrogarse.

El criterio de efecto de Ottawa y la pendiente de estabilidad de la suficiencia
son otra cosa, porque ambos fallan allí donde sí aplican: el primero tiene un
denominador que tiende a cero precisamente sobre los metanálisis nulos a los
que su método apunta, y el segundo carece por completo de distribución nula
válida. `rcma` es implementable tal como se describe y se comporta sin
sorpresas.

No son métodos nuevos y esto no es una propuesta para reemplazarlos. La
contribución es que la pregunta ahora puede formularse siquiera, de forma
repetida, por cualquiera, sobre cualquier cuerpo de evidencia — y que cuando se
formula, las respuestas son en su mayoría poco halagadoras.

**La cobertura merece ir junto a la sensibilidad y la especificidad, no en una
nota al pie.** La forma convencional de comparar detectores es una tabla de
contingencia, y una tabla de contingencia no tiene casilla para «este método no
tenía derecho a responder». Nuestros resultados sugieren que la primera
pregunta no es cuál detector acierta más, sino *en qué fracción de situaciones
un detector está siquiera definido* — y que esas dos preguntas pueden tener
respuestas opuestas, porque un método que responde rara vez puede lucir
excelente en las pocas ocasiones en que lo hace. Por eso el paquete devuelve
una fila con `n = 0` y métricas `NA` para un detector que nunca aplicó, en vez
de omitirla: así una ausencia es un valor en la tabla y no algo que el lector
tenga que notar.

El resultado de aplicabilidad merece énfasis porque cambia cómo debe leerse la
única comparación existente. Pattanittum et al. no hallaron que `barrowman` y
`simulation` no marquen nada porque esos métodos sean insensibles; lo hallaron
sobre una cohorte seleccionada por la única condición bajo la cual esos
métodos están definidos. Sobre evidencia tomada como viene, el problema es más
básico: normalmente ni siquiera pueden interrogarse.

La ruta hacia una afirmación más fuerte es visible y está en parte construida.
El brazo B adjunta un desenlace que un humano registró, que es justo lo que
los objetivos operativos del brazo A no son, y alcanza miles de pares en vez de
diecisiete. Lo que no alcanza, porque las tablas de análisis están tras el muro
de pago, son los datos por estudio que dos de los cinco detectores necesitan.
Una evaluación anclada a desenlaces de los cinco necesitaría o bien esas
tablas, o bien un corpus equivalente con efectos y fechas por estudio.

# 7. Limitaciones

Reenunciadas de forma compacta, porque cada una ya ha calificado un resultado
más arriba.

El brazo A son 17 revisiones, seleccionadas por disponibilidad de datos y no
muestreadas, sin nada retenido; sus cifras son exploratorias. Sus objetivos de
evaluación observan el movimiento de la estimación combinada y no ninguna
decisión, de modo que sus tasas miden concordancia con un criterio declarado.
Dos de los cinco detectores se apartan de su procedimiento publicado, y ambos
llevan el nombre de esa desviación. El brazo B alcanza solo la estimación
combinada, dando soporte a tres detectores y no a cinco; su codificación de
desenlaces se realizó una vez, por un codificador automático, y su tasa del
52 % frente a un 9 % publicado sigue sin explicación y se reporta como
provisional. El criterio de comparabilidad entre estimaciones combinadas es un
proxy de similitud textual de la identidad del desenlace, no una coincidencia
verificada.

# 8. Disponibilidad

El paquete es público y tiene licencia MIT en
<https://github.com/jnverbel/staleness>. El brazo A lo regenera
`inst/applicability/applicability.R`, las figuras de calibración
`inst/calibration/calibration.R`, y el brazo B los tres guiones de `corpus/`.
La búsqueda fechada en CRAN que respalda la afirmación del §1 está en
`inst/cran-search/`.

# Referencias
