#import "@preview/elsearticle:3.1.0": *

#show: elsearticle.with(
  title: [_Embeddings_ locales de _protein language models_ y un baseline clásico BLOSUM
    para estimar el efecto funcional de mutaciones puntuales en TP53],
  abstract: [
    La predicción del efecto funcional de mutaciones missense constituye un problema
    central en bioinformática. En este trabajo se evalúa si una representación local
    de mutaciones puntuales en TP53, obtenida a partir de modelos de lenguaje
    proteico de la familia ESM-2, captura mejor el efecto funcional de variantes
    missense que una representación global de toda la proteína, y cómo se compara
    frente a un baseline clásico basado en BLOSUM62. El análisis se llevó a cabo en
    dos fases. En la primera se generó un panel sintético de 228 mutaciones puntuales
    en 12 posiciones del dominio de unión a DNA de TP53 y se cuantificó la perturbación
    inducida en el espacio de embeddings. En la segunda se integraron scores funcionales
    experimentales de MaveDB y se evaluó la asociación entre distintas métricas
    derivadas de embeddings y el efecto funcional de los mutantes. Los resultados
    muestran que el promedio global diluye el efecto de mutaciones puntuales, mientras
    que una métrica local centrada en el sitio mutado recupera una señal más fuerte e
    interpretable. Además, la perturbación local del embedding mostró asociación
    significativa con el score funcional experimental, y esta asociación mejoró al
    emplear un modelo ESM de mayor capacidad. La comparación frente a BLOSUM62 reveló
    que la ventaja relativa de ESM depende del subconjunto analizado: ESM superó a
    BLOSUM en un subconjunto comparativo del dominio de unión a DNA, mientras que
    BLOSUM fue más competitivo en hotspots canónicos. En conjunto, estos resultados
    sugieren que el valor de los modelos de lenguaje proteico no reside en reemplazar
    automáticamente a los enfoques clásicos, sino en aportar contexto en regiones
    donde la química local de la sustitución no basta para explicar el efecto funcional.
  ],
  authors: (
    (
      name: "Javier Roberto Rubalcava Cortes",
      email: "rubal820[at]ciencias[dot]unam[dot]mx",
    ),
  ),
)

// ─────────────────────────────────────────────────────────────────────────────
= Introducción
// ─────────────────────────────────────────────────────────────────────────────

La estimación del efecto funcional de mutaciones missense es un problema de gran
relevancia en biología computacional, genómica funcional y medicina de precisión.
Una sustitución de un solo aminoácido puede alterar la estabilidad, el plegamiento,
la afinidad por ligandos, las interacciones proteína–proteína o la actividad
molecular; sin embargo, inferir dicho efecto a partir de la secuencia sigue siendo
una tarea difícil. Históricamente, muchos métodos de predicción se han apoyado en
señales evolutivas o fisicoquímicas relativamente simples: matrices de sustitución
aminoacídica, conservación por posición y perfiles derivados de alineamientos
múltiples de secuencias.

Entre estos enfoques clásicos, BLOSUM62 ocupa un lugar central. Aunque es una
representación compacta y conceptualmente simple, captura información evolutiva
agregada sobre qué tan aceptables son distintas sustituciones aminoacídicas y, por
ello, continúa siendo un baseline sorprendentemente fuerte en problemas de mutación
puntual. No obstante, BLOSUM62 no incorpora de manera explícita el contexto
secuencial específico en el que ocurre la mutación, lo que limita su capacidad para
representar efectos dependientes del entorno local del residuo.

En años recientes, los modelos de lenguaje proteico (_protein language models_,
PLMs) han ofrecido una alternativa prometedora para representar secuencias
proteicas. Modelos como ESM-2 producen embeddings contextuales por residuo y por
secuencia, lo que sugiere que podrían capturar restricciones implícitas de la
familia proteica y dependencias locales que los métodos clásicos no modelan
directamente. No obstante, el uso de embeddings no garantiza por sí solo una mejora
automática. En particular, para el problema de mutaciones puntuales surge una
pregunta metodológica inmediata: ¿conviene representar una mutación mediante una
descripción global de toda la proteína o mediante una representación local centrada
en el sitio alterado?

TP53 constituye un sistema ideal para explorar esta pregunta. Es el gen más
frecuentemente mutado en cáncer humano @hainaut2016, con abundantes mutaciones
funcionalmente caracterizadas mediante ensayos de mutagénesis en paralelo a gran
escala (_multiplex assays of variant effect_, MAVE). En particular, el estudio de
Giacomelli et al. generó scores funcionales para más de 8,000 alelos de TP53 en
líneas celulares isogénicas mediante pantallas de selección fenotípica, y estos
datos han sido depositados en MaveDB @giacomelli2018. Además, el dominio de unión
a DNA de TP53 concentra múltiples hotspots mutacionales y posiciones de relevancia
estructural y funcional bien documentadas, lo que permite contrastar distintos
tipos de sitios dentro de una misma proteína.

En este trabajo se comparó un enfoque basado en embeddings locales de ESM-2 con un
baseline clásico basado en BLOSUM62 para estimar el efecto funcional de mutaciones
puntuales en TP53. Se plantearon tres preguntas principales: (1) si una
representación local supera a una representación global para describir mutaciones
puntuales; (2) si la perturbación local del embedding se asocia con scores
funcionales experimentales; y (3) si la ventaja relativa de un PLM frente a BLOSUM62
es uniforme o depende del subconjunto de mutaciones considerado. La hipótesis de
partida fue que una mutación puntual en TP53 se reflejaría mejor en una
representación local que en una global y que, aunque los embeddings contextuales
aportarían información adicional, su ventaja frente a BLOSUM62 dependería del
contexto biológico del sitio mutado.

// ─────────────────────────────────────────────────────────────────────────────
= Métodos
// ─────────────────────────────────────────────────────────────────────────────

== Diseño general del estudio

El estudio se llevó a cabo en dos fases. En la primera se ejecutó un análisis
geométrico con mutaciones sintéticas para determinar qué tipo de representación
preserva mejor la señal inducida por una sustitución puntual. En la segunda se
integraron datos funcionales experimentales de TP53 para evaluar si las métricas
derivadas de embeddings se asocian con función y cómo se comparan con un baseline
clásico basado en BLOSUM62.

== Secuencia de referencia y panel de posiciones

Se trabajó con la secuencia canónica humana de TP53. Para mantener un análisis
controlado y centrado en una región funcionalmente relevante, se seleccionaron
12 posiciones dentro del dominio de unión a DNA (DBD) @cho1994, organizadas en tres
subconjuntos con distinto nivel de caracterización funcional conocida:

- *_Hotspots_ canónicos* (R175, Y220, G245, R248): posiciones con alta frecuencia
  de mutación en bases de datos de cáncer y cuya pérdida de función está
  extensamente documentada @giacomelli2018.
- *Estructurales* (C176, H179, R249, R282): residuos críticos para la integridad
  estructural del DBD, entre ellos coordinadores del ion Zn#super[2+] y contactos
  de la lámina beta central @bauer2020@leroy2017.
- *Comparativos del DBD* (T125, A138, P151, V157): posiciones dentro del mismo
  dominio pero sin anotación funcional consolidada en la literatura. Se incluyeron
  como controles negativos internos: si las métricas derivadas de embeddings
  mostraran perturbaciones equivalentes en este subconjunto y en los _hotspots_,
  ello indicaría que el método carece de especificidad funcional. Su presencia
  permite además evaluar si la señal del embedding refleja relevancia biológica
  o simplemente variación fisicoquímica genérica de la sustitución.

Esta estratificación permitió contrastar distintos niveles de importancia funcional
dentro de una misma región estructural, controlando así la contribución del dominio
como variable de confusión.

== Generación de mutantes sintéticos

Para cada una de las 12 posiciones seleccionadas se generaron todas las
sustituciones posibles hacia los otros 19 aminoácidos estándar, excluyendo el
residuo de tipo silvestre (WT). En total se construyó un panel de 228 mutantes
sintéticos. Cada mutante se etiquetó con la notación convencional `WTposMUT`
(p. ej., R175H) y se registró junto con su grupo de pertenencia.

== Modelos de embeddings

Se utilizaron dos modelos de la familia ESM-2 @lin2023. El primero,
`esm2_t33_650M_UR50D` (650 millones de parámetros, 33 capas de transformer),
se empleó para implementar el pipeline base y explorar la geometría del problema.
El segundo, `esm2_t36_3B_UR50D` (3 mil millones de parámetros, 36 capas), se
utilizó para evaluar si una mayor capacidad del modelo mejoraba la asociación
entre embeddings y scores funcionales experimentales. En ambos casos se extrajeron
los embeddings de la última capa del transformer, siguiendo la práctica habitual
de usar representaciones de las capas finales como descriptores de secuencia
@meier2021.

== Extracción de embeddings y definición de representaciones

A partir de cada secuencia se extrajeron embeddings por residuo. Sobre estas
representaciones se definieron cuatro métricas de mutación mediante comparación
entre WT y mutante:

- *Representación global*: promedio del embedding sobre toda la secuencia de TP53.
- *Representación local k = 0*: embedding del residuo mutado únicamente.
- *Representación local k = 3*: promedio en una ventana de ±3 residuos alrededor
  del sitio mutado.
- *Representación local k = 5*: promedio en una ventana de ±5 residuos alrededor
  del sitio mutado.

La comparación entre WT y mutante se efectuó principalmente mediante distancia
coseno, elegida por la mayor estabilidad numérica que ofrece en espacios latentes
de alta dimensionalidad frente a métricas basadas en norma. Como análisis auxiliar
también se calcularon distancias euclidianas. Tras la exploración inicial, la
métrica principal adoptada fue `k5_cosine_distance_to_WT`.

== Análisis geométrico del panel sintético

Para caracterizar la perturbación inducida por las mutaciones sintéticas se
efectuaron varios análisis descriptivos: reducción de dimensión mediante análisis
de componentes principales (PCA, del inglés _principal component analysis_) sobre
embeddings globales, diagramas de caja por grupo de posiciones, ranking de sitios
según sensibilidad promedio y mapas de calor posición × aminoácido mutado. El
objetivo de esta fase fue determinar si la señal inducida por mutaciones puntuales
se organizaba mejor con una descripción global o local, y si dicha señal dependía más
del sitio o del grupo funcional asignado.

== Anotación fisicoquímica de mutaciones

Cada sustitución se anotó mediante un conjunto de descriptores binarios y
categóricos que incluyó: cambio de carga, cambio de clase fisicoquímica, cambio de
tamaño, cambio de grupo hidrofóbico, introducción de glicina, introducción de
prolina, introducción de cisteína, cambio de aromaticidad y categoría
`conservative-like`. Estas anotaciones se emplearon para evaluar qué propiedades
del cambio aminoacídico se asocian con mayores perturbaciones locales en el espacio
de embeddings.

== Pruebas estadísticas sobre mutaciones sintéticas

Para comparar la métrica `k5_cosine_distance_to_WT` entre los tres grupos de
posiciones se aplicó inicialmente una prueba de ANOVA de un factor. Dado que las
distribuciones no cumplían los supuestos de normalidad y homocedasticidad, se empleó
la prueba no paramétrica de Kruskal-Wallis. Las comparaciones entre pares de grupos
y las de tipo «uno contra todos» se llevaron a cabo mediante pruebas de
Mann-Whitney. Adicionalmente, se evaluó el efecto de cada característica
fisicoquímica binaria sobre la distancia coseno con pruebas de Mann-Whitney. En
todos los casos, los valores $p$ se corrigieron por comparaciones múltiples mediante
el método de Holm, y los tamaños de efecto se cuantificaron con la delta de Cliff.

== Integración de datos funcionales experimentales

Para evaluar el valor biológico de las métricas derivadas de embeddings se integró
el conjunto de scores funcionales experimentales de MaveDB con identificador
`urn:mavedb:00000068-0-1`, derivado del estudio de mutagénesis por saturación de
TP53 publicado por Giacomelli et al. @giacomelli2018. En dicho estudio se
construyó una biblioteca de 8,258 alelos mutantes de TP53 y se introdujo en dos
líneas celulares isogénicas de A549: una con TP53 de tipo silvestre (p53#super[WT])
y una con TP53 eliminado (p53#super[NULL]). Las células se sometieron a tres
condiciones de selección —p53#super[WT] + nutlin-3, p53#super[NULL] + nutlin-3 y
p53#super[NULL] + etopósido— durante 12 días, y la abundancia de cada alelo se
cuantificó mediante secuenciación masiva en paralelo. El score funcional combinado
de cada variante se calculó a partir del cambio en log#sub[2] de las lecturas
relativas al punto de tiempo inicial, promediado y estandarizado en las tres
condiciones. Valores negativos del score corresponden a alelos con actividad
dominante-negativa (DN) o de pérdida de función (LOF), mientras que valores
positivos son compatibles con actividad similar al WT.

Para la integración en el presente análisis, las variantes del dataset se
normalizaron a una notación uniforme tipo `R175H` y se filtraron para conservar
únicamente mutantes missense simples presentes en el panel de interés. El score
experimental se interpretó como una medida de tipo _deleteriousness-like_,
consistente con pérdida de función y/o efecto dominante-negativo. Para facilitar
ciertas visualizaciones se definió además una variable reorientada de actividad
como `activity_score = −score_raw`.

== Correlaciones con el score funcional

La asociación entre los scores experimentales y las métricas derivadas de embeddings
se evaluó principalmente mediante correlación de Spearman, dado que no se asumió
una relación estrictamente lineal entre perturbación geométrica y efecto funcional.
Como análisis complementario se calcularon también correlaciones de Pearson y Kendall.

== Baseline clásico: BLOSUM62

Como baseline pre-_deep learning_ se calculó para cada mutación el score
correspondiente de BLOSUM62. Para alinear la interpretación de ambos métodos
respecto del target funcional, se definió `blosum_badness = −BLOSUM62`, de modo
que valores mayores reflejaran sustituciones potencialmente más deletéreas.

== Comparación entre ESM-2 y BLOSUM62

La comparación entre ambos enfoques se llevó a cabo de dos maneras. Primero, de forma
global, midiendo la correlación de cada score con el target experimental. Segundo,
por subconjuntos, analizando el total de mutantes y varios grupos biológicamente
relevantes: _hotspot_, _structural_ y _comparison DBD_. Para cuantificar la ventaja
relativa de ESM-2 frente a BLOSUM62 se definió la diferencia entre correlaciones de
Spearman
$Delta rho = rho_"ESM" - rho_"BLOSUM"$,
y se estimaron intervalos de confianza mediante bootstrap.

// ─────────────────────────────────────────────────────────────────────────────
= Resultados
// ─────────────────────────────────────────────────────────────────────────────
== Las representaciones globales diluyen el efecto de mutaciones puntuales

El primer resultado del proyecto fue que la representación global de TP53
resultó poco informativa para estudiar mutaciones puntuales. Cuando se promedió
el embedding sobre toda la secuencia, las distancias entre WT y mutantes fueron
pequeñas y la separación entre los grupos de posiciones fue débil. En otras
palabras, el efecto de una sola sustitución quedó fuertemente diluido por el
tamaño total de la proteína. Este resultado motivó la transición hacia
representaciones locales centradas en el sitio mutado.

#figure(
  image("../figures/compare_rank.svg"),
  caption: [Cambio de sensibilidad de los embeddings al pasar de representaciones globales a representaciones locales.
    Se muestra el desplazamiento promedio en el espacio latente para cada residuo.],
)<fig:rank_k5>

== Las representaciones locales recuperan una señal más fuerte e interpretable

Al restringir el análisis a ventanas locales, la perturbación inducida por
mutaciones aumentó de forma clara. La representación `k=0`, basada únicamente
en el residuo mutado, produjo diferencias muy marcadas, pero resultó demasiado
extrema y poco estable como métrica principal. En contraste, las ventanas `k=3`
y `k=5` mostraron una señal más robusta y mejor alineada con la idea de
capturar el contexto secuencial inmediato del sitio mutado. De ellas, `k=5`
emergió como la opción más razonable para continuar el análisis,
como se observa en la figura @fig:rank_k5.

#figure(
  image("../figures/metric_comparison_boxplot.svg"),
  caption: [ Boxplots de las distancias coseno de los distintos tipos de embeddings, globales y locales],
)

== Las representaciones locales preservan la estructura funcional del espacio de mutaciones

Al aplicar reducción de dimensionalidad mediante PCA a las representaciones
generadas por el modelo de 650M de parámetros, las secuencias se agruparon
de acuerdo con el tipo de posición: las mutaciones en sitios _hotspot_ se
concentraron en regiones del espacio distintas a las de los sitios
_structural_, con separación visible entre grupos. De manera notable,
las mutaciones en las posiciones _comparativas del DBD_ se agruparon en
torno a la representación de la secuencia de tipo silvestre (WT), lo que
es consistente con su menor perturbación funcional esperada. Este patrón
puede observarse en la figura @fig:pca_embed.

#figure(
  scope: "parent",
  placement: auto,
  image("../figures/tsne-global-embeddings.svg"),
  caption: [
    Reducción de dimensionalidad de las representaciones de las secuencias sintéticas
    mediante PCA a dos dimensiones, empleando el modelo de 33 capas y 650M de parámetros.
    Se puede observar cómo la representación local (k=5) muestra agrupaciones más concentradas
    en comparación con la representación global. Asimismo, las secuencias del grupo
    _comparativo del DBD_ se agrupan en torno a la representación de la secuencia
    de tipo silvestre (WT).
  ],
)<fig:pca_embed>

Además de la reducción de dimensionalidad, se aplicaron pruebas de hipótesis sobre
las distribuciones de distancia coseno obtenidas con la representación local $k=5$.
En primer lugar, las mutaciones se separaron según el grupo al que pertenece la
posición mutada. Se efectuó una prueba de ANOVA, que arrojó un valor
$p = 5.96 dot 10^{-8}$; sin embargo, las distribuciones de las distancias no
cumplían los supuestos de ANOVA, por lo que se recurrió a una prueba de
Kruskal-Wallis. Esta última obtuvo un valor $p = 6.7 dot 10^{-6}$ y confirmó una
diferencia significativa entre los tres grupos de mutaciones.

Posteriormente, se compararon los grupos entre sí mediante pruebas de Mann-Whitney
pareadas, cuyos resultados se presentan en la @table:results_mw1. También se
llevaron a cabo pruebas de tipo _uno contra todos_ para cada grupo
(@table:results_mw2). Por último, se emplearon las anotaciones fisicoquímicas de
las mutaciones y se calculó el delta de Cliff para cuantificar el efecto de cada
tipo de cambio sobre la distancia coseno (@table:results_mw3).




#figure(
  table(
    align: center,
    columns: 6,
    table.header[*Grupo 1*][*Grupo 2*][*W*][*$p$ ajustado*][*Cliff's $delta$*][*Magnitud*],
    [Hotspots canónicos], [Estructurales], [3364], [7.974070e-02], [0.1648199], [small],
    [Hotspots canónicos], [Comparativos del DBD], [4099], [2.451583e-05], [0.4193213], [medium],
    [Estructurales], [Comparativos del DBD], [3864], [6.496151e-04], [0.3379501], [medium],
  ),
  caption: [
    Resultados de pruebas de Mann-Whitney a pares entre los distintos grupos de posición de mutaciones.
    Se reportan el estadístico $W$, el $p$-valor ajustado y el tamaño de efecto (Cliff's $delta$).
  ],
)<table:results_mw1>


#figure(
  table(
    align: center,
    columns: 5,
    table.header[*Grupo*][*W*][*$p$ ajustado*][*Cliff's $delta$*][*Magnitud*],
    [Hotspots canónicos], [7463], [6.564183e-04], [0.2920706], [small],
    [Estructurales], [6276], [2.873973e-01], [0.0865651], [negligible],
    [Comparativos del DBD], [3589], [9.631762e-06], [-0.3786357], [medium],
  ),
  caption: [
    Resultados de pruebas de Mann-Whitney (uno contra todos) entre los distintos grupos de posición de mutaciones.
    Se reportan el estadístico $W$, el $p$-valor ajustado y el tamaño de efecto (Cliff's $delta$).
  ],
)<table:results_mw2>


#figure(
  table(
    columns: 5,
    table.header[*Característica*][*W*][*$p$ ajustado*][*Cliff's $delta$*][*Magnitud*],
    [Sustitución conservativa], [346], [0.0002353704], [-0.7330247], [large],
    [Elimina glicina], [3078], [0.0008011229], [0.5502392], [large],
    [Cambio de clase], [5953], [0.0052617509], [0.3189321], [small],
    [Cambio de grupo hidrofóbico], [5953], [0.0052617509], [0.3189321], [small],
  ),
  caption: [
    Principales resultados de pruebas de Mann-Whitney (uno contra todos) para características fisicoquímicas de las mutaciones.
    Se reportan el estadístico $W$, el $p$-valor ajustado y el tamaño de efecto (Cliff's $delta$).
    La tabla completa con las 12 características se encuentra en el Anexo (Tabla A3).
  ],
)<table:results_mw3>





== La métrica local se asocia con el score funcional experimental, la global no

Al integrar los scores experimentales de MaveDB, la diferencia entre
representaciones globales y locales se volvió más clara. La métrica global
mostró una asociación débil y no significativa con el score funcional. En
cambio, la métrica local `k=5` mostró una asociación significativa, aunque de
magnitud moderada. Este resultado indicó que la representación local conservó
información biológicamente relevante sobre el efecto funcional de mutaciones
puntuales que se perdió al promediar sobre toda la proteína.

#figure(
  image("../figures/embedding-vs-functional-score-global-k5.svg"),
  caption: [Scatter plots del score funcional contra el desplazamiento (_shift_) en el espacio latente
    para la representación global (arriba) y la representación local (abajo). Se reportan las correlaciones
    entre el desplazamiento y el score funcional para ambos casos junto con su valor p; la correlación
    es mayor en la representación local y el valor p es menor],
)

== Modelos ESM más grandes mejoran la asociación con función

Cuando se probó una versión de ESM con más parámetros, la correlación entre la
métrica local y el score funcional experimental aumentó sustancialmente. Este
resultado fortaleció dos conclusiones principales: primero, que las
representaciones locales eran preferibles a las globales; y segundo, que la
capacidad del modelo influyó en la calidad biológica de la representación
obtenida.

#figure(
  image("../figures/correlation-models.svg"),
  caption: [Correlaciones entre distintos modelos de la familia ESM-2. En la figura de arriba
    se muestra el diagrama de dispersión de la representación local (_k_=5) para
    el modelo de 6 capas y 8M de parámetros; abajo, el del modelo de 33 capas y
    650M de parámetros. Se observa cómo la correlación es mayor para el modelo de mayor capacidad.
  ],
)

#figure(
  table(
    align: center,
    columns: 3,
    table.header[Método][Spearman $rho$][$p"-valor"$],
    [BLOSUM62], [0.298836], [0.000004],
    [Distancia coseno global], [0.234914], [0.000346],
    [Distancia coseno local (k=5)], [*0.312365*], [0.000001],
  ),
  caption: [Correlaciones de Spearman para los scores calculados con cada uno de los métodos contra
    el score funcional experimental. Se observa cómo las representaciones ESM-2 locales (k=5)
    tienen una mayor correlación y un valor p menor.
  ],
)

== BLOSUM62 fue un baseline más competitivo de lo esperado

La comparación global entre ESM y BLOSUM mostró que la diferencia entre ambos
métodos era menor de lo esperado inicialmente. Lejos de ser un baseline débil,
BLOSUM62 capturó una fracción importante de la señal funcional, lo que sugirió
que en mutaciones missense puntuales la química local de la sustitución ya
contenía mucha información relevante. Esto no invalidó el valor de ESM, pero sí
indicó que el beneficio de un PLM no debe asumirse automáticamente.

== La ventaja relativa de ESM y BLOSUM depende del subconjunto

El principal hallazgo surgió al efectuar las comparaciones por subconjuntos.
En el conjunto completo de mutantes no hubo evidencia clara de que ESM superara
a BLOSUM. Sin embargo, al separar por grupo, se observó un patrón contrastante.
En el subconjunto _comparison\_dbd_, ESM mostró una ventaja clara, con un $Delta rho$
positivo y un intervalo bootstrap que no incluyó cero. En cambio, en _hotspot_,
BLOSUM fue superior, con un $Delta rho$ negativo y un intervalo bootstrap
completamente por debajo de cero, como se observa en la figura @fig:boots_group.

#figure(
  image("../figures/tp53_esm_vs_blosum_selected_subsets.svg"),
  caption: [Intervalos de confianza de la diferencia entre correlaciones de Spearman ($Delta rho$) para ESM-2 frente a BLOSUM62, estimados mediante bootstrap por subconjunto de posiciones],
)<fig:boots_group>



// ─────────────────────────────────────────────────────────────────────────────
= Discusión
// ─────────────────────────────────────────────────────────────────────────────

Este trabajo muestra que la forma de representar una mutación puntual es tan
importante como el modelo utilizado para generar la representación. En TP53, el
promedio global del embedding sobre toda la secuencia resultó inadecuado, pues
diluyó el efecto de una sola sustitución. En contraste, las representaciones
locales permitieron recuperar una señal más fuerte, más interpretable y más
alineada con el score funcional experimental. Esto confirma que, para proteínas
relativamente largas, el contexto inmediato del residuo mutado debe tener
prioridad sobre descriptores globales cuando se busca estudiar el efecto de
variantes puntuales.

Estos resultados también ponen en perspectiva el valor de los _protein language
models_ frente a enfoques clásicos. BLOSUM62 fue un baseline sorprendentemente
competitivo, lo que sugiere que una parte importante de la señal funcional de
mutaciones missense sigue estando contenida en la severidad química local del
cambio aminoacídico. Este hallazgo es relevante porque evita una interpretación
simplista según la cual un PLM necesariamente supera a métodos más antiguos por
el mero hecho de ser m/ás complejo. En realidad, la comparación depende de qué
tan bien se explota la capacidad del modelo y de qué tan exigente es el problema.

En ese sentido, cabe notar que la métrica basada en ESM-2 empleada aquí fue
deliberadamente simple: una distancia entre embeddings locales del WT y el
mutante. Esta definición probablemente no agota toda la información disponible
en el modelo. Es razonable pensar que scores basados en log-probabilidades
condicionadas o en _masked marginals_ (probabilidades marginales bajo
enmascaramiento de residuos) podrían explotar mejor la capacidad del PLM y
ampliar la brecha frente a BLOSUM62. Aun así, incluso con esta formulación
austera, ESM-2 mostró asociación significativa con función y mejoró claramente
sobre la versión global, lo que constituye ya un resultado informativo.

El hallazgo más interesante desde el punto de vista biológico y metodológico
fue que la ventaja relativa entre ESM-2 y BLOSUM62 dependió del subconjunto de
mutaciones. En el grupo _comparison\_dbd_, ESM-2 superó a BLOSUM62, lo que
sugiere que el contexto secuencial aporta información útil en sitios donde la
química local de la sustitución no basta. En cambio, en el grupo _hotspot_,
BLOSUM62 fue superior, lo que sugiere que en esos sitios canónicos la gravedad
funcional de una mutación podría estar más fuertemente determinada por la
severidad intrínseca de la sustitución aminoacídica. Esta observación apunta a
una visión más matizada del problema: los PLMs no reemplazan uniformemente a los
enfoques clásicos, sino que pueden complementarlos en regiones donde el contexto
local es especialmente relevante.

El estudio presenta varias limitaciones. En primer lugar, el análisis se centró
exclusivamente en TP53, por lo que no es posible determinar hasta qué punto las
conclusiones se generalizan a otras proteínas. En segundo lugar, se empleó una
métrica simple basada en distancia en el espacio de embeddings, sin explotar
scores probabilísticos más ricos. En tercer lugar, aunque BLOSUM62 es un
baseline importante, no se incluyeron otros enfoques clásicos como Grantham,
perfiles evolutivos, PSSM o modelos coevolutivos como EVmutation. Por último,
algunos subconjuntos definidos por propiedades fisicoquímicas fueron pequeños,
lo que limita la precisión de ciertas comparaciones bootstrap.

Como líneas futuras de trabajo, sería natural ampliar la comparación hacia otros
baselines clásicos y hacia scores más ricos derivados de ESM-2. También sería
relevante evaluar si la ventaja relativa de ESM-2 se asocia con propiedades
estructurales del sitio mutado, como exposición al solvente, cercanía al DNA o
rol en la estabilidad local. Finalmente, extender este análisis a otras proteínas
con datos MAVE disponibles permitiría determinar si el patrón observado en TP53
es específico de este sistema o refleja un fenómeno más general.

---

= Conclusión

Los resultados de este estudio muestran que las representaciones locales de
mutaciones derivadas de _protein language models_ son considerablemente más
informativas que las representaciones globales para estudiar mutaciones puntuales
en TP53. La perturbación local del embedding se asoció significativamente con
scores funcionales experimentales, y esta asociación mejoró al emplear un modelo
ESM-2 de mayor capacidad. Sin embargo, la comparación frente a BLOSUM62 reveló
que no existe una superioridad global uniforme de ESM-2. En el conjunto completo
de mutantes no hubo un ganador claro, pero la comparación por subconjuntos mostró
que ESM-2 supera a BLOSUM62 en el grupo _comparison\_dbd_, mientras que BLOSUM62
sigue siendo más competitivo en el grupo _hotspot_.

Por lo tanto, el mensaje principal de este trabajo no es que los _protein language
models_ sustituyan automáticamente a los enfoques clásicos, sino que su utilidad
depende del contexto biológico del sitio mutado. En TP53, los embeddings
contextuales aportan valor especialmente en subconjuntos donde la química local
de la sustitución no basta para explicar el efecto funcional, mientras que los
métodos clásicos permanecen notablemente competitivos en posiciones altamente
canónicas.



#bibliography("references.bib")


