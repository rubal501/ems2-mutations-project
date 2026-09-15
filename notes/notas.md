
# Título tentativo

**Embeddings locales de protein language models y un baseline clásico BLOSUM para estimar el efecto funcional de mutaciones puntuales en TP53**

---

# Resumen / Abstract

La predicción del efecto funcional de mutaciones missense sigue siendo un problema central en bioinformática. En este trabajo evaluamos si una representación local de mutaciones puntuales en **TP53**, obtenida a partir de **protein language models** de la familia **ESM-2**, captura mejor el efecto funcional de variantes missense que una representación global de toda la proteína y cómo se compara frente a un baseline clásico pre-deep-learning basado en **BLOSUM62**. Para ello, desarrollamos un análisis en dos fases. En la primera, generamos un panel sintético de 228 mutaciones puntuales en 12 posiciones del dominio de unión a DNA de TP53 y cuantificamos la perturbación inducida por cada mutación en el espacio de embeddings. En la segunda, integramos un dataset funcional experimental de MaveDB y evaluamos la asociación entre distintas métricas derivadas de embeddings y el score funcional de los mutantes. Nuestros resultados muestran que el promedio global de embeddings sobre toda la secuencia diluye el efecto de mutaciones puntuales, mientras que una métrica local basada en una ventana alrededor del sitio mutado recupera una señal más fuerte e interpretable. Además, la perturbación local del embedding mostró asociación significativa con el score funcional experimental, y esta asociación mejoró al usar un modelo ESM de mayor capacidad. Sin embargo, la comparación contra BLOSUM62 reveló que no existe un ganador global uniforme: la ventaja relativa de ESM dependió del subconjunto analizado. En particular, ESM superó a BLOSUM en un subconjunto comparativo del dominio de unión a DNA, mientras que BLOSUM fue más competitivo en hotspots canónicos. En conjunto, estos resultados sugieren que el valor de los protein language models no reside en reemplazar automáticamente a los enfoques clásicos, sino en aportar contexto en regiones donde la química local de la sustitución no basta para explicar el efecto funcional.

---

# Introduction

La estimación del efecto funcional de mutaciones missense es un problema de gran relevancia en biología computacional, genómica funcional y medicina de precisión. Una sustitución de un solo aminoácido puede alterar estabilidad, plegamiento, afinidad por ligandos, interacciones proteína-proteína o actividad molecular, pero inferir dicho efecto a partir de la secuencia sigue siendo una tarea difícil. Históricamente, muchos métodos de predicción se han basado en señales evolutivas o fisicoquímicas relativamente simples, como matrices de sustitución aminoacídica, conservación por posición y perfiles derivados de alineamientos múltiples.

Entre estos enfoques clásicos, **BLOSUM62** ocupa un lugar central. Aunque es una representación compacta y conceptualmente simple, captura información evolutiva agregada sobre qué tan aceptables son distintas sustituciones aminoacídicas y, por ello, continúa siendo un baseline sorprendentemente fuerte en problemas de mutación puntual. Sin embargo, BLOSUM no incorpora de manera explícita el contexto secuencial específico en el que ocurre la mutación, lo que limita su capacidad para representar efectos dependientes del entorno local del residuo.

En años recientes, los **protein language models (PLMs)** han ofrecido una alternativa prometedora para representar secuencias proteicas. Modelos como **ESM-2** producen embeddings contextuales por residuo y por secuencia, lo que sugiere que podrían capturar restricciones implícitas de la familia proteica y dependencias locales que los métodos clásicos no modelan directamente. No obstante, el simple uso de embeddings no garantiza una mejora automática. En particular, para el problema de mutaciones puntuales surge una pregunta metodológica inmediata: ¿conviene representar una mutación mediante una descripción **global** de toda la proteína o mediante una representación **local** centrada en el sitio alterado?

**TP53** constituye un sistema ideal para explorar esta pregunta. Es una proteína ampliamente estudiada, con abundantes mutaciones funcionalmente caracterizadas y con datasets experimentales disponibles a gran escala. Además, su dominio de unión a DNA concentra múltiples hotspots mutacionales y posiciones de relevancia estructural y funcional bien documentadas, lo que permite contrastar distintos tipos de sitios dentro de una misma proteína.

En este trabajo nos propusimos comparar un enfoque basado en **embeddings locales de ESM-2** con un baseline clásico basado en **BLOSUM62** para estimar el efecto funcional de mutaciones puntuales en TP53. Planteamos tres preguntas principales. Primero, si una representación local supera a una representación global para describir mutaciones puntuales. Segundo, si la perturbación local del embedding se asocia con scores funcionales experimentales. Tercero, si la ventaja relativa de un protein language model frente a BLOSUM es uniforme o depende del subconjunto de mutaciones considerado. Nuestra hipótesis fue que una mutación puntual en TP53 se reflejaría mejor en una representación local que en una global y que, aunque los embeddings contextuales aportarían información adicional, su ventaja frente a BLOSUM dependería del contexto biológico del sitio mutado.

---

# Methods

## Diseño general del estudio

El estudio se llevó a cabo en dos fases. En la primera fase realizamos un análisis geométrico con mutaciones sintéticas para determinar qué tipo de representación preserva mejor la señal inducida por una sustitución puntual. En la segunda fase integramos un dataset funcional experimental de TP53 para evaluar si las métricas derivadas de embeddings se asocian con función y cómo se comparan con un baseline clásico basado en BLOSUM62.

## Secuencia de referencia y panel de posiciones

Trabajamos con la secuencia canónica humana de **TP53**. Para mantener un análisis controlado y centrado en una región funcionalmente relevante, seleccionamos 12 posiciones dentro del dominio de unión a DNA. Estas posiciones se organizaron en tres subconjuntos:  
1. **Hotspots**: R175, Y220, G245 y R248.  
2. **Structural**: C176, H179, R249 y R282.  
3. **Comparison DBD**: T125, A138, P151 y V157.  

Esta selección buscó equilibrar sitios canónicos y ampliamente estudiados con otros residuos del mismo dominio que funcionaran como grupo comparativo.

## Generación de mutantes sintéticos

Para cada una de las 12 posiciones seleccionadas generamos todas las sustituciones posibles hacia los otros 19 aminoácidos estándar, excluyendo el residuo WT. En total se construyó un panel de **228 mutantes sintéticos**. Cada mutante se etiquetó con la notación convencional `WTposMUT`, por ejemplo `R175H`, y se registró junto con su grupo de pertenencia.

## Modelos de embeddings

Se utilizaron modelos de la familia **ESM-2**. Primero se trabajó con una versión pequeña para implementar el pipeline base y explorar la geometría del problema. Posteriormente se probó una versión con mayor número de parámetros con el fin de evaluar si la capacidad del modelo mejoraba la asociación entre embeddings y scores funcionales experimentales.

## Extracción de embeddings y representaciones globales/locales

A partir de cada secuencia se extrajeron embeddings por residuo. Sobre estas representaciones se definieron varias métricas de mutación mediante comparación entre WT y mutante:

- **Representación global**: promedio del embedding sobre toda la secuencia de TP53.
- **Representación local k=0**: embedding del residuo mutado solamente.
- **Representación local k=3**: promedio en una ventana de ±3 residuos alrededor del sitio mutado.
- **Representación local k=5**: promedio en una ventana de ±5 residuos alrededor del sitio mutado.

La comparación entre WT y mutante se realizó principalmente mediante **distancia coseno**, aunque también se calcularon distancias euclidianas como análisis auxiliar. Tras la exploración inicial, la métrica principal adoptada fue `k5_cosine_distance_to_WT`.

## Análisis geométrico del panel sintético

Para caracterizar la perturbación inducida por las mutaciones sintéticas se realizaron varios análisis descriptivos:

- reducción de dimensión mediante PCA sobre embeddings globales,
- boxplots por grupo de posiciones,
- ranking de sitios según sensibilidad promedio,
- heatmaps posición × aminoácido mutado.

El objetivo de esta fase fue identificar si la señal inducida por mutaciones puntuales se organiza mejor con una descripción global o local y si dicha señal depende más del sitio o del grupo funcional asignado.

## Anotación fisicoquímica de mutaciones

Cada sustitución se anotó mediante un conjunto de descriptores binarios y categóricos, incluyendo:

- cambio de carga,
- cambio de clase fisicoquímica,
- cambio de tamaño,
- cambio de grupo hidrofóbico,
- introducción de glicina,
- introducción de prolina,
- introducción de cisteína,
- cambio de aromaticidad,
- categoría `conservative-like`.

Estas anotaciones se utilizaron para evaluar qué propiedades del cambio aminoacídico se asociaban con mayores perturbaciones locales en el espacio de embeddings.

## Pruebas estadísticas sobre mutaciones sintéticas

Se comparó la distribución de la métrica `k5_cosine_distance_to_WT` entre distintas categorías binarias usando pruebas de **Mann–Whitney**. Los valores p se corrigieron por múltiples pruebas mediante el procedimiento de **Benjamini–Hochberg**, y se calcularon tamaños de efecto con **rank-biserial correlation**.

## Integración de datos funcionales experimentales

Para evaluar el valor biológico de las métricas derivadas de embeddings se integró un score set experimental de **MaveDB** para TP53. Las variantes se normalizaron a una notación uniforme tipo `R175H` y se filtraron para conservar únicamente mutantes missense simples presentes en el panel de interés. El score experimental se interpretó como una medida tipo **deleteriousness-like**, consistente con pérdida de función y/o efecto dominante-negativo. Para facilitar ciertas visualizaciones también se definió una variable reorientada de actividad como `activity_score = -score_raw`.

## Correlaciones con el score funcional

La asociación entre scores experimentales y métricas derivadas de embeddings se evaluó principalmente mediante **correlación de Spearman**, dado que no se asumió una relación estrictamente lineal entre perturbación geométrica y efecto funcional. Como análisis complementario se calcularon también correlaciones de Pearson y Kendall.

## Baseline clásico: BLOSUM62

Como baseline pre-deep-learning se calculó para cada mutación el score correspondiente de **BLOSUM62**. Para alinear la interpretación de ambos métodos respecto del target funcional, se definió `blosum_badness = -BLOSUM62`, de modo que valores mayores reflejaran sustituciones potencialmente más deletéreas.

## Comparación ESM vs BLOSUM

La comparación entre ambos enfoques se realizó de dos maneras. Primero, de forma global, midiendo la correlación entre cada score y el target experimental. Segundo, por subconjuntos, analizando el total de mutantes y varios grupos biológicamente relevantes, incluyendo `hotspot`, `structural` y `comparison_dbd`. Para cuantificar la ventaja relativa de ESM frente a BLOSUM se definió `Δρ = ρ_ESM - ρ_BLOSUM`, y se estimaron intervalos de confianza mediante bootstrap.

---

# Results

## Las representaciones globales diluyen el efecto de mutaciones puntuales

El primer resultado del proyecto fue que la representación global de TP53
resulta poco informativa para estudiar mutaciones puntuales. Cuando se promedió
el embedding sobre toda la secuencia, las distancias entre WT y mutantes fueron
pequeñas y la separación entre los grupos de posiciones fue débil. En otras
palabras, el efecto de una sola sustitución quedó fuertemente diluido por el
tamaño total de la proteína. Este resultado motivó la transición hacia
representaciones locales centradas en el sitio mutado.

## Las representaciones locales recuperan una señal más fuerte e interpretable

Al restringir el análisis a ventanas locales, la perturbación inducida por
mutaciones aumentó de forma clara. La representación `k=0`, basada únicamente
en el residuo mutado, produjo diferencias muy marcadas, pero resultó demasiado
extrema y poco estable como métrica principal. En contraste, las ventanas `k=3`
y `k=5` mostraron una señal más robusta y mejor alineada con la idea de
capturar el contexto secuencial inmediato del sitio mutado. De ellas, `k=5`
emergió como la opción más razonable para continuar el análisis.

## La señal se organiza mejor por sitio que por grupo

Aunque inicialmente se esperaba cierta separación entre grupos gruesos de
posiciones, los análisis mostraron que la señal del embedding se organizaba
mejor a nivel de **sitio individual**. Algunas posiciones presentaron una
sensibilidad elevada a mutación, mientras que otras no, independientemente de
pertenecer a categorías como `hotspot` o `structural`. Este resultado sugirió
que el embedding local captura propiedades específicas del sitio mutado que no
pueden resumirse adecuadamente con etiquetas gruesas de grupo.

## El tipo de sustitución influye en la perturbación local del embedding

Las anotaciones fisicoquímicas mostraron que la perturbación local del
embedding depende del tipo de sustitución. Las mutaciones clasificadas como
`conservative-like` presentaron desplazamientos menores, mientras que los
cambios de clase fisicoquímica tendieron a inducir perturbaciones mayores. En
particular, la introducción de glicina produjo una señal fuerte y consistente.
En cambio, cambios de carga, tamaño e introducción de prolina no mostraron un
patrón igualmente robusto.

## Resultados estadísticos de la fase sintética

Las pruebas de Mann–Whitney confirmaron estos patrones. Después de corregir por
múltiples pruebas, observamos evidencia significativa para `class_change`,
`hydro_group_change`, `introduces_glycine` e `is_conservative_like`. Estos
resultados apoyan la idea de que la métrica local basada en embeddings refleja
principalmente la severidad fisicoquímica de la sustitución y la no
conservatividad del cambio.

## La métrica local se asocia con el score funcional experimental, la global no

Al integrar los scores experimentales de MaveDB, la diferencia entre
representaciones globales y locales se volvió más clara. La métrica global
mostró una asociación débil y no significativa con el score funcional. En
cambio, la métrica local `k=5` mostró una asociación significativa, aunque de
magnitud moderada. Este resultado indica que la representación local conserva
información biológicamente relevante sobre el efecto funcional de mutaciones
puntuales que se pierde al promediar sobre toda la proteína.

## Modelos ESM más grandes mejoran la asociación con función

Cuando se probó una versión de ESM con más parámetros, la correlación entre la
métrica local y el score funcional experimental aumentó sustancialmente. Este
resultado fortaleció dos conclusiones principales: primero, que las
representaciones locales son preferibles a las globales; y segundo, que la
capacidad del modelo influye en la calidad biológica de la representación
obtenida.

## BLOSUM62 fue un baseline más competitivo de lo esperado

La comparación global entre ESM y BLOSUM mostró que la diferencia entre ambos
métodos era menor de lo esperado inicialmente. Lejos de ser un baseline débil,
BLOSUM62 capturó una fracción importante de la señal funcional, lo que sugiere
que en mutaciones missense puntuales la química local de la sustitución ya
contiene mucha información relevante. Esto no invalida el valor de ESM, pero sí
indica que el beneficio de un PLM no debe asumirse automáticamente.

## La ventaja relativa de ESM y BLOSUM depende del subconjunto

El hallazgo más interesante emergió al realizar comparaciones por subconjuntos.
En el conjunto completo de mutantes no hubo evidencia clara de que ESM superara
a BLOSUM. Sin embargo, al separar por grupo, se observó un patrón contrastante.
En el subconjunto `comparison_dbd`, ESM mostró una ventaja clara, con un `Δρ`
positivo y un intervalo bootstrap que no incluyó cero. En cambio, en `hotspot`,
BLOSUM fue superior, con un `Δρ` negativo y un intervalo bootstrap
completamente por debajo de cero. Para el grupo `structural`, la comparación no
fue concluyente.

## Resultado conceptual principal

En conjunto, los resultados indican que la pregunta relevante no es si ESM
“gana” globalmente frente a BLOSUM, sino **en qué subconjuntos el contexto
modelado por ESM aporta información adicional**. En TP53, los embeddings
contextuales fueron especialmente útiles en un subconjunto comparativo del
dominio de unión a DNA, mientras que en hotspots canónicos la severidad química
de la sustitución, bien capturada por BLOSUM, fue suficiente o incluso
superior.

---

# Discussion

Este trabajo muestra que la forma de representar una mutación puntual es tan
importante como el modelo utilizado para generar la representación. En TP53, el
promedio global del embedding sobre toda la secuencia resultó inadecuado, pues
diluyó el efecto de una sola sustitución. En contraste, las representaciones
locales permitieron recuperar una señal más fuerte, más interpretable y más
alineada con el score funcional experimental. Esto confirma que, para proteínas
relativamente largas, el contexto inmediato del residuo mutado debe tener
prioridad sobre descriptores globales cuando se busca estudiar el efecto de
variantes puntuales.

Nuestros resultados también ponen en perspectiva el valor de los protein
language models frente a enfoques clásicos. BLOSUM62 fue un baseline
sorprendentemente competitivo, lo que sugiere que una parte importante de la
señal funcional de mutaciones missense sigue estando contenida en la severidad
química local del cambio aminoacídico. Este hallazgo es importante porque evita
una interpretación simplista según la cual un PLM necesariamente supera a
métodos más antiguos por el mero hecho de ser más complejo. En realidad, la
comparación depende de qué tan bien explota uno la capacidad del modelo y de
qué tan exigente es el problema.

En ese sentido, es importante notar que la métrica basada en ESM utilizada aquí
fue deliberadamente simple: una distancia entre embeddings locales de WT y
mutante. Esto probablemente no agota toda la información disponible en el
modelo. Es razonable pensar que scores basados en log-probabilidades
condicionadas o masked marginals podrían explotar mejor la capacidad del PLM y
ampliar la brecha frente a BLOSUM. Aun así, incluso con esta definición
relativamente austera, ESM mostró asociación significativa con función y mejoró
claramente sobre la versión global, lo que ya representa un resultado
relevante.

El hallazgo más interesante desde el punto de vista biológico y metodológico
fue que la ventaja relativa entre ESM y BLOSUM dependió del subconjunto de
mutaciones. En `comparison_dbd`, ESM superó a BLOSUM, lo que sugiere que el
contexto secuencial aporta información útil en sitios donde la química local de
la sustitución no basta. En cambio, en `hotspot`, BLOSUM fue superior, lo que
sugiere que en esos sitios canónicos la gravedad funcional de una mutación
podría estar más fuertemente determinada por la severidad intrínseca de la
sustitución aminoacídica. Esta observación sugiere una visión más matizada del
problema: los PLMs no reemplazan uniformemente a los enfoques clásicos, sino
que pueden complementarlos en regiones o mutaciones donde el contexto local es
especialmente relevante.

El estudio presenta varias limitaciones. En primer lugar, se centró
exclusivamente en TP53, por lo que no sabemos hasta qué punto las conclusiones
se generalizan a otras proteínas. En segundo lugar, se utilizó una métrica
simple basada en distancia en el espacio de embeddings, sin explotar scores
probabilísticos más ricos. En tercer lugar, aunque BLOSUM es un baseline
importante, no incluimos todavía otros enfoques clásicos como Grantham,
perfiles evolutivos, PSSM o modelos coevolutivos como EVmutation. Por último,
algunos subconjuntos definidos por flags fisicoquímicos fueron pequeños, lo que
limita la precisión de ciertas comparaciones bootstrap.

Como líneas futuras de trabajo, sería natural ampliar la comparación hacia
otros baselines clásicos y hacia scores más fuertes derivados de ESM. También
sería relevante evaluar si la ventaja relativa de ESM se asocia con propiedades
estructurales del sitio mutado, como exposición al solvente, cercanía al DNA o
rol en estabilidad local. Finalmente, extender este análisis a otras proteínas
con datasets MAVE disponibles permitiría determinar si el patrón observado en
TP53 es específico de este sistema o refleja un fenómeno más general.

---

# Conclusion

En este estudio mostramos que las representaciones locales de mutaciones
derivadas de **protein language models** son mucho más informativas que las
representaciones globales para estudiar mutaciones puntuales en **TP53**. La
perturbación local del embedding se asoció significativamente con scores
funcionales experimentales, y esta asociación mejoró al utilizar un modelo ESM
de mayor capacidad. Sin embargo, la comparación frente a **BLOSUM62** reveló
que no existe una superioridad global uniforme de ESM. En el conjunto completo
de mutantes no hubo un ganador claro, pero la comparación por subconjuntos
mostró que **ESM supera a BLOSUM en ciertas regiones del dominio de unión a
DNA**, mientras que **BLOSUM sigue siendo más competitivo en hotspots
canónicos**.

Por lo tanto, el principal mensaje del trabajo no es que los protein language
models sustituyan automáticamente a los enfoques clásicos, sino que **su
utilidad depende del contexto biológico del sitio mutado**. En TP53, los
embeddings contextuales aportan valor especialmente en subconjuntos donde la
química local de la sustitución no basta para explicar el efecto funcional,
mientras que los métodos clásicos permanecen notablemente fuertes en posiciones
altamente canónicas.

---

# Texto corto para cada figura

## Figura 1. Esquema general del estudio
Diagrama del pipeline completo: selección de TP53, generación del panel sintético, extracción de embeddings ESM, definición de métricas globales y locales, integración con scores funcionales de MaveDB y comparación frente a BLOSUM62.

## Figura 2. Comparación entre representación global y local
Gráfica que muestra que el promedio global de embeddings diluye el efecto de mutaciones puntuales, mientras que la representación local `k=5` recupera una señal más clara.

## Figura 3. Efecto del tipo de sustitución sobre la perturbación local
Boxplots o resúmenes estadísticos para `class_change`, `introduces_glycine` e `is_conservative_like`, mostrando que la perturbación local depende del tipo fisicoquímico de la mutación.

## Figura 4. Asociación entre embedding local y score funcional
Scatter plot entre `k5_cosine_distance_to_WT` y el score experimental, mostrando una asociación significativa y su mejora al usar un modelo ESM de mayor capacidad.

## Figura 5. Comparación global entre ESM y BLOSUM
Figura o tabla resumen mostrando que BLOSUM62 es un baseline competitivo y que la diferencia global frente a ESM es menor a la esperada.

## Figura 6. Ventaja relativa de ESM vs BLOSUM por subconjunto
Gráfica bootstrap de `Δρ = ρ_ESM - ρ_BLOSUM` para `all_mutants`, `comparison_dbd` y `hotspot`, destacando el carácter subconjunto-dependiente de la comparación.

---

# Texto corto para cada tabla

## Tabla 1. Panel sintético de posiciones en TP53
Resumen de las 12 posiciones seleccionadas, aminoácido WT, grupo asignado y justificación biológica.

## Tabla 2. Resultados estadísticos del análisis por tipo de mutación
Comparación de `k5_cosine_distance_to_WT` entre flags fisicoquímicos, incluyendo tamaño de muestra, medianas, p ajustado y tamaño de efecto.

## Tabla 3. Correlaciones entre métricas derivadas de embeddings y score funcional
Correlaciones de Spearman, Pearson y Kendall para las representaciones global y locales, con comparación entre modelo pequeño y modelo grande.

## Tabla 4. Comparación entre ESM y BLOSUM por subconjunto
Valores de `Δρ`, tamaños de muestra e intervalos bootstrap para `all_mutants`, `comparison_dbd` y `hotspot`.

---

# Cierre operativo

Este borrador ya puede servir como base para empezar a redactar el reporte completo. Lo más probable es que, al pasar a la versión final, sólo necesites:

1. insertar los valores numéricos definitivos en las secciones de Results,  
2. añadir referencias bibliográficas,  
3. incorporar las figuras y tablas en el lugar correspondiente,  
4. pulir el tono según el formato del curso, reporte o artículo.
