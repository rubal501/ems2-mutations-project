#set page(paper: "a4", flipped: true, margin: 12mm)
#set text(size: 8pt, font: "Libertinus Serif")

// Make header row bold
#show table.cell.where(y: 0): set text(weight: "bold")

= Anexo: Figuras extras


#figure(
  image("../figures/violin-plots.svg"),
  caption: [
    Distribucion de las distancias coseno en los distintos grupos de mutaciones.
  ],
)



= Anexo: Tablas completas de pruebas estadísticas

== Tabla A1. Mann-Whitney a pares entre grupos de posición de mutaciones <table:anexo_mw1>

#figure(
  table(
    align: center,
    columns: 12,
    table.header[group_1][group_2][n_1][n_2][median_1][median_2][W][p_value][cliff_delta][cliff_magnitude][p_adjusted][direction],
    [Hotspots canónicos],
    [Estructurales],
    [76],
    [76],
    [0.003640831],
    [0.003335714],
    [3364],
    [7.974070e-02],
    [0.1648199],
    [small],
    [7.974070e-02],
    [Hotspots canónicos > Estructurales],

    [Hotspots canónicos],
    [Comparativos del DBD],
    [76],
    [76],
    [0.003640831],
    [0.002249002],
    [4099],
    [8.171943e-06],
    [0.4193213],
    [medium],
    [2.451583e-05],
    [Hotspots canónicos > Comparativos del DBD],

    [Estructurales],
    [Comparativos del DBD],
    [76],
    [76],
    [0.003335714],
    [0.002249002],
    [3864],
    [3.248075e-04],
    [0.3379501],
    [medium],
    [6.496151e-04],
    [Estructurales > Comparativos del DBD],
  ),
  caption: [
    Resultados completos de pruebas de Mann-Whitney a pares entre los grupos de posición de mutaciones.
    $n_i$: tamaño muestral, $W$: estadístico de Mann-Whitney, Cliff's $delta$: tamaño de efecto.
  ],
)

#pagebreak()

== Tabla A2. Mann-Whitney uno contra todos entre grupos de posición de mutaciones <table:anexo_mw2>

#figure(
  table(
    align: center,
    columns: 10,
    table.header[group][n_group][n_rest][median_group][median_rest][W][p_value][cliff_delta][cliff_magnitude][p_adjusted],
    [Hotspots canónicos],
    [76],
    [152],
    [0.003640831],
    [0.002776474],
    [7463],
    [3.282091e-04],
    [0.2920706],
    [small],
    [6.564183e-04],

    [Estructurales],
    [76],
    [152],
    [0.003335714],
    [0.002775788],
    [6276],
    [2.873973e-01],
    [0.0865651],
    [negligible],
    [2.873973e-01],

    [Comparativos del DBD],
    [76],
    [152],
    [0.002249002],
    [0.003445983],
    [3589],
    [3.210587e-06],
    [-0.3786357],
    [medium],
    [9.631762e-06],
  ),
  caption: [
    Resultados completos de pruebas de Mann-Whitney uno contra todos entre los grupos de posición de mutaciones.
    $n_"group"$/$n_"rest"$: tamaños muestrales del grupo focal y del resto.
  ],
)

#pagebreak()

== Tabla A3. Mann-Whitney uno contra todos para características fisicoquímicas <table:anexo_mw3>

#figure(
  table(
    align: center,
    columns: 11,
    table.header[feature][n_true][n_false][median_true][median_false][W][p_value][cliff_delta][cliff_magnitude][p_adjusted][direction],
    [Sustitución conservativa],
    [12],
    [216],
    [0.001283050],
    [0.003111392],
    [346],
    [1.961420e-05],
    [-0.7330247],
    [large],
    [0.0002353704],
    [TRUE < FALSE],

    [Elimina glicina],
    [19],
    [209],
    [0.004851520],
    [0.002862394],
    [3078],
    [7.282936e-05],
    [0.5502392],
    [large],
    [0.0008011229],
    [TRUE > FALSE],

    [Cambio de clase],
    [177],
    [51],
    [0.003166020],
    [0.001979709],
    [5953],
    [5.261751e-04],
    [0.3189321],
    [small],
    [0.0052617509],
    [TRUE > FALSE],

    [Cambio de grupo hidrofóbico],
    [177],
    [51],
    [0.003166020],
    [0.001979709],
    [5953],
    [5.261751e-04],
    [0.3189321],
    [small],
    [0.0052617509],
    [TRUE > FALSE],

    [Cambio de carga],
    [120],
    [108],
    [0.003335714],
    [0.002595782],
    [8143],
    [8.288788e-04],
    [0.2566358],
    [small],
    [0.0066310305],
    [TRUE > FALSE],

    [Elimina prolina],
    [19],
    [209],
    [0.003967524],
    [0.002818704],
    [2747],
    [5.702012e-03],
    [0.3835306],
    [medium],
    [0.0399140825],
    [TRUE > FALSE],

    [Cambio aromático],
    [50],
    [178],
    [0.002623260],
    [0.003192931],
    [3468],
    [1.723777e-02],
    [-0.2206742],
    [small],
    [0.1034266436],
    [TRUE < FALSE],

    [Introduce glicina],
    [11],
    [217],
    [0.004028022],
    [0.002896070],
    [1682],
    [2.222647e-02],
    [0.4093004],
    [medium],
    [0.1111323682],
    [TRUE > FALSE],

    [Introduce prolina],
    [11],
    [217],
    [0.004059970],
    [0.002896070],
    [1634],
    [3.924812e-02],
    [0.3690825],
    [medium],
    [0.1569924853],
    [TRUE > FALSE],

    [Elimina cisteína],
    [19],
    [209],
    [0.003736556],
    [0.002862394],
    [2507],
    [5.841002e-02],
    [0.2626542],
    [small],
    [0.1752300745],
    [TRUE > FALSE],

    [Cambio de tamaño],
    [156],
    [72],
    [0.003149003],
    [0.002625674],
    [6380],
    [9.912218e-02],
    [0.1360399],
    [negligible],
    [0.1797282840],
    [TRUE > FALSE],

    [Introduce cisteína],
    [11],
    [217],
    [0.002531826],
    [0.003073454],
    [831],
    [8.986414e-02],
    [-0.3037285],
    [small],
    [0.1797282840],
    [TRUE < FALSE],
  ),
  caption: [
    Resultados completos de pruebas de Mann-Whitney uno contra todos para características fisicoquímicas.
    $n_"true"$/$n_"false"$: cantidad de mutaciones con y sin la característica.
  ],
)
