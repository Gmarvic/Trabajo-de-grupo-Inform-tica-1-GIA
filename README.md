
[Web - Optimización Comercial](https://marcoaias.github.io/Optimizacion-Comercial/)

# Optimizacion-Comercial
GIA 1º Informática - Trabajo de Equipo 20/21


La compañía de trenes de Alta Velocidad, Correcaminos Mic-Mic, ofrece servicios punto a punto
entre Madrid y 6 ciudades peninsulares (1 convoy diario a cada destino).

Vosotros formáis parte del equipo de Optimización Comercial y vuestra responsabilidad es
proporcionar las recomendaciones para comercializar los productos que la compañía pone a la
venta en su web y en los sistemas de Distribución a los que acceden las agencias de viaje, con el
fin de poder maximizar los ingresos en cada convoy.

(Nota: Recomendar como deben comercializarse los productos, es fijar el número máximo de
asientos que se puede poner a la venta en un convoy para cada producto ofertado por la
compañía)

Los 6 destinos servidos son (vamos a hacer un par de supuestos: 1- que la red de Alta Velocidad
ya está muy extendida llegando a muchos puntos de España y 2: que el producto ofertado es
exclusivamente punto a punto, es decir va de origen a destino):
Barcelona
Sevilla
Valencia
Malaga
Bilbao
Santiago de Compostela

La compañía dispone de tres máquinas de Alta Velocidad, que emplea diariamente para llegar a
cada uno de los destinos (3 servicios de ida y vuelta en horario de mañana y 3 en horario de
tarde), configurando el convoy con 2, 3 o 4 vagones (en función de la demanda), en una sola
cabina, cabina Turista, con 80 plazas por vagón.

Es asumido por el Departamento de Comercial que la demanda está perfectamente segmentada
en 5 diferentes perfiles o clases de clientes, siendo cada perfil o clase independientes del resto.
Los precios son comunes para todos los destinos desde Madrid y para cada uno de estos
segmentos de potenciales clientes, la compañía ofrece un producto ajustado a sus necesidades,
y una clase de reserva para que el cliente pueda realizar una reserva de asiento con las
características del servicio que quiere obtener.

Las características de los productos, vienen resumidas en la siguiente tabla:
Clase de reserva Tarifa Cambios permitidos
Sala
Vip
Fast
track Elección asiento
A 180 Si Si Si Si
B 130 Solo 1 (penalización 25€) No Si Si
C 100 Solo 1 (penalización 60€) No No Si
D 80 No No No No
E 40 No No No No


Objetivo del trabajo:
Se ha dividido el curso en 6 equipos y se ha asignado a cada uno de los 6 equipos un destino,
eligiendo al azar un team leader para cada equipo.

Para cada uno de los destinos, el departamento de forecasting ha proporcionado una demanda
media esperada y un error o desviación esperado para la operación del día 12 de Julio de 2021.
(Nota: En algunos casos esta demanda y error serán las soluciones de un sistema de ecuaciones
no lineales que relacionan consumo y precio)

Se pide a cada uno de los equipos, que proporcione la comercialización óptima que maximizaría
los ingresos del convoy que les haya correspondido para el 12 de Julio de 2021.

La técnica que se utilizará para la optimización será una aproximación al algoritmo EMSR B
desarrollado por Peter Bellobaba (PRINCIPAL RESEARCH SCIENTIST, Department of Aeronautics
and Astronautics, MIT) y para su desarrollo utilizaremos los materiales utilizados a lo largo del
curso, de los que destacamos como principales:
- Lenguaje Fortran:
o Bucles, sentencias condicionales.
o Sentencias de lectura, escritura en archivos.
o Estructuras vectores, matrices.
o Programación Modular.
- Codeblocks
- Algoritmos de resolución de Sistemas de Ecuaciones Lineales y no Lineales.
- Método Trapecio para cálculo de integrales simples.
- Algoritmos clasificación / selección.
