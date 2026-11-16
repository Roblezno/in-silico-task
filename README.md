# in-silico-task

Reanálisis de RNA-seq del obelisco-S.s en *Streptococcus sanguinis* SK36 y benchmark de
mapeadores para su cuantificación. Trabajo final de la asignatura *Estudios in silico en
Biomedicina* (máster en Bioinformática).

## Datos

- **PRJNA1129866**: RNA-seq de Zheludev et al. (2024), subcepas ObP1 y ObN1.
- **PRJNA937727**: RNA-seq de SK36 en monocultivo (diseño oxígeno × azúcar).

Los datos no se versionan: los scripts los descargan en `data/`.

## Entorno

```bash
mamba env create -f environment.yml --strict-channel-priority
mamba activate in-silico-task
```

## Estructura

```
config/    tablas de muestras y parámetros
scripts/   pasos del análisis, numerados en orden de ejecución
report/    informe en Quarto
data/      datos descargados y referencias (no versionado)
results/   resultados intermedios y finales (no versionado)
```

## Estado

En desarrollo.
