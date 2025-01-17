# PCA_reads_kmers
How to use Principal Component Analysis in Computational Biology

Source: https://medium.com/computational-biology/pca-for-bioinformaticians-4cf2cc7e792f

Principal Component Analysis (PCA) is one of the most popular dimension reduction techniques in Machine Learning and Statistics. PCA is very useful in order to visualize data in human-readable formats, such as 2D plots for better understanding of high-dimensional data. We can obtain the same set of advantages in the domain of Bioinformatics by applying PCA on our genomic/proteomic datasets. In this article, I will explain how one might use PCA to differentiate between two species using the species’ oligonucleotide frequency vectors.

# Simulating Reads

In order to demonstrate the use of PCA, first, we shall download two species. In this article I will be using the two species;

Descargar las siguientes secuencias del NCBI

Staphylococcus aureus (Nuccore ID NC_007795.1)
Pseudomonas aeruginosa (Nuccore ID NC_002516.2)

Instalar simlord y dependencias dentro de un ambiente de conda segun esta pagina: https://bitbucket.org/genomeinformatics/simlord/src/master/
y ejecutar en la terminal las siguientes lineas
```
#!/bin/bash
simlord --no-sam -rr SA.fasta -n 100 -mr 1500 SA
simlord --no-sam -rr PA.fasta -n 100 -mr 1500 PA
```

# Vectorize sequences
source: https://medium.com/computational-biology/vectorization-of-dna-sequences-2acac972b2ee


I compile the above program using the following command;

```
#!/bin/bash
sudo apt install g++
g++ vectorize.cpp -o vectorize -fopenmp --std=c++17
```

I run the following commands to get the vectors for each species’ fastq reads.

```
#!/bin/bash
./vectorize ./SA.fastq 3 8 100 > ./SA.txt
./vectorize ./PA.fastq 3 8 100 > ./PA.txt
```
Vamos a descomponer cada elemento:

```./vectorize```: Esto ejecuta el programa compilado (el ejecutable vectorize).
```./SA.fastq``` y ```./PA.fastq```: Estos son los archivos FASTQ que contienen las lecturas genómicas para dos especies diferentes, posiblemente Staphylococcus aureus (SA) y Pseudomonas aeruginosa (PA). El programa parece procesar estos archivos para generar vectores basados en las secuencias de ADN contenidas en ellos.
```3```: Es un parámetro pasado al programa que probablemente corresponde a la longitud de los k-mers que se van a generar. En este caso, parece que el programa generará k-mers de longitud 3 (trinucleótidos).
```8```: Otro parámetro pasado al programa, que podría estar relacionado con el número de procesos paralelos (hilos) que se deben usar con OpenMP, lo que indicaría que el programa podría estar aprovechando 8 núcleos del procesador para realizar la tarea de manera más eficiente.
```100```: Podría ser un parámetro que controla cuántas secuencias o lecturas del archivo FASTQ se deben procesar, o quizás el tamaño de la ventana en la que se aplican los k-mers.
```> ./SA.txt``` y ```> ./PA.txt```: Estos redirigen la salida del programa a archivos de texto


# Computing PCA and Plotting on 2D Plane (python de preferencia en Jupyter notebook)

Now that we have vectorized the sequences from long reads, we can use Python’s Sci-Kit learn library to obtain the PCA decomposition. I will use the Seaborn library for Python to visualize them on the 2D plane.

```
# Imports
import numpy as np
import seaborn as sns
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt# Reading data
sa = np.loadtxt("SA.txt")
pa = np.loadtxt("PA.txt")# Preparing labels
data = np.append(pa, sa, axis=0)
labels  = ["Staphylococcus aureus" for x in range(100)]
labels += ["Pseudomonas aeruginosa" for x in range(100)]
labels = np.array(labels)# Run PCA 
pca = PCA(n_components=2)
data_2d = pca.fit_transform(data)# Plotting using Seaborn
fig = plt.figure(figsize=(10, 10))
sns.scatterplot(data_2d[:,0], data_2d[:,1], hue=labels)
plt.xlabel("PCA 1")
plt.ylabel("PCA 2")
plt.savefig("PCA.png")

The above code will plot us the following diagram.
```

# PCA Versión corregida

```
# Load data
sa = np.loadtxt("SA.txt")
pa = np.loadtxt("PA.txt")

# Prepare data and labels
data = np.append(pa, sa, axis=0)
labels = ["Staphylococcus aureus" for x in range(100)] + ["Pseudomonas aeruginosa" for x in range(100)]
labels = np.array(labels)

# Run PCA
pca = PCA(n_components=2)
data_2d = pca.fit_transform(data)

# Create a DataFrame for plotting
df = pd.DataFrame(data_2d, columns=["PCA 1", "PCA 2"])
df["Label"] = labels

# Plotting using Seaborn
fig = plt.figure(figsize=(10, 10))
sns.scatterplot(x="PCA 1", y="PCA 2", hue="Label", data=df)
plt.xlabel("PCA 1")
plt.ylabel("PCA 2")
plt.title("PCA of Staphylococcus aureus and Pseudomonas aeruginosa")
plt.savefig("PCA.png")
plt.show()

```


