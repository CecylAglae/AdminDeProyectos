from fastaReader import fastaReader
import random
import numpy
import numpy as np
import copy
from evaluadorBlosum import evaluadorBlosum

class bacteria():
    
    def __init__(self, path):
        self.matrix = fastaReader(path)  # Carga las secuencias desde un archivo
        self.blosumScore = 0
        self.fitness = 0
        self.interaction = 0
        self.NFE = 0

    def showGenome(self):
        for seq in self.matrix.seqs:
            print(seq)  # Muestra todas las secuencias

    def clonar(self, path):
        newBacteria = bacteria(path)  # Crea una nueva instancia de bacteria
        newBacteria.matrix.seqs = np.array(copy.deepcopy(self.matrix.seqs))  # Copia las secuencias
        return newBacteria

    def tumboNado(self, numGaps):
        self.cuadra()  # Asegura que todas las secuencias tengan la misma longitud
        matrixCopy = copy.deepcopy(self.matrix.seqs)  # Crea una copia de las secuencias
        matrixCopy = matrixCopy.tolist()  # Convierte la matriz a lista
        gapRandomNumber = random.randint(0, numGaps)  # Selecciona un número aleatorio de gaps
        for i in range(gapRandomNumber):
            seqnum = random.randint(0, len(matrixCopy)-1)  # Selecciona una secuencia aleatoria
            pos = random.randint(0, len(matrixCopy[0])-1)  # Selecciona una posición aleatoria
            part1 = matrixCopy[seqnum][:pos]  # Parte antes de la posición
            part2 = matrixCopy[seqnum][pos:]  # Parte después de la posición
            temp = "-".join([part1, part2])  # Inserta un gap
            matrixCopy[seqnum] = temp  # Actualiza la secuencia con el gap
        matrixCopy = np.array(matrixCopy)  # Convierte de nuevo a matriz
        self.matrix.seqs = matrixCopy  # Actualiza la matriz de secuencias
        self.cuadra()  # Rellena las secuencias cortas
        self.limpiaColumnas()  # Elimina columnas con solo gaps

    def cuadra(self):
        """Rellena con gaps las secuencias más cortas"""
        seq = self.matrix.seqs
        maxLen = len(max(seq, key=len))  # Encuentra la longitud máxima
        for i in range(len(seq)):
            if len(seq[i]) < maxLen:
                seq[i] = seq[i] + "-"*(maxLen-len(seq[i]))  # Añade gaps para igualar longitud
        self.matrix.seqs = np.array(seq)  # Actualiza la matriz de secuencias

    def gapColumn(self, col):
        """Verifica si la columna tiene gaps en todas las filas"""
        for i in range(len(self.matrix.seqs)):
            if self.matrix.seqs[i][col] != "-":
                return False  # Devuelve False si encuentra un elemento que no es gap
        return True  # Devuelve True si todos son gaps

    def limpiaColumnas(self):
        """Elimina las columnas que tienen gaps en todas las filas"""
        i = 0
        while i < len(self.matrix.seqs[0]):
            if self.gapColumn(i):  # Comprueba si la columna tiene solo gaps
                self.deleteCulmn(i)  # Elimina la columna
            else:
                i += 1  # Avanza a la siguiente columna

    def deleteCulmn(self, pos):
        """Elimina un elemento específico en cada secuencia"""
        for i in range(len(self.matrix.seqs)):
            self.matrix.seqs[i] = self.matrix.seqs[i][:pos] + self.matrix.seqs[i][pos+1:]  # Actualiza la secuencia

    def getColumn(self, col):
        """Obtiene los elementos de una columna específica"""
        column = []
        for i in range(len(self.matrix.seqs)):
            column.append(self.matrix.seqs[i][col])  # Añade el elemento de la columna
        return column  # Devuelve la lista de elementos de la columna

    def calculate_penalty(self):
        gap_penalty = sum(seq.count("-") for seq in self.matrix.seqs) * 0.01  # Penaliza por la cantidad de gaps

        scores = []  # Lista para almacenar puntajes
        evaluador = evaluadorBlosum()  # Crea un evaluador de Blosum
        for i in range(len(self.matrix.seqs[0])):
            column = self.getColumn(i)  # Obtiene la columna actual
            column = [x for x in column if x != "-"]  # Filtra los gaps
            pares = self.obtener_pares_unicos(column)  # Obtiene pares únicos de la columna
            for par in pares:
                score = evaluador.getScore(par[0], par[1])  # Calcula el puntaje para cada par
                scores.append(score)  # Añade el puntaje a la lista

        score_variability_penalty = np.std(scores) * 0.3 if scores else 0  # Penaliza por variabilidad en puntajes
        total_penalty = gap_penalty + score_variability_penalty  # Calcula la penalización total
        return total_penalty  # Devuelve la penalización total

    def autoEvalua(self):
        evaluador = evaluadorBlosum()  # Crea un evaluador de Blosum
        score = 0  # Inicializa el puntaje

        # Calcula el puntaje de Blosum
        for i in range(len(self.matrix.seqs[0])):
            column = self.getColumn(i)  # Obtiene la columna actual
            gapCount = column.count("-")  # Cuenta los gaps en la columna
            column = [x for x in column if x != "-"]  # Filtra los gaps
            pares = self.obtener_pares_unicos(column)  # Obtiene pares únicos de la columna
            for par in pares:
                score += evaluador.getScore(par[0], par[1])  # Suma el puntaje de cada par
            score -= gapCount * 2  # Penaliza por los gaps

        penalty = self.calculate_penalty()  # Calcula la penalización total

        # Calcula el fitness final
        self.blosumScore = score - penalty  # Ajusta el puntaje con la penalización
        self.fitness = self.blosumScore  # Actualiza el fitness
        self.NFE += 1  # Incrementa el contador de evaluaciones

    def obtener_pares_unicos(self, columna):
        """Obtiene pares únicos de una lista"""
        pares_unicos = set()  # Usa un conjunto para evitar duplicados
        for i in range(len(columna)):
            for j in range(i+1, len(columna)):
                par = tuple(sorted([columna[i], columna[j]]))  # Ordena los elementos
                pares_unicos.add(par)  # Añade el par al conjunto
        return list(pares_unicos)  # Devuelve la lista de pares únicos
