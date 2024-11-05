import math
import random
import numpy as np
from bacteria import bacteria

class chemiotaxis():
    def __init__(self):
        self.parcialNFE = 0  # Contador de evaluaciones parciales
    
    def compute_cell_interaction(self, bacteria, poblacion, d, w):
        """Calcula la interacción total con otras bacterias en función de su puntaje de Blosum."""
        total = 0.0
        for other in poblacion:
            diff = (bacteria.blosumScore - other.blosumScore)
            total += d * math.exp(w * diff)  # Suma la interacción
        return total

    def attract_repel(self, bacteria, poblacion, d_attr, w_attr, h_rep, w_rep):
        """Calcula la atracción y repulsión de la bacteria en la población."""
        attract = self.compute_cell_interaction(bacteria, poblacion, -d_attr, -w_attr)  # Atracción
        repel = self.compute_cell_interaction(bacteria, poblacion, h_rep, -w_rep)  # Repulsión
        return attract + repel  # Resultado de la interacción
    
    def chemio(self, bacteria, poblacion, d_attr, w_attr, h_rep, w_rep):
        """Calcula la interacción y fitness de la bacteria."""
        bacteria.interaction = self.attract_repel(bacteria, poblacion, d_attr, w_attr, h_rep, w_rep)
        bacteria.fitness = bacteria.blosumScore + bacteria.interaction  # Actualiza fitness

    def doChemioTaxis(self, poblacion, d_attr, w_attr, h_rep, w_rep):
        """Aplica quimiotaxis a toda la población de bacterias."""
        for bacteria in poblacion:
            bacteria.autoEvalua()  # Evalúa la bacteria
            self.chemio(bacteria, poblacion, d_attr, w_attr, h_rep, w_rep)  # Calcula quimiotaxis

            self.parcialNFE += bacteria.NFE  # Acumula NFE
            bacteria.NFE = 0  # Reinicia NFE de la bacteria

        # Imprime el fitness final de cada bacteria
        for bacteria in poblacion:
            print(f"Fitness final de la bacteria: {bacteria.fitness}")

    def eliminarClonar(self, path, poblacion):
        """Elimina bacterias de menor fitness y clona las restantes."""
        poblacion.sort(key=lambda x: x.fitness)  # Ordena por fitness
        num_eliminar = int(len(poblacion) * 0.7)  # Calcula cuántas eliminar
        del poblacion[:num_eliminar]  # Elimina el 70% de menor fitness
        
        # Clona las bacterias restantes y las añade a la población
        clones = self.clonacion(path, poblacion)
        poblacion.extend(clones)
    
    def clonacion(self, path, poblacion):
        """Crea clones de las bacterias de la población."""
        poblacionClones = []
        best = max(poblacion, key=lambda x: x.fitness)  # Mejor bacteria
        fitness_stddev = np.std([b.fitness for b in poblacion])  # Desviación estándar del fitness
        
        for bacteria in poblacion:
            newBacteria = bacteria.clonar(path)  # Clona la bacteria
            
            # Calcula mutación si la desviación estándar es mayor que cero
            if fitness_stddev > 0:
                mutacion = int((best.fitness - bacteria.fitness) / (10 * fitness_stddev))
            else:
                mutacion = 0  # Si no hay desviación estándar
            
            newBacteria.tumboNado(mutacion)  # Aplica mutación
            newBacteria.autoEvalua()  # Evalúa el clon
            poblacionClones.append(newBacteria)  # Añade el clon a la lista
        return poblacionClones

    def randomBacteria(self, path):
        """Crea una bacteria aleatoria."""
        bact = bacteria(path)  # Crea una nueva bacteria
        bact.tumboNado(random.randint(1, 10))  # Aplica un nado aleatorio
        return bact
   
    def insertRamdomBacterias(self, path, num, poblacion):
        """Inserta bacterias aleatorias en la población."""
        for i in range(num):
            poblacion.append(self.randomBacteria(path))  # Añade una bacteria aleatoria
            # Elimina la bacteria con menor fitness
            poblacion.sort(key=lambda x: x.fitness)  # Ordena por fitness
            del poblacion[0]  # Elimina la de menor fitness
