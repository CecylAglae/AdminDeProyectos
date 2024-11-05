from bacteria import bacteria
from chemiotaxis import chemiotaxis
import numpy

poblacion = []  # Lista para almacenar la población de bacterias
path = r"C:\Users\Cecyl\Documents\Uni\AdministracioinDeProyectos\BFOA-Mejorado\multiFasta.fasta"
temp_bacteria = bacteria(path)  # Instancia temporal de bacteria
tumbo = 1  # Parámetro para el movimiento de nado
numeroDeBacterias = 8  # Número total de bacterias en la población
numRandomBacteria = 2  # Número de bacterias aleatorias a insertar
iteraciones = 5  # Número de iteraciones del proceso
nado = 4  # Parámetro para el nado
chemio = chemiotaxis()  # Instancia del objeto de quimiotaxis
veryBest = bacteria(path)  # Mejores bacterias encontradas
tempBacteria = bacteria(path)  # Instancia para validación
original = bacteria(path)  # Instancia de las bacterias originales
globalNFE = 0  # Contador global de evaluaciones
dAttr = 0.15  # Atributo de distancia para la quimiotaxis
wAttr = 0.15  # Peso para el atributo de quimiotaxis
hRep = dAttr  # Repulsión basada en distancia
wRep = 0.08  # Peso para la repulsión

def clonaBest(veryBest, best):
    """Clona la mejor bacteria encontrada"""
    if veryBest is None:
        veryBest = bacteria(path)  # Inicializa veryBest si es None
    
    veryBest.matrix.seqs = numpy.array(best.matrix.seqs)  # Copia las secuencias
    veryBest.blosumScore = best.blosumScore  # Copia el puntaje de Blosum
    veryBest.fitness = best.fitness  # Copia el fitness
    veryBest.interaction = best.interaction  # Copia la interacción
    return veryBest  # Devuelve la mejor bacteria

def validaSecuencias(path, veryBest):
    """Valida que las secuencias originales coincidan con las de veryBest"""
    tempBacteria.matrix.seqs = numpy.array(veryBest.matrix.seqs)  # Clona las secuencias
    for i in range(len(tempBacteria.matrix.seqs)):
        tempBacteria.matrix.seqs[i] = tempBacteria.matrix.seqs[i].replace("-", "")  # Elimina gaps

    # Compara las secuencias originales con las de tempBacteria
    for i in range(len(tempBacteria.matrix.seqs)):
        if tempBacteria.matrix.seqs[i] != original.matrix.seqs[i]:
            print("*****************Secuencias no coinciden********************")
            return  # Sale si hay diferencias
      
# Inicializa la población de bacterias
for i in range(numeroDeBacterias):
    poblacion.append(bacteria(path))  # Añade una nueva bacteria a la población

# Inicialización de veryBest antes del bucle
veryBest = bacteria(path)  # Crea la primera bacteria

for _ in range(iteraciones):  # Iteraciones del proceso
    for bacteria in poblacion:
        bacteria.tumboNado(tumbo)  # Introduce cambios en las bacterias
        bacteria.autoEvalua()  # Evalúa el rendimiento de cada bacteria
    
    chemio.doChemioTaxis(poblacion, dAttr, wAttr, hRep, wRep)  # Aplica quimiotaxis a la población
    globalNFE += chemio.parcialNFE  # Actualiza el contador global de evaluaciones
    
    best = max(poblacion, key=lambda x: x.fitness)  # Encuentra la mejor bacteria en la población
    if best.fitness > veryBest.fitness:  # Compara fitness con veryBest
        veryBest = clonaBest(veryBest, best)  # Clona la mejor bacteria
    
    print("Interacción: ", veryBest.interaction, " NFE:", globalNFE)  # Imprime interacción y NFE
    
    chemio.eliminarClonar(path, poblacion)  # Elimina bacterias clonadas
    chemio.insertRamdomBacterias(path, numRandomBacteria, poblacion)  # Inserta bacterias aleatorias
    print("poblacion: ", len(poblacion))  # Imprime la cantidad actual de bacterias

veryBest.showGenome()  # Muestra el genoma de la mejor bacteria
validaSecuencias(path, veryBest)  # Valida las secuencias de veryBest
