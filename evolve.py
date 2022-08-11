from gapy.genome import Genome
from gapy.geneticAlgoritm import GeneticAlgorithm
import random
import polycubes.solve.py.utils as utils
import libpolycubes as lp
import json

class PolycubeGenome(Genome):
    def __init__(self, rule, mutationRate=0.9, mutationWeights=[
        1, # Change orientation
        5, # Color replacement
        1, # Remove species
        1, # Update color
    ]):
        self.rule = utils.parseDecRule(utils.ruleToDec(rule))
        self.mutationRate = mutationRate
        self.mutationWeights = mutationWeights
        super().__init__()

    def clone(self):
        other = PolycubeGenome(utils.parseDecRule(utils.ruleToDec(self.rule)))
        return other

    def mutateRemoveSpecies(self):
        species = random.choice(self.rule)
        self.rule = [s for s in self.rule if s is not species]

    def mutateChangeOrientation(self):
        species = random.choice(self.rule)
        face = random.choice(species)
        face['orientation'] = random.choice(range(4))

    def mutateUpdateColor(self):
        maxColor = countColors(self.rule)
        species = random.choice(self.rule)
        face = random.choice(species)
        face['color'] = round(random.uniform(-maxColor, maxColor))

    def mutateColorReplace(self):
        colors = set(abs(face['color']) for species in self.rule for face in species)
        oldColor = random.choice(list(colors))
        newColor = random.choice(range(oldColor + 2))
        for species in self.rule:
            for face in species:
                if face['color'] == oldColor:
                    face['color'] = newColor
                elif face['color'] == -oldColor:
                    face['color'] = -newColor

    def mutate(self):
        if random.random() < self.mutationRate:
            mutations = random.choices([
                self.mutateChangeOrientation,
                self.mutateColorReplace,
                self.mutateRemoveSpecies,
                self.mutateUpdateColor
            ], weights=self.mutationWeights, k=max(self.mutationRate, 1))
            for m in mutations:
                m()
        self.rule = utils.simplifyRule(self.rule)

    def __str__(self) -> str:
        return utils.ruleToDec(self.rule)

def onGenerationStep(generation, maxFitness, bestGenome):
    nColors = countColors(bestGenome.rule)
    nSpecies = len(bestGenome.rule)
    print("\nGeneration {}\n  Max fitness: {}\n ({} colors and {} species)  Best genome: {}\n".format(
        generation, maxFitness, nColors, nSpecies, bestGenome))


populationSize = 25
nGenerations = 100

def countColors(rule):
    return max(abs(face['color']) for species in rule for face in species)


def countNonZero(rule):
    return sum(1 for species in rule for face in species if face['color'] != 0)

def run(solveSpecPath):
    with open(solveSpecPath, 'r') as f:
        data = f.read()
    solveSpec = json.loads(data)
    top = solveSpec['bindings']
    maxRule = utils.getFullyAdressableRule(top)
    maxRule = utils.simplifyRule(maxRule) # Why is this neccesary?
    decMaxRule = utils.ruleToDec(maxRule)
    maxComplexity = countNonZero(maxRule)
    print("Fully adressable rule: "+decMaxRule)
    coords = utils.calcCoordsFromTop(top)[0]

    def fitnessFunc(genome):
        assembleRatio = lp.assembleRatio(
            coords,
            utils.ruleToDec(genome.rule),
            isHexString=False,
            nTries=250,
            assemblyMode='stochastic'
        )
        if assembleRatio < 1:
            print(assembleRatio, end=' ', flush=True)
            return assembleRatio
        else:
            #nColors = max(face['color'] for species in genome.rule for face in species)
            #nCubeTypes = len(genome.rule)
            #return 1 / (nColors * nCubeTypes)
            #nFewerSpecies = len(maxRule) - len(genome.rule)
            nFewerNonZero = maxComplexity - countNonZero(genome.rule)
            print(nFewerNonZero+1, end=' ', flush=True)
            #return nFewerSpecies + nFewerColors + 1
            return nFewerNonZero + 1


    # Initialize population
    population = [PolycubeGenome(maxRule) for _ in range(populationSize)]

    # Introduce some genetic diversity
    for i in population:
        if random.random() < 0.9:
            i.mutate()

    # Initialise genetic algorithm
    evolver = GeneticAlgorithm(population, fitnessFunc)

    # Evolve for a given number of generations
    evolver.run(nGenerations, onGenerationStep, nProcesses=5)

if __name__ == "__main__":
    solveSpecPath = 'polycubes/solve/shapes/scaling/cube5.json'
    run(solveSpecPath)