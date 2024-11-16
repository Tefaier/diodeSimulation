import math

import matplotlib.pyplot as plt
import matplotlib.colors as colors
import numpy as np
from scipy.integrate import odeint

from Simulation.simulationConstants import electron_mass, k, electron_charge

VoltFunction = lambda t: 1
functionToCheck = [
    lambda t: 1,
    lambda t: 0.2,
]
a_k_offset = 0.01
initialOffset = a_k_offset / 100000

def getData(voltFunction, timeLimit):
    def model(p, t):
        dp = -1 * k * electron_charge * voltFunction(t) / (electron_mass * (a_k_offset - (p[0] ** 2) * 0.5))
        return dp

    timeDots = np.linspace(0, timeLimit, 5000)
    solution = odeint(model, math.sqrt(initialOffset * 2), timeDots)
    solution = np.power(solution, 2) * 0.5
    solution = [x[0] for x in solution.tolist()]
    cutOffIndex = next((i for i, v in enumerate(solution) if v > a_k_offset - initialOffset), len(solution))
    print(cutOffIndex)
    solution = solution[:cutOffIndex]
    return (timeDots[:cutOffIndex], solution)

def simple():
    timeLimit = 7e-25
    fig = plt.figure(1, (9, 6))
    axes = plt.axes()
    axes.set_ylim([0, a_k_offset])
    axes.set_xlim([0, timeLimit])
    voltFunction = lambda t: 2
    time, pos = getData(voltFunction, timeLimit)
    plt.plot(time, pos)
    plt.show()

def simpleVoltage():
    timeLimit = 7e-24
    fig = plt.figure(1, (9, 6))
    axes = plt.axes()
    axes.set_ylim([0, a_k_offset])
    axes.set_xlim([0, timeLimit])
    for val in np.linspace(0.2, 5, 50):
        voltFunction = lambda t: val
        time, pos = getData(voltFunction, timeLimit)
        time = [t * 1e9 for t in time]
        pos = [p * 1e3 for p in pos]
        plt.plot(time, pos, color=colors.hsv_to_rgb([(val - 0.2) / 10, 1, 1]))
    plt.plot([], [], color=colors.hsv_to_rgb([1., 1., 1.]), label="0,2 V")
    plt.plot([], [], color=colors.hsv_to_rgb([(5 - 0.2) / 10, 1., 1.]), label="5 V")
    plt.legend()
    plt.title("Зависимость движения электрона при разном анодном напряжении")
    plt.xlabel("Время с начала полета, нс")
    plt.ylabel("Расстояние от катода, мм")
    plt.show()

def oscillation():
    timeLimit = 4e-24
    fig = plt.figure(1, (9, 6))
    axes = plt.axes()
    axes.set_ylim([0, a_k_offset * 1e3])
    axes.set_xlim([0, timeLimit * 1e9])
    for val in np.linspace(0.1, 3, 20):
        voltFunction = lambda t: 2 * math.sin(t * (1 / 2.966655551850617e-25) * val)
        time, pos = getData(voltFunction, timeLimit)
        time = [t * 1e9 for t in time]
        pos = [p * 1e3 for p in pos]
        plt.plot(time, pos, color=colors.hsv_to_rgb([(val - 0.1) / 5, 1, 1]))

    # fly time 2.966655551850617e-25
    voltFunction = lambda t: 2
    time, pos = getData(voltFunction, timeLimit)
    time = [t * 1e9 for t in time]
    pos = [p * 1e3 for p in pos]
    plt.plot(time, pos, color=[0, 0, 0], label="постоянный ток")
    plt.plot([], [], label="угловая скорость w = 1 / время пролета при постоянном токе")
    plt.plot([], [], color=colors.hsv_to_rgb([0., 1., 1.]), label="w * 0,1")
    plt.plot([], [], color=colors.hsv_to_rgb([(3 - 0.1) / 5, 1., 1.]), label="w * 3")
    plt.legend()
    plt.title("Зависимость движения электрона при разной частоте переменного тока")
    plt.xlabel("Время с начала полета, нс")
    plt.ylabel("Расстояние от катода, мм")
    plt.show()

if __name__ == '__main__':
    oscillation()
