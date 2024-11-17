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
initialOffset = a_k_offset / 1000
area = 1e-4

def getData(voltFunction, timeLimit):
    def model(p, t):
        dp = electron_mass * voltFunction(t) / (electron_mass * a_k_offset)
        #  dp = -2 * electron_charge * voltFunction(t) * area / (electron_mass * 4 * math.pi * a_k_offset * ((a_k_offset - (p[0] ** 2) * 0.5) ** 2))
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
    timeLimit = 7e-4
    fig = plt.figure(1, (9, 6))
    axes = plt.axes()
    axes.set_ylim([0, a_k_offset])
    axes.set_xlim([0, timeLimit])
    voltFunction = lambda t: 5
    time, pos = getData(voltFunction, timeLimit)
    plt.plot(time, pos)
    plt.show()

def simpleVoltage():
    timeLimit = 7e-4
    fig = plt.figure(1, (9, 6))
    axes = plt.axes()
    axes.set_ylim([0, a_k_offset * 1e3])
    axes.set_xlim([0, timeLimit * 1e6])
    for val in np.linspace(0.2, 10, 50):
        voltFunction = lambda t: val
        time, pos = getData(voltFunction, timeLimit)
        time = [t * 1e6 for t in time]
        pos = [p * 1e3 for p in pos]
        plt.plot(time, pos, color=colors.hsv_to_rgb([(val - 0.2) / 15, 1, 1]))
    plt.plot([], [], color=colors.hsv_to_rgb([1., 1., 1.]), label="0,2 V")
    plt.plot([], [], color=colors.hsv_to_rgb([(10 - 0.2) / 15, 1., 1.]), label="10 V")
    plt.legend()
    plt.title("Зависимость движения электрона при разном анодном напряжении")
    plt.xlabel("Время с начала полета, us")
    plt.ylabel("Расстояние от катода, мм")
    plt.show()

def binarySearch(function, precision: float, left: float, right: float):
    while right - left > precision:
        mid = (right + left) / 2
        if function(mid):
            left = mid
        else:
            right = mid
    return right

def oscillation():
    timeLimit = 2e-3
    fig = plt.figure(1, (9, 6))
    axes = plt.axes()
    axes.set_ylim([0, a_k_offset * 1e3])
    axes.set_xlim([0, timeLimit * 1e6])
    for val in np.linspace(0.1, 3, 20):
        voltFunction = lambda t: 5 * math.sin(t * (1 / 0.0002737547509501901) * val)
        time, pos = getData(voltFunction, timeLimit)
        time = [t * 1e6 for t in time]
        pos = [p * 1e3 for p in pos]
        plt.plot(time, pos, color=colors.hsv_to_rgb([(val - 0.1) / 5, 1, 1]))

    # fly time 0.0002737547509501901
    voltFunction = lambda t: 5
    time, pos = getData(voltFunction, timeLimit)
    time = [t * 1e6 for t in time]
    pos = [p * 1e3 for p in pos]
    plt.plot(time, pos, color=[0, 0, 0], label="постоянный ток")
    plt.plot([], [], label="угловая скорость w = 1 / время пролета при постоянном токе")
    plt.plot([], [], color=colors.hsv_to_rgb([0., 1., 1.]), label="w * 0,1")
    plt.plot([], [], color=colors.hsv_to_rgb([(3 - 0.1) / 5, 1., 1.]), label="w * 3")
    plt.legend()
    plt.title("Зависимость движения электрона при разной частоте переменного тока")
    plt.xlabel("Время с начала полета, us")
    plt.ylabel("Расстояние от катода, мм")
    plt.show()

def binarySearchForOscillation():
    timeLimit = 2e-3
    def boolFunc(val):
        voltFunction = lambda t: 5 * math.sin(t * (1 / 0.0002737547509501901) * val)
        time, pos = getData(voltFunction, timeLimit)
        return len(time) < 5000
    result = binarySearch(boolFunc, 0.0000000001, 0.1, 3)
    print(result)
    # result 1.9999828 appr

if __name__ == '__main__':
    binarySearchForOscillation()
