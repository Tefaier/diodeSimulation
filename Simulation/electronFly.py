import math

import matplotlib.pyplot as plt
import matplotlib.colors as colors
import numpy as np
from scipy.integrate import odeint, cumulative_trapezoid

from Simulation.simulationConstants import electron_mass, k, electron_charge

VoltFunction = lambda t: 1
functionToCheck = [
    lambda t: 1,
    lambda t: 0.2,
]
a_k_offset = 0.01
initialOffset = a_k_offset / 1000

def getData(voltFunction, timeLimit):
    def model(v_pos, t):
        dv = -1 * electron_charge * voltFunction(t) / (electron_mass * a_k_offset)
        dp = v_pos[0]
        #  dv = -2 * electron_charge * voltFunction(t) * area / (electron_mass * 4 * math.pi * a_k_offset * ((a_k_offset - (p[0] ** 2) * 0.5) ** 2))
        return [dv, dp]

    timeDots = np.linspace(0, timeLimit, 5000)
    solution = odeint(model, [0, initialOffset], timeDots)
    speed = solution[:, 0]
    position = solution[:, 1]
    cutOffIndex = next((i for i, v in enumerate(position) if v > a_k_offset), len(position))
    print(cutOffIndex)
    return (timeDots[:cutOffIndex], position[:cutOffIndex], speed[:cutOffIndex])

def getDataSpeed(voltFunction, timeLimit):
    def model(v, t):
        dv = -1 * electron_charge * voltFunction(t) / (electron_mass * a_k_offset)
        #  dp = -2 * electron_charge * voltFunction(t) * area / (electron_mass * 4 * math.pi * a_k_offset * ((a_k_offset - (p[0] ** 2) * 0.5) ** 2))
        return dv

    timeDots = np.linspace(0, timeLimit, 5000)
    solution = odeint(model, 0, timeDots)
    solution = [s[0] for s in solution]
    return (timeDots, solution)

def simple():
    timeLimit = 2e-8
    fig = plt.figure(1, (9, 6))
    axes = plt.axes()
    axes.set_ylim([0, a_k_offset])
    axes.set_xlim([0, timeLimit])
    voltFunction = lambda t: 5
    time, pos, _ = getData(voltFunction, timeLimit)
    plt.plot(time, pos)
    plt.show()

def simpleVoltage():
    timeLimit = 2e-8
    fig = plt.figure(1, (9, 6))
    axes = plt.axes()
    axes.set_ylim([0, a_k_offset * 1e3])
    axes.set_xlim([0, timeLimit * 1e9])
    for val in np.linspace(0.2, 20, 50):
        voltFunction = lambda t: val
        time, pos, _ = getData(voltFunction, timeLimit)
        time = [t * 1e9 for t in time]
        pos = [p * 1e3 for p in pos]
        plt.plot(time, pos, color=colors.hsv_to_rgb([(val - 0.2) / 30, 1, 1]))
    plt.plot([], [], color=colors.hsv_to_rgb([1., 1., 1.]), label="0,2 V")
    plt.plot([], [], color=colors.hsv_to_rgb([(20 - 0.2) / 30, 1., 1.]), label="20 V")
    plt.legend()
    plt.title("Зависимость движения электрона при разном анодном напряжении")
    plt.xlabel("Время с начала полета, ns")
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
    timeLimit = 8e-8
    fig = plt.figure(1, (9, 6))
    axes = plt.axes()
    axes.set_ylim([0, a_k_offset * 1e3])
    axes.set_xlim([0, timeLimit * 1e9])
    for val in np.linspace(0.1, 5, 20):
        voltFunction = lambda t: 5 * math.cos(t * (1 / 1.507101420284057e-08) * val)
        time, pos, _ = getData(voltFunction, timeLimit)
        pos = [p * 1e3 for p in pos]
        time = [t * 1e9 for t in time]
        plt.plot(time, pos, color=colors.hsv_to_rgb([(val - 0.1) / 10, 1, 1]))

    # fly time 1.507101420284057e-08
    voltFunction = lambda t: 5
    time, pos, _ = getData(voltFunction, timeLimit)
    pos = [p * 1e3 for p in pos]
    time = [t * 1e9 for t in time]
    plt.plot(time, pos, color=[0, 0, 0], label="постоянный ток")
    plt.plot([], [], color=None, label="угловая скорость w = 1 / время пролета при постоянном токе")
    plt.plot([], [], color=colors.hsv_to_rgb([0., 1., 1.]), label="w * 0,1")
    plt.plot([], [], color=colors.hsv_to_rgb([(5 - 0.1) / 10, 1., 1.]), label="w * 5")
    plt.legend()
    plt.title("Зависимость движения электрона при разной частоте переменного тока")
    plt.xlabel("Время с начала полета, ns")
    plt.ylabel("Расстояние от катода, мм")
    plt.show()

def binarySearchForOscillation():
    timeLimit = 8e-8
    def boolFunc(val):
        voltFunction = lambda t: 5 * math.cos(t * (1 / 1.507101420284057e-08) * val)
        time, pos, _ = getData(voltFunction, timeLimit)
        return len(time) < 5000
    result = binarySearch(boolFunc, 0.0000000001, 0.1, 5)
    print(result)
    # result 1.9997255 appr
    # theory is 4.186121254
    # due to fly time being different (they have multiplication by 3 instead of 2)

def speedCheck():
    timeLimit = 1e-7
    voltFunction = lambda t: 5
    time, pos, speed = getData(voltFunction, timeLimit)
    print(speed)
    # expected 1326205

if __name__ == '__main__':
    binarySearchForOscillation()
