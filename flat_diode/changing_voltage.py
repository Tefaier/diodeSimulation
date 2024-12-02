import math

from flat_diode.model import FlatDiode1DModel, ELECTRON_CHARGE
import numpy as np
import matplotlib.pyplot as plt

if __name__ == '__main__':
    plt.figure(figsize=(9, 5))
    plt.title("Анодный ток при переменном токе")
    plt.xlabel("Время, c")
    plt.ylabel("Ток, A")
    delta_time = 2e-10
    diode = FlatDiode1DModel(0, 1e-3, 3e-3, 100e-3, delta_time, 10 ** 3, ELECTRON_CHARGE * 1e18, 10000, 0.01)
    x = []
    y = []
    max_voltage = 3
    make_changes = 5
    iterations = 5000
    frequency = make_changes / (delta_time * iterations)
    v_function = lambda t: max_voltage * math.cos(t * frequency * 2 * math.pi)
    for i in range(iterations):
        diode.set_voltage(v_function(i * delta_time))
        x.append(i * delta_time)
        y.append(abs(diode.iterate()))
    window_size = 20
    y2 = np.convolve(y, np.ones((window_size,)) / window_size, 'same')
    plt.plot(x[:int(window_size / -2)], y2[:int(window_size / -2)])
    plt.show()