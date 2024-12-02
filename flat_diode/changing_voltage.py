import math

from matplotlib import colors

from flat_diode.model import FlatDiode1DModel, ELECTRON_CHARGE
import numpy as np
import matplotlib.pyplot as plt

if __name__ == '__main__':
    plt.figure(figsize=(9, 5))
    plt.title("Анодный ток при переменном токе")
    plt.xlabel("Время, c")
    plt.ylabel("Ток, A")
    delta_time = 1e-10
    max_voltage = 3
    make_changes = 500
    iterations = 5000
    for frequency in np.linspace(5e6, 4e7, 5):
        diode = FlatDiode1DModel(0, 1e-3, 3e-3, 100e-3, delta_time, 10 ** 4, ELECTRON_CHARGE * 1e18, 10000, 0.001)
        x = []
        y = []
        print(f"{frequency} is in progress")
        v_function = lambda t: max_voltage * math.cos(t * frequency * 2 * math.pi)
        for i in range(iterations):
            diode.set_voltage(v_function(i * delta_time))
            x.append(i * delta_time)
            y.append(abs(diode.iterate()))
        window_size = 200
        y2 = np.convolve(y, np.ones((window_size,)) / window_size, 'same')
        plt.plot(x[:int(window_size / -2)], y2[:int(window_size / -2)], color=colors.hsv_to_rgb([(frequency - 5e6) / (1.2*4e7), 1, 1]))
    plt.plot([], [], color=colors.hsv_to_rgb([0., 1., 1.]), label="5МГц")
    plt.plot([], [], color=colors.hsv_to_rgb([(4e7 - 5e6) / (1.2*4e7), 1., 1.]), label="40Мгц")
    plt.legend(loc="upper left")
    plt.show()