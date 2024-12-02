import matplotlib.pyplot as plt
import numpy as np

from flat_diode.model import FlatDiode1DModel, ELECTRON_CHARGE

delta_time = 1e-9
make_iterations = 300
calculate_mean_by = 50
x = []
y = []
for v in np.linspace(-5, 5, 40):
    print(f"{v} is in progress")
    diode = FlatDiode1DModel(v, 1e-3, 3e-3, 100e-3, delta_time, 10 ** 3, ELECTRON_CHARGE * 1e18, 10000, 0.01)
    sum_of_amperages = 0
    for i in range(0, make_iterations - calculate_mean_by):
        diode.iterate()
    for i in range(0, calculate_mean_by):
        sum_of_amperages += abs(diode.iterate())

    x.append(v)
    y.append(sum_of_amperages / calculate_mean_by)

plt.figure(figsize=(9, 5))
plt.title("Зависимость установившегося тока через диод при разных анодных напряжениях")
plt.xlabel("Напряжение, V")
plt.ylabel("Установившейся ток, A")
plt.plot(x, y)
plt.show()
