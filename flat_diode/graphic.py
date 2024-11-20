import matplotlib.pyplot as plt
import flat_diode.main as main


diode = main.FlatDiode(square=0.00005, distance=0.01, work_function=3, temperature=500)
x = []
y = []
for i in range(5, 50):
    x.append(i)
    y.append(diode.GetAmperage(i))
plt.plot(x, y)
plt.xlabel('Ось U, В')
plt.ylabel('Ось I, А')
plt.title('График зависимости анодного тока от анодного напряжения')
plt.show()
