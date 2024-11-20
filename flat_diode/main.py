from typing import Callable
from math import exp, sqrt

# универсальная термоэлектронная постоянная Зоммерфельда
UNIVERSAL_ZOMMERFELD_CONSTANT = 1200000  # Амп/m^2K^2
ELECTRON_CHARGE_CONSTANT = -1.602176634 * 10 ** -19
BOLTZMANN_CONSTANT = 1.38064852e-23


class FlatDiode:
    def __init__(self, square: float | None = None, distance: float | None = None,
                 perveance: float | None = None,
                 # work_function = Работа выхода материала катода
                 work_function: float | None = None,
                 # в кельвинах
                 temperature: float | None = None):
        assert distance or perveance
        assert square and work_function
        assert temperature > 0  # проверка на кельвины
        if perveance is not None:
            self.perveance = perveance
        else:
            self.perveance = 2.33 * (10 ** -6) * square / distance ** 2
        self.square = square
        self.work_function = work_function
        self.temperature = temperature
        self.saturated_amperage = square * UNIVERSAL_ZOMMERFELD_CONSTANT * temperature ** 2 * exp(
            - ELECTRON_CHARGE_CONSTANT * work_function / BOLTZMANN_CONSTANT / temperature
        )

    def GetAmperage(self, voltage: float):
        if voltage <= 0:
            return 0
        amperage_using_3_2_law = self.perveance * sqrt(voltage ** 3)
        if amperage_using_3_2_law > self.saturated_amperage:
            return self.saturated_amperage
        return amperage_using_3_2_law

    def GetAmperageFunc(self) -> Callable:
        return self.GetAmperage
