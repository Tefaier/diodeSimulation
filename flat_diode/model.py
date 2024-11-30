import math
from typing import List, Tuple

import numpy as np
import matplotlib.pyplot as plt
import scipy

ELECTRON_MASS = 9.109383719e-31
ELECTRON_CHARGE = -1.602176634e-19
ELECTRIC_CONSTANT = 8.85418781762039e-12 # E0
COULOMB_CONSTANT = 8.9875517923e9  # 1/4piE0
MAX_CHARGE_AT_CATHODE = -1e-10  # think about value


class FlatDiode1DModel:
    zero_potential: float = 0  # todo
    nets: List[Tuple[int, float]] = []

    def __init__(
            self,
            voltage: float,  # positive value makes charge flow
            cathode_radius: float,
            anode_radius: float,
            height: float,
            delta_time: float,
            amount_of_cells_between_cathode_and_anode: int,
            charge_increase_rate: float, # should be negative -> charge per second
            electron_leave_speed: float,
            appearance_part: float
    ):
        self.dt = delta_time
        self.n = amount_of_cells_between_cathode_and_anode
        self.r1 = cathode_radius
        self.r2 = anode_radius
        self.h = height
        self.big_charges = np.zeros(self.n)
        self.electric_field = np.zeros(self.n)
        self.velocities = np.zeros(self.n)
        self.increase_rate = charge_increase_rate
        self.start_speed = electron_leave_speed
        self.appearance_part = appearance_part
        self.set_voltage(voltage)

    # nets with first value being radius and second value being its charge
    def set_nets(self, nets: List[Tuple[float, float]]):
        withSlots = [(self.round_to_cell(net[0]), net[1]) for net in nets]
        withSlots = filter(lambda net: net[0] >= 0 and net[0] < self.n, withSlots)
        self.nets = sorted(withSlots, key=lambda net: net[0])

    def set_voltage(self, voltage: float):
        self.q_cathode = (-1 * (voltage + self.zero_potential) / (2 * COULOMB_CONSTANT) /
                          math.log(self.r2 / self.r1))

    def create_charge_by_poisson(self):
        n_max = max(1, int(self.n * self.appearance_part)) * 2
        return scipy.stats.poisson.pmf(np.arange(0, n_max), mu=n_max, loc=-n_max) * 2 * (self.increase_rate * self.dt)

    def calculate_pref_sum(self):
        scipy.stats.poisson.pmf(1, mu=10)
        # pref sum is n+1 size with 0 showing cathode!
        pref_sum = np.zeros(self.n + 1)
        pref_sum[0] = self.q_cathode
        net_pointer = 0
        for i in range(1, self.n + 1):
            pref_sum[i] = pref_sum[i-1] + self.big_charges[i-1]
            if net_pointer < len(self.nets) and self.nets[net_pointer][0] == i - 1:
                pref_sum[i] += self.nets[net_pointer][1]
                net_pointer += 1
        return pref_sum

    def iterate(self) -> float:
        electric_field = np.zeros(self.n)
        dr = (self.r2 - self.r1) / self.n  # delta r per cell
        start_velocity = max(self.start_speed, dr / self.dt)
        charges = self.create_charge_by_poisson()
        for i in range(0, len(charges)):
            new_q, new_v = self.compose([(self.big_charges[i], self.velocities[i]), (charges[i], start_velocity)])
            self.big_charges[i] = new_q
            self.velocities[i] = new_v
        # self.big_charges[0] += self.increase_rate * self.dt
        # if self.big_charges[0] < MAX_CHARGE_AT_CATHODE:  # because negative charge
        #     self.big_charges[0] = MAX_CHARGE_AT_CATHODE
        # self.velocities[0] = max(self.start_speed, dr / self.dt)
        pref_sum = self.calculate_pref_sum()
        new_cords = {}
        for i in range(0, self.n):
            if self.big_charges[i] != 0:
                r_i = (self.r1 + i * dr)
                electric_field[i] = (pref_sum[i] - pref_sum[-1] + pref_sum[i+1]) * (2 * COULOMB_CONSTANT) / (r_i * self.h)
                self.velocities[i] += electric_field[i] * ELECTRON_CHARGE / ELECTRON_MASS * self.dt
                new_i = self.round_to_cell(r_i + self.velocities[i] * self.dt)
                if new_i in new_cords:
                    new_cords[new_i].append((self.big_charges[i], self.velocities[i]))
                else:
                    new_cords[new_i] = [(self.big_charges[i], self.velocities[i])]
                self.big_charges[i] = 0
                self.velocities[i] = 0

        amperage = 0
        for new_i, arr in new_cords.items():
            new_q, new_v = self.compose(arr)
            if new_i >= self.n:
                amperage += new_q
            elif new_i < 0:
                continue
            else:
                self.big_charges[new_i] = new_q
                self.velocities[new_i] = new_v
        return amperage / self.dt

    def round_to_cell(self, x):
        if x <= self.r1:
            return -1  # fell on cathode
        if x >= self.r2:
            return self.n  # self.r2
        a = (x - self.r1) / (self.r2 - self.r1) * self.n
        if int(a + 0.5) == int(a):
            return int(a)  # self.r1 + int(a) * ((self.r2 - self.r1) / self.n)
        return int(a) + 1  # self.r1 + (int(a) + 1) * ((self.r2 - self.r1) / self.n)

    def compose(self, arr: List[Tuple[float, float]]) -> Tuple[float, float]:
        new_q = 0
        sum_imp = 0
        for q, v in arr:
            new_q += q
            sum_imp += v*q
        return new_q, sum_imp / new_q


delta_time = 1e-9
diode = FlatDiode1DModel(0,1e-3, 3e-3, 100e-3, delta_time, 10**3, ELECTRON_CHARGE * 1e15, 10000, 0.01)
x = []
y = []
for i in range(200):
    x.append(i * delta_time)
    y.append(diode.iterate())
window_size = 40
y2 = np.convolve(y, np.ones((window_size,)) / window_size, 'same')
plt.plot(x, y2)
plt.show()
