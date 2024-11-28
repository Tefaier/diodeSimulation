from typing import List, Tuple

import numpy as np

ELECTRON_MASS = 9.109383719e-31
ELECTRON_CHARGE = -1.602176634e-19
COULOMB_CONSTANT = 8.9875517923e9  # 1/4piE0
MAX_CATHODE_DT_CHARGE = 1  # think about value


class FlatDiode1DModel:
    def __init__(self, voltage: float, cathode_radius: float, anode_radius: float, height: float, delta_time: float,
                 amount_of_cells_between_cathode_and_anode: int):
        zero_potential = 0  # todo
        self.q_cathode = (-1 * (voltage + zero_potential) / (2 * COULOMB_CONSTANT) /
                          np.log(anode_radius / cathode_radius))
        self.dt = delta_time
        self.n = amount_of_cells_between_cathode_and_anode
        self.r1 = cathode_radius
        self.r2 = anode_radius
        self.h = height
        self.big_charges = np.zeros(self.n)
        self.electric_field = np.zeros(self.n)
        self.velocities = np.zeros(self.n)

    def iterate(self) -> float:
        self.big_charges[0] += -1e-6 * self.dt
        electric_field = np.zeros(self.n)
        dr = (self.r2 - self.r1) / self.n  # delta r per cell
        self.velocities[0] = start_velocity = dr / self.dt
        if self.big_charges[0] > MAX_CATHODE_DT_CHARGE:
            self.big_charges[0] = MAX_CATHODE_DT_CHARGE
        pref_sum = 0
        new_cords = {}
        for i in range(0, self.n):
            if self.big_charges[i]:
                r_i = (self.r1 + i * dr)
                electric_field[i] = (self.q_cathode * 2 * COULOMB_CONSTANT / r_i / self.h +
                                     2 * COULOMB_CONSTANT / r_i / self.h * pref_sum)
                pref_sum += self.big_charges[i]
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
            else:
                self.big_charges[new_i] = new_q
                self.velocities[new_i] = new_v
        return amperage / self.dt

    def round_to_cell(self, x):
        assert (x > self.r1)
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


diode = FlatDiode1DModel(2e-2,1e-3, 3e-3, 100e-3, 1e-9, 10**4)
for i in range(200):
    print(diode.iterate())