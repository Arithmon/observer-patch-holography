"""Producer's exact Pauli arithmetic: i**phase X**x Z**z, little-endian bits."""

from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class Word:
    x: int = 0
    z: int = 0
    phase: int = 0

    def __mul__(self, other):
        return Word(self.x ^ other.x, self.z ^ other.z,
                    (self.phase + other.phase + 2*(self.z & other.x).bit_count()) % 4)

    def times_i(self, power):
        return Word(self.x, self.z, (self.phase + power) % 4)

    def commutes(self, other):
        return ((self.x & other.z).bit_count() + (self.z & other.x).bit_count()) % 2 == 0

    def row(self):
        return [str(self.x), str(self.z), self.phase]

    def weight(self):
        return (self.x | self.z).bit_count()

    def matrix(self, qubits):
        size = 1 << qubits
        out = np.zeros((size, size), complex)
        for j in range(size):
            out[j ^ self.x, j] = 1j**self.phase * (-1)**((j & self.z).bit_count())
        return out


def rotation(word, theta, qubits):
    """exp(i theta P) for a Hermitian Pauli involution."""
    if word*word != Word():
        raise ValueError('rotation requires a Hermitian involution')
    return np.cos(theta)*np.eye(1 << qubits) + 1j*np.sin(theta)*word.matrix(qubits)
