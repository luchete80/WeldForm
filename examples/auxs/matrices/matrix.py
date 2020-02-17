import numpy as np

B = np.array([[1, 2, 3, 4], [5, 6, 7, 8]])

print(f'Original Array:\n{B}')

C = (B.transpose()).dot(B)
B_transpose = B.transpose()

print(f'Transposed Array:\n{B_transpose}')

print(f'BtB:\n{C}')