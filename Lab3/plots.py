# График для задания 3.3: табличная функция и приближающие многочлены МНК 1-й и 2-й степени.
# Запуск из корня репозитория: python Lab3/plots.py
# Картинка сохраняется в Lab3/images
import os

import matplotlib.pyplot as plt
import numpy as np

LAB_DIR = os.path.dirname(os.path.abspath(__file__))
DATA_FILE = os.path.join(LAB_DIR, "..", "data", "lab3", "input_3_3.txt")
IMAGES_DIR = os.path.join(LAB_DIR, "images")
os.makedirs(IMAGES_DIR, exist_ok=True)

# те же данные, что читает программа
with open(DATA_FILE) as f:
    lines = [line for line in f.read().splitlines() if line.strip()]
x = np.array([float(v) for v in lines[0].split()])
y = np.array([float(v) for v in lines[1].split()])


# нормальная система (3.17), как в least_squares.cpp
def least_squares(n):
    A = [[np.sum(x ** (k + i)) for i in range(n + 1)] for k in range(n + 1)]
    b = [np.sum(y * x**k) for k in range(n + 1)]
    return np.linalg.solve(A, b)


a1 = least_squares(1)
a2 = least_squares(2)
t = np.linspace(x.min() - 0.3, x.max() + 0.3, 300)

plt.figure(figsize=(8, 5))
plt.plot(x, y, "ko", label="табличные значения")
plt.plot(t, a1[0] + a1[1] * t, "b-", label=f"F1(x) = {a1[0]:.4f} + {a1[1]:.4f}x")
sign = "-" if a2[2] < 0 else "+"
plt.plot(t, a2[0] + a2[1] * t + a2[2] * t**2, "r--", label=f"F2(x) = {a2[0]:.4f} + {a2[1]:.4f}x {sign} {abs(a2[2]):.4f}x^2")
plt.xlabel("x")
plt.ylabel("y")
plt.title("Задание 3.3. Метод наименьших квадратов")
plt.grid(True)
plt.legend()
plt.savefig(os.path.join(IMAGES_DIR, "task_3_3.png"), dpi=100)
plt.close()

print("График сохранен в", IMAGES_DIR)
