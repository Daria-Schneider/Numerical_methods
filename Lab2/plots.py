# Графики для лабораторной 2: по ним выбирается начальное приближение.
# Запуск из корня репозитория: python Lab2/plots.py
# Картинки сохраняются в Lab2/images
import os

import matplotlib.pyplot as plt
import numpy as np

IMAGES_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "images")
os.makedirs(IMAGES_DIR, exist_ok=True)


# 2.1: f(x) = x^6 - 5x^3 - 2
x = np.linspace(-1.2, 2.1, 1000)
y = x**6 - 5 * x**3 - 2

plt.figure(figsize=(8, 5))
plt.plot(x, y, label="f(x) = x^6 - 5x^3 - 2")
plt.axhline(0, color="black", linewidth=0.8)
plt.axvspan(1.5, 2, color="orange", alpha=0.25, label="отрезок [1.5, 2]")
plt.ylim(-12, 12)
plt.xlabel("x")
plt.ylabel("y")
plt.title("Задание 2.1")
plt.grid(True)
plt.legend()
plt.savefig(os.path.join(IMAGES_DIR, "task_2_1.png"), dpi=100)
plt.close()


# 2.2: a*x1^2 - x2 + x2^2 - a = 0 и x1 - sqrt(x2 + a) + 1 = 0
a = 3

# первое уравнение рисуем как линию уровня f1 = 0 (получается эллипс)
x1, x2 = np.meshgrid(np.linspace(-1.5, 1.6, 600), np.linspace(-3.2, 2.6, 600))
f1 = a * x1**2 - x2 + x2**2 - a

# из второго уравнения x1 выражается явно: x1 = sqrt(x2 + a) - 1
t = np.linspace(-a, 2.6, 600)
x1_curve = np.sqrt(t + a) - 1

plt.figure(figsize=(7, 7))
plt.contour(x1, x2, f1, levels=[0], colors="tab:blue")
plt.plot([], [], color="tab:blue", label="a*x1^2 - x2 + x2^2 - a = 0")
plt.plot(x1_curve, t, color="tab:red", label="x1 - sqrt(x2 + a) + 1 = 0")

# область G и начальное приближение, как в data/lab2/input_2_2.txt
plt.plot([0.95, 1.45, 1.45, 0.95, 0.95], [0.95, 0.95, 1.45, 1.45, 0.95], "k--", linewidth=1, label="область G")
plt.plot(1.2, 1.2, "ko", label="x0 = (1.2, 1.2)")

plt.axhline(0, color="black", linewidth=0.8)
plt.axvline(0, color="black", linewidth=0.8)
plt.xlabel("x1")
plt.ylabel("x2")
plt.title("Задание 2.2, a = 3")
plt.grid(True)
plt.legend(loc="lower right")
plt.savefig(os.path.join(IMAGES_DIR, "task_2_2.png"), dpi=100)
plt.close()

print("Графики сохранены в", IMAGES_DIR)
