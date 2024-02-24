# coding: utf-8
import matplotlib
matplotlib.use('QtAgg')
import matplotlib.pyplot as plt

fig, ax = plt.subplots(2,2)

ax1 = plt.subplot(2,2,1)
ax1.plot(1,2)

plt.show()


