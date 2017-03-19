import matplotlib.pyplot as plt

y = (
    0.009796129721,
    -0.01043321873,
    -0.05622882224,
    -0.08271575748,
    -0.08197990399,
    -0.1209313062,
    -0.1159941332,
    -0.1159432548,
    -0.1280633738,
    -0.1444161366,
    -0.1178353527,
    -0.1395782112,
    -0.1147615488,
    -0.11454177
)

x = (
    2,
    3,
    4,
    5,
    6,
    7,
    8,
    9,
    10,
    11,
    12,
    13,
    14,
    15,
)

large_pars_std = (
    0.01215827226,
    0.0174884351,
    0.06114779544,
    0.06218118886,
    0.04512208822,
    0.05012998331,
    0.05012998331,
    0.05867380453,
    0.03607619141,
    0.004052747056,
    0.05019767933,
    0.00246034457,
    0.04826998854,
    0.04991872104
)
import numpy as np

line =plt.plot(x, y, color='k')
plt.errorbar(x, y, large_pars_std, marker='^', capthick=1,elinewidth=1, ecolor='r')
plt.setp(line,color='k')
plt.xticks(np.arange(2, 15, 1))
plt.xlabel('Number of clusters')
plt.ylabel('Silhouette Score')
plt.title('Silhouette Score variation with number of clusters \nusing single-cell sequencing data')
plt.grid(True)
plt.show()
