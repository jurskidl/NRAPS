import matplotlib.pyplot as plt
import csv
import numpy as np

y1 = y2 = y3 = y4 = []

count = 1
for line in open('k_eff.csv', 'r'):
    k_eff = np.array(line.strip().split(',')).astype(float)

for line in open('interface.csv', 'r'):
    if count == 1:
        flux0 = np.array(line.strip().split(',')).astype(float)
    elif count == 2:
        flux1 = np.array(line.strip().split(',')).astype(float)
    else:
        fission = np.array(line.strip().split(',')).astype(float)

    count += 1

print(flux0[10])
print(flux1[10])

plt.plot(k_eff, label='multplication factor')
plt.title('Multiplication Factor')
# plt.show()
plt.savefig('./k_eff.svg', bbox_inches='tight')
plt.clf()

plt.plot(flux0, label='fast flux')
plt.title('Fast Flux')
# plt.show()
plt.savefig('./fast_flux.svg', bbox_inches='tight')
plt.clf()

plt.plot(flux1, label='thermal flux')
plt.title('Thermal Flux')
# plt.show()
plt.savefig('./thermal_flux.svg', bbox_inches='tight')
plt.clf()

plt.plot(fission, label='fission')
plt.title('Fission Density')
# plt.show()
plt.savefig('./fission.svg', bbox_inches='tight')
plt.clf()
