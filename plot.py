import matplotlib.pyplot as plt
import csv
import numpy as np

count = 1
for line in open('vars.csv', 'r'):
    if count == 1:
        length = float(line.strip())
    elif count == 2:
        meshed = int(line.strip())
    elif count == 3:
        generations = int(line.strip())

    count += 1

k_fund = []
count = 1
for line in open('k_eff.csv', 'r'):
    if count == 1:
        k = np.array(line.strip().split(',')).astype(float)
    elif count == 2:
        k_fund = np.array(line.strip().split(',')).astype(float)

    count += 1

count = 1
for line in open('interface.csv', 'r'):
    if count == 1:
        flux0 = np.array(line.strip().split(',')).astype(float)
    elif count == 2:
        flux1 = np.array(line.strip().split(',')).astype(float)
    elif count == 3:
        average0 = np.array(line.strip().split(',')).astype(float)
    elif count == 4:
        average1 = np.array(line.strip().split(',')).astype(float)
    else:
        fission = np.array(line.strip().split(',')).astype(float)

    # if count == 1:
    #     flux0 = np.array(line.strip().split(',')).astype(float)
    # elif count == 2:
    #     flux1 = np.array(line.strip().split(',')).astype(float)
    # elif count == 3:
    #     flux2 = np.array(line.strip().split(',')).astype(float)
    # elif count == 4:
    #     flux3 = np.array(line.strip().split(',')).astype(float)
    # elif count == 5:
    #     average0 = np.array(line.strip().split(',')).astype(float)
    # elif count == 6:
    #     average1 = np.array(line.strip().split(',')).astype(float)
    # elif count == 7:
    #     average2 = np.array(line.strip().split(',')).astype(float)
    # elif count == 8:
    #     average3 = np.array(line.strip().split(',')).astype(float)
    # else:
    #     fission = np.array(line.strip().split(',')).astype(float)

    count += 1

x = np.linspace(0.0, length, num=meshed)

if type(k_fund) != float:
    plt.plot(k, 'o', markerfacecolor='none', label='multplication factor')
    plt.plot(k_fund, color='orange')
    plt.title('Multiplication Factor')
    plt.xlim(0, generations)
    plt.ylim(0, 2)

# plt.show()
plt.savefig('./k_eff.png', bbox_inches='tight')
plt.clf()

plt.plot(x, flux0, label='fast flux')
plt.plot(x, average0, color='orange')
plt.title('Fast Flux')
plt.xlim(0, length)
# plt.show()
plt.savefig('./fast_flux.png', bbox_inches='tight')
plt.clf()

plt.plot(x, flux1, label='thermal flux')
plt.plot(x, average1, color='orange')
plt.title('Thermal Flux')
plt.xlim(0, length)
# plt.show()
plt.savefig('./thermal_flux.png', bbox_inches='tight')
plt.clf()

# plt.plot(x, flux0, label='Group 1 flux')
# plt.plot(x, average0, color='orange')
# plt.title('Group 1 Flux')
# plt.xlim(0, length)
# # plt.show()
# plt.savefig('./group_1_flux.png', bbox_inches='tight')
# plt.clf()

# plt.plot(x, flux1, label='Group 2 flux')
# plt.plot(x, average1, color='orange')
# plt.title('Group 2 Flux')
# plt.xlim(0, length)
# # plt.show()
# plt.savefig('./group_2_flux.png', bbox_inches='tight')
# plt.clf()

# plt.plot(x, flux2, label='Group 3 flux')
# plt.plot(x, average2, color='orange')
# plt.title('Group 3 Flux')
# plt.xlim(0, length)
# # plt.show()
# plt.savefig('./group_3_flux.png', bbox_inches='tight')
# plt.clf()

# plt.plot(x, flux3, label='Group 4 flux')
# plt.plot(x, average3, color='orange')
# plt.title('Group 4 Flux')
# plt.xlim(0, length)
# # plt.show()
# plt.savefig('./group_4_flux.png', bbox_inches='tight')
# plt.clf()

# plt.plot(x, fission, label='fission')
# plt.title('Fission Density')
# plt.xlim(0, length)
# # plt.show()
# plt.savefig('./fission.png', bbox_inches='tight')
# plt.clf()
