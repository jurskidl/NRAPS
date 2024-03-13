import matplotlib
import csv

y1 = y2 = y3 = y4 = []

count = 1
for line in open('k_eff.csv', 'r'):
    k_eff = line
for line in open('interface.csv', 'r'):
    if count == 1:
        flux0 = line.strip().split(',')
    elif count == 2:
        flux1 = line.strip().split(',')
    else:
        fission = line.strip().split(',')

    count += 1

for item in k_eff:
    float(item)

print(k_eff)