#!/usr/bin/env python
import numpy
import pylab

import csv
with open('size_params.csv','r') as file0:
	table = csv.reader(file0)
	for row in table:
		Mx = int(row[0])
		My = int(row[1])
file0.close()


x = numpy.zeros([Mx,My])
y = numpy.zeros([Mx,My])

p0 = numpy.zeros([Mx,My])
p1 = numpy.zeros([Mx,My])
p2 = numpy.zeros([Mx,My])
p3 = numpy.zeros([Mx,My])
p4 = numpy.zeros([Mx,My])
p5 = numpy.zeros([Mx,My])
p6 = numpy.zeros([Mx,My])
p7 = numpy.zeros([Mx,My])
p8 = numpy.zeros([Mx,My])
p9 = numpy.zeros([Mx,My])
p10 = numpy.zeros([Mx,My])

i = 0
j = 0
row_num = 0
import csv
with open('xU_Final.csv','r') as file1:
    table = csv.reader(file1)
    for row in table:
        p10[i, j] = row[12]

        j = (j+1) % (My)
        if (j == 0):
            i = (i+1)%(Mx)
     
file1.close()

i = 0
j = 0
row_num = 0
import csv
with open('xU_t0.csv','r') as file1:
    table = csv.reader(file1)
    for row in table:
	
        x[i , j] = row[0]
        y[i , j] = row[1]
        p0[i, j] = row[6]
        j = (j+1) % (My)
        if (j == 0):
            i = (i+1)%(Mx)
     
file1.close()

i = 0
j = 0
row_num = 0
import csv
with open('xU_t1.csv','r') as file1:
    table = csv.reader(file1)
    for row in table:
        p1[i, j] = row[6]
        j = (j+1) % (My)
        if (j == 0):
            i = (i+1)%(Mx)
     
file1.close()

i = 0
j = 0
row_num = 0
import csv
with open('xU_t2.csv','r') as file1:
    table = csv.reader(file1)
    for row in table:

        p2[i, j] = row[6]
        j = (j+1) % (My)
        if (j == 0):
            i = (i+1)%(Mx)
     
file1.close()

i = 0
j = 0
row_num = 0
import csv
with open('xU_t3.csv','r') as file1:
    table = csv.reader(file1)
    for row in table:
	
        p3[i, j] = row[6]
        j = (j+1) % (My)
        if (j == 0):
            i = (i+1)%(Mx)
     
file1.close()

i = 0
j = 0
row_num = 0
import csv
with open('xU_t4.csv','r') as file1:
    table = csv.reader(file1)
    for row in table:

        p4[i, j] = row[6]
        j = (j+1) % (My)
        if (j == 0):
            i = (i+1)%(Mx)
     
file1.close()

i = 0
j = 0
row_num = 0
import csv
with open('xU_t5.csv','r') as file1:
    table = csv.reader(file1)
    for row in table:

        p5[i, j] = row[6]
        j = (j+1) % (My)
        if (j == 0):
            i = (i+1)%(Mx)
     
file1.close()

i = 0
j = 0
row_num = 0
import csv
with open('xU_t6.csv','r') as file1:
    table = csv.reader(file1)
    for row in table:

        p6[i, j] = row[6]
        j = (j+1) % (My)
        if (j == 0):
            i = (i+1)%(Mx)
     
file1.close()

i = 0
j = 0
row_num = 0
import csv
with open('xU_t7.csv','r') as file1:
    table = csv.reader(file1)
    for row in table:

        p7[i, j] = row[6]
        j = (j+1) % (My)
        if (j == 0):
            i = (i+1)%(Mx)
     
file1.close()

i = 0
j = 0
row_num = 0
import csv
with open('xU_t8.csv','r') as file1:
    table = csv.reader(file1)
    for row in table:

        p8[i, j] = row[6]
        j = (j+1) % (My)
        if (j == 0):
            i = (i+1)%(Mx)
     
file1.close()

i = 0
j = 0
row_num = 0
import csv
with open('xU_t9.csv','r') as file1:
    table = csv.reader(file1)
    for row in table:

        p9[i, j] = row[6]
        j = (j+1) % (My)
        if (j == 0):
            i = (i+1)%(Mx)
     
file1.close()

fig100 = pylab.figure(100)
pylab.pcolor(pylab.array(x),pylab.array(y),pylab.array(p0),cmap = 'RdBu')
pylab.title('pressure initial')
pylab.colorbar()

fig1 = pylab.figure(1)
pylab.pcolor(pylab.array(x),pylab.array(y),pylab.array(p1),cmap = 'RdBu')
pylab.title('pressure t1')
pylab.colorbar()

fig2 = pylab.figure(2)
pylab.pcolor(pylab.array(x),pylab.array(y),pylab.array(p2),cmap = 'RdBu')
pylab.title('pressure t2')
pylab.colorbar()


fig3 = pylab.figure(3)
pylab.pcolor(pylab.array(x),pylab.array(y),pylab.array(p3),cmap = 'RdBu')
pylab.title('pressure t3')
pylab.colorbar()

fig4 = pylab.figure(4)
pylab.pcolor(pylab.array(x),pylab.array(y),pylab.array(p4),cmap = 'RdBu')
pylab.title('pressure t4')
pylab.colorbar()

fig5 = pylab.figure(5)
pylab.pcolor(pylab.array(x),pylab.array(y),pylab.array(p5),cmap = 'RdBu')
pylab.title('pressure t5')
pylab.colorbar()

fig6 = pylab.figure(6)
pylab.pcolor(pylab.array(x),pylab.array(y),pylab.array(p6),cmap = 'RdBu')
pylab.title('pressure t6')
pylab.colorbar()

fig7 = pylab.figure(7)
pylab.pcolor(pylab.array(x),pylab.array(y),pylab.array(p7),cmap = 'RdBu')
pylab.title('pressure t7')
pylab.colorbar()

fig8 = pylab.figure(8)
pylab.pcolor(pylab.array(x),pylab.array(y),pylab.array(p8),cmap = 'RdBu')
pylab.title('pressure t8')
pylab.colorbar()

fig9 = pylab.figure(9)
pylab.pcolor(pylab.array(x),pylab.array(y),pylab.array(p9),cmap = 'RdBu')
pylab.title('pressure t9')
pylab.colorbar()

fig10 = pylab.figure(10)
pylab.pcolor(pylab.array(x),pylab.array(y),pylab.array(p10),cmap = 'RdBu')
pylab.title('pressure t10')
pylab.colorbar()

pylab.show()


