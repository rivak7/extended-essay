import math
import csv

g=9.80065
Cd=0.67
A=0.003425
p=1.229
T=73.4
mB=0.625
mC=0.610
tB=0.85

a0=(T)/(mB*g)-1
b0=(mB*g)/(0.5*Cd*A*p)
aC=-1
bC=(mC*g)/(0.5*Cd*A*p)

vB=0
a=0
v=0
s=0
t=0
dt=0.001

time=[]
altitude=[]
velocity=[]
acceleration=[]


while t<=tB:
    v=math.sqrt(b0*a0)*math.tanh(g*t*math.sqrt(a0/b0))
    velocity.append(v)
    s=(b0/g)*math.log(math.cosh(g*t*math.sqrt(a0/b0)))
    altitude.append(s)
    a=T/mB-g-(0.5*Cd*A*p*v*v/mB)
    acceleration.append(a)
    time.append(t)
    vB=v
    sB=s
    t+=dt

tC=(math.sqrt(bC)/g)*math.atan(vB/math.sqrt(bC))+tB
sC=(bC/(2*g))*math.log(1+math.pow(vB,2)/bC)

while t<=tC:
    v=math.sqrt(bC)*math.tan(math.atan(vB/(math.sqrt(bC)))-(g*(t-tB))/(math.sqrt(bC)))
    velocity.append(v)
    s=sC+(bC/g)*math.log(math.cos((g/math.sqrt(bC))*(tC-t)))+sB
    altitude.append(s)
    a=-g-(0.5*Cd*A*p*v*v/mC)
    acceleration.append(a)
    time.append(t)
    t+=dt

output_file = "analytical.csv"
with open(output_file, mode='w', newline='') as file:
    writer = csv.writer(file)
    writer.writerow(["Time (s)", "Altitude (m)", "Velocity (m/s)", "Acceleration (m/s^2)"])
    for i in range(len(time)):
        writer.writerow([time[i], altitude[i], velocity[i], acceleration[i]])
    
