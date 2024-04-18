import math
import random
Event_name = [random.randint(1, 99999) for _ in range(1000)]
Element = [14] * 1000
Energy = [10000] * 1000
X = [0]*1000
Y = [0]*1000
Z = [0]*1000
cosX = [1]*1000
cosY = [0]*1000
cosZ = [0]*1000
l = 1000
with open("TRIM.DAT", "w") as ofs:
    ofs.truncate()  # clear the file if it exists
    for j in range(10):
        ofs.write("#\r\n")
    for i in range(l):
        s = str(Event_name[i])
        if len(s) == 1:
            ofs.write(f"{Event_name[i]}      {Element[i]} {Energy[i]} {X[i]} {Y[i]} {Z[i]} {cosX[i]} {cosY[i]} {cosZ[i]}\r\n")
        elif len(s) == 2:
            ofs.write(f"{Event_name[i]}     {Element[i]} {Energy[i]} {X[i]} {Y[i]} {Z[i]} {cosX[i]} {cosY[i]} {cosZ[i]}\r\n")
        elif len(s) == 3:
            ofs.write(f"{Event_name[i]}    {Element[i]} {Energy[i]} {X[i]} {Y[i]} {Z[i]} {cosX[i]} {cosY[i]} {cosZ[i]}\r\n")
        elif len(s) == 4:
            ofs.write(f"{Event_name[i]}   {Element[i]} {Energy[i]} {X[i]} {Y[i]} {Z[i]} {cosX[i]} {cosY[i]} {cosZ[i]}\r\n")
        elif len(s) == 5:
            ofs.write(f"{Event_name[i]}  {Element[i]} {Energy[i]} {X[i]} {Y[i]} {Z[i]} {cosX[i]} {cosY[i]} {cosZ[i]}\r\n")
