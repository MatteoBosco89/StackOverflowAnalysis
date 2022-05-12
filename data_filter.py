from datetime import datetime

s = open("sx-stackoverflow-a2q.txt", "r")
t = open("filtered-2016.txt", "w")
count = 0;
for line in s.readlines():
    count += 1
    timestamp = line.split()[2]
    src = line.split()[0]
    dst = line.split()[1]
    dt = datetime.fromtimestamp(int(timestamp))
    if(int(dt.year) > 2016): break
    if(int(dt.year) > 2015): t.write(line)
    if(count % 10000 == 0): print("Linea: " + str(count))