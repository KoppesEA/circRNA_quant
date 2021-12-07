circRNAInNeg=[]
circRNAInPos=[]
inFile=open("./STARcustom/Control_1_S1/Control_1_S1Chimeric.out.junction", "r")
for aLine in inFile:
    fields=aLine.rstrip("\n").split("\t")
    if fields[6]<0:
        continue
    if fields[0]!=fields[3]:
        continue
    if fields[2]!=fields[5]:
        continue
    
    if (fields[2]=="-") and (int(fields[4])>int(fields[1])) and (int(fields[4])-int(fields[1])<1000000):
        circRNAInNeg.append((fields[9], aLine))
    elif (fields[2]=="+") and (int(fields[1])>int(fields[4])) and (int(fields[1])-int(fields[4])<1000000):
        circRNAInPos.append((fields[9], aLine))
inFile.close()

outFile=open("circRNAInNeg.csv", "w")
for name, aLine in circRNAInNeg:
    outFile.write(aLine)
outFile.close()
outFile=open("circRNAInPos.csv", "w")
for name, aLine in circRNAInPos:
    outFile.write(aLine)
outFile.close()
	