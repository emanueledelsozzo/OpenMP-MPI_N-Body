import sys

filename = sys.argv[1]

initGrid = "InitGrid"
initGridVal = []
genField = "GenerationField"
genFieldVal = []
partGen = "ParticleGeneration"
partGenVal = []
sysEv = "SystemEvolution"
sysEvVal = []

f = open(filename, "r")
for i in f.readlines():
	tok = i.split(" ")
	if initGrid == tok[0]:
		initGridVal.append(float(tok[-2]))
	elif genField == tok[0]:
		genFieldVal.append(float(tok[-2]))
	elif partGen == tok[0]:
		partGenVal.append(float(tok[-2]))
	elif sysEv == tok[0]:
		sysEvVal.append(float(tok[-2]))

print "initGrid list = %s" % (str(initGridVal))
print "genField list = %s" % (str(genFieldVal))
print "partGen list  = %s" % (str(partGenVal))
print "sysEv list    = %s" % (str(sysEvVal))

print "initGrid = %.10f" % (sum(initGridVal)/len(initGridVal))
print "genField = %.10f" % (sum(genFieldVal)/len(genFieldVal))
print "partGen  = %.10f" % (sum(partGenVal)/len(partGenVal))
print "sysEv    = %.10f" % (sum(sysEvVal)/len(sysEvVal))
